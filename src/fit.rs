// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev coefficient fitting via the DCT formula.
//!
//! Given function values at `N` Chebyshev nodes, computes the
//! Chebyshev coefficients `c_0, …, c_{N-1}` using the discrete
//! cosine transform identity:
//!
//! ```text
//! c_0 = (1/N) Σ_{k=0}^{N-1} f(ξ_k)
//! c_j = (2/N) Σ_{k=0}^{N-1} f(ξ_k) · cos(jπ(2k+1) / (2N))   for j ≥ 1
//! ```

use crate::scalar::ChebyScalar;
use crate::nodes;

/// Compute Chebyshev coefficients from function values at the
/// canonical Chebyshev nodes.
///
/// `values[k]` must correspond to the function evaluated at the `k`-th
/// node returned by [`nodes::<N>()`](crate::nodes).
///
/// # Example
///
/// ```
/// let xi: [f64; 9] = cheby::nodes();
/// let vals: [f64; 9] = std::array::from_fn(|k| xi[k].sin());
/// let coeffs = cheby::fit_coeffs(&vals);
/// // coeffs can now be passed to cheby::evaluate()
/// ```
#[inline]
pub fn fit_coeffs<T: ChebyScalar, const N: usize>(values: &[T; N]) -> [T; N] {
    let mut coeffs = [T::zero(); N];
    let n = N as f64;

    for j in 0..N {
        let mut sum = T::zero();
        for k in 0..N {
            let arg = std::f64::consts::PI * (j as f64) * (2.0 * k as f64 + 1.0) / (2.0 * n);
            sum = sum + values[k] * arg.cos();
        }
        coeffs[j] = if j == 0 { sum / n } else { sum * (2.0 / n) };
    }
    coeffs
}

/// Sample a function at `N` Chebyshev nodes on `[start, end]` and fit
/// Chebyshev coefficients.
///
/// This is a convenience wrapper combining [`nodes_mapped`](crate::nodes_mapped),
/// function sampling, and [`fit_coeffs`].
///
/// # Example
///
/// ```
/// let coeffs: [f64; 9] = cheby::fit_from_fn(|t| t.sin(), 0.0, 3.14);
/// let val = cheby::evaluate(&coeffs, 0.0); // evaluate at midpoint
/// ```
#[inline]
pub fn fit_from_fn<T: ChebyScalar, const N: usize>(
    f: impl Fn(f64) -> T,
    start: f64,
    end: f64,
) -> [T; N] {
    let mapped: [f64; N] = nodes::nodes_mapped(start, end);
    let values: [T; N] = std::array::from_fn(|k| f(mapped[k]));
    fit_coeffs(&values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evaluate;

    #[test]
    fn test_fit_constant() {
        // Fitting a constant function: all coeffs except c_0 should be ≈ 0
        let values = [5.0_f64; 9];
        let coeffs = fit_coeffs(&values);
        assert!((coeffs[0] - 5.0).abs() < 1e-14);
        for &c in &coeffs[1..] {
            assert!(c.abs() < 1e-14, "non-zero higher coeff: {c}");
        }
    }

    #[test]
    fn test_fit_linear() {
        // f(x) = x on [-1, 1] → c_0 = 0, c_1 = 1, rest = 0
        let xi: [f64; 9] = crate::nodes();
        let values: [f64; 9] = std::array::from_fn(|k| xi[k]);
        let coeffs = fit_coeffs(&values);
        assert!(coeffs[0].abs() < 1e-14, "c_0 = {}", coeffs[0]);
        assert!((coeffs[1] - 1.0).abs() < 1e-14, "c_1 = {}", coeffs[1]);
        for &c in &coeffs[2..] {
            assert!(c.abs() < 1e-14, "non-zero higher coeff: {c}");
        }
    }

    #[test]
    fn test_fit_sin_roundtrip() {
        // Fit sin(x) on [-1, 1] with degree 14 (15 coefficients)
        let coeffs: [f64; 15] = fit_from_fn(f64::sin, -1.0, 1.0);

        for &x in &[-0.9, -0.5, 0.0, 0.3, 0.8] {
            let approx = evaluate(&coeffs, x);
            let exact = x.sin();
            assert!(
                (approx - exact).abs() < 1e-12,
                "sin({x}): approx={approx}, exact={exact}"
            );
        }
    }

    #[test]
    fn test_fit_from_fn_mapped() {
        // Fit sin(t) on [0, π], evaluate at midpoint
        let coeffs: [f64; 15] = fit_from_fn(f64::sin, 0.0, std::f64::consts::PI);

        // Evaluate at t = π/2 → tau = 0.0 (midpoint maps to tau=0)
        let approx = evaluate(&coeffs, 0.0);
        assert!(
            (approx - 1.0).abs() < 1e-12,
            "sin(π/2) ≈ {approx}"
        );

        // Evaluate at t = π/4 → tau = -0.5
        let approx2 = evaluate(&coeffs, -0.5);
        let exact2 = (std::f64::consts::PI / 4.0).sin();
        assert!(
            (approx2 - exact2).abs() < 1e-10,
            "sin(π/4) ≈ {approx2}, exact = {exact2}"
        );
    }

    #[test]
    fn test_fit_quantity_type() {
        use qtty::Quantity;
        type Km = qtty::Kilometer;
        type Kilometers = Quantity<Km>;

        let xi: [f64; 9] = crate::nodes();
        let values: [Kilometers; 9] = std::array::from_fn(|k| Kilometers::new(xi[k].sin() * 1000.0));
        let coeffs = fit_coeffs(&values);
        // Evaluate and compare
        let val = evaluate(&coeffs, 0.0);
        let f64_vals: [f64; 9] = std::array::from_fn(|k| xi[k].sin() * 1000.0);
        let f64_coeffs = fit_coeffs(&f64_vals);
        let f64_val = evaluate(&f64_coeffs, 0.0);
        assert!(
            (val.value() - f64_val).abs() < 1e-10,
            "quantity val = {}, f64 val = {f64_val}",
            val.value()
        );
    }
}
