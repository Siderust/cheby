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
//!
//! # Convergence and how many coefficients to keep
//!
//! For a function `f` analytic on `[-1, 1]` (and analytically continuable
//! to a Bernstein ellipse `E_ρ`, foci `±1`, semi-axis sum `ρ > 1`) the
//! Chebyshev coefficients decay geometrically:
//!
//! ```text
//! |a_n| ≤ M · ρ^{-n}.
//! ```
//!
//! Once the tail is well into the geometric regime the truncation error is
//! dominated by the first dropped coefficient,
//! `|f(τ) − Σ_{k<N} a_k T_k(τ)| ≈ |a_N|`. A practical recipe:
//!
//! 1. Fit with `N` somewhat larger than you expect to need.
//! 2. Look at the trailing coefficients; once they are at the level
//!    `≈ ε · max_k |a_k|` with `ε` your target tolerance, you can either
//!    keep the full set (cheap to evaluate via Clenshaw) or truncate.
//!
//! For non-analytic but smooth `f`, decay is algebraic rather than
//! geometric and `N` may need to be much larger.
//!
//! # Domain rescaling
//!
//! [`fit_from_fn`] and [`fit_from_fn_t`] sample on a physical interval
//! `[start, end]` after applying the standard affine map
//! `τ = (2x − (start + end)) / (end − start)`. The user is responsible for
//! later evaluating only at points inside `[start, end]` — outside that
//! interval the Chebyshev basis polynomials grow like
//! `cosh(k · acosh|τ|)` and extrapolation diverges exponentially in the
//! degree.

use crate::nodes;
use crate::scalar::ChebyScalar;

/// Compute Chebyshev coefficients from function values at the
/// canonical Chebyshev nodes.
///
/// `values[k]` must correspond to the function evaluated at the `k`-th
/// node returned by [`nodes::<N>()`](crate::nodes).
///
/// # Complexity
///
/// O(N²) multiplications. For the typical degrees used in ephemeris work
/// (N ≤ 30) this is negligible. For significantly larger N an FFT-based
/// DCT (O(N log N)) would be needed.
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

    for (j, coeff) in coeffs.iter_mut().enumerate() {
        let mut sum = T::zero();
        for (k, value) in values.iter().enumerate() {
            let arg = std::f64::consts::PI * (j as f64) * (2.0 * k as f64 + 1.0) / (2.0 * n);
            sum = sum + *value * arg.cos();
        }
        *coeff = if j == 0 { sum / n } else { sum * (2.0 / n) };
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
    fit_from_fn_t(&f, start, end)
}

/// Sample a function at `N` Chebyshev nodes on `[start, end]` and fit
/// Chebyshev coefficients. Generic over the time domain type `Tt`.
///
/// This variant accepts typed time values (e.g., `Quantity<Second>`),
/// removing the need to call `.value()` at call sites.
#[inline]
pub fn fit_from_fn_t<T: ChebyScalar, Tt: crate::scalar::ChebyTime, const N: usize>(
    f: &impl Fn(Tt) -> T,
    start: Tt,
    end: Tt,
) -> [T; N] {
    let mapped: [Tt; N] = nodes::nodes_mapped_t(start, end);
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
        assert!((approx - 1.0).abs() < 1e-12, "sin(π/2) ≈ {approx}");

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
        use qtty::Kilometer;

        let xi: [f64; 9] = crate::nodes();
        let values: [Kilometer; 9] = std::array::from_fn(|k| Kilometer::new(xi[k].sin() * 1000.0));
        let coeffs = fit_coeffs(&values);
        let val = evaluate(&coeffs, 0.0);
        let f64_vals: [f64; 9] = std::array::from_fn(|k| xi[k].sin() * 1000.0);
        let f64_coeffs = fit_coeffs(&f64_vals);
        let f64_val = evaluate(&f64_coeffs, 0.0);
        assert!(
            (val - Kilometer::new(f64_val)).abs() < Kilometer::new(1e-10),
            "quantity val = {val}, f64 val = {f64_val}"
        );
    }

    /// Validates the practical truncation-error bound for analytic
    /// functions: the error of the degree-`(N-1)` Chebyshev approximation
    /// is dominated by the magnitude of the *first dropped* coefficient.
    ///
    /// We fit `cos(x)` on `[-1, 1]` with `N` coefficients, then refit with
    /// `N + 1` coefficients to obtain `a_N` (the first one we would have
    /// dropped). The max-abs error on a fine grid should not exceed
    /// `3 · |a_N|`. The factor of 3 comfortably covers the geometric tail
    /// `Σ_{k≥N} |a_k| ≈ |a_N| / (1 − 1/ρ)` plus accumulated round-off.
    #[test]
    fn test_truncation_error_bound_cos() {
        const N: usize = 8;
        let coeffs_n: [f64; N] = fit_from_fn(f64::cos, -1.0, 1.0);
        let coeffs_np1: [f64; N + 1] = fit_from_fn(f64::cos, -1.0, 1.0);
        let a_next = coeffs_np1[N].abs();
        assert!(a_next > 0.0, "expected a non-trivial first dropped coeff");

        let mut max_err: f64 = 0.0;
        let samples = 2001;
        for i in 0..samples {
            let x = -1.0 + 2.0 * (i as f64) / ((samples - 1) as f64);
            let approx = evaluate(&coeffs_n, x);
            let err = (approx - x.cos()).abs();
            if err > max_err {
                max_err = err;
            }
        }

        let bound = 3.0 * a_next;
        assert!(
            max_err <= bound,
            "max error {max_err:e} exceeded 3·|a_N| = {bound:e} (|a_N| = {a_next:e})"
        );
    }

    /// Sanity check: a Chebyshev series fitted from samples at the
    /// canonical Chebyshev nodes must reproduce those samples *exactly*
    /// (to round-off) at the same nodes — this is the interpolation
    /// property of the DCT-based fit.
    #[test]
    fn test_reconstruction_at_nodes() {
        const N: usize = 16;
        let xi: [f64; N] = crate::nodes();
        let f = |x: f64| (2.0 * x).cos() + 0.5 * x.sin();
        let values: [f64; N] = std::array::from_fn(|k| f(xi[k]));
        let coeffs = fit_coeffs(&values);

        for k in 0..N {
            let reconstructed = evaluate(&coeffs, xi[k]);
            let diff = (reconstructed - values[k]).abs();
            assert!(
                diff < 1e-13,
                "node {k}: reconstructed {reconstructed} vs sample {} (diff = {diff:e})",
                values[k]
            );
        }
    }
}
