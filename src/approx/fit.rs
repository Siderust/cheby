//! Coefficient fitting from samples and functions.
//!
//! The fixed-size fitter samples at roots of `T_N` and applies the direct
//! discrete-cosine coefficient formula. With `T_j(x_k) = cos(j θ_k)` evaluated
//! by Chebyshev recurrence in `j`, the work is `O(N²)` multiplies-and-adds
//! and only `O(N)` calls to `cos` (used to build the node grid once).

use crate::core::{nodes, ChebyError, ChebyScalar, ChebySeries, ChebySeriesOn, Domain, NodeKind};

#[cfg(feature = "alloc")]
use crate::core::ChebySeriesDyn;

/// Compute Chebyshev coefficients from values sampled at roots of `T_N`.
///
/// Uses the Chebyshev three-term recurrence to generate `T_j(x_k)` for each
/// sample `k` in a single pass over `j`, avoiding `N²` `cos()` calls.
#[inline]
pub fn fit_coeffs<T: ChebyScalar, const N: usize>(values: &[T; N]) -> [T; N] {
    let mut coeffs = [T::zero(); N];
    if N == 0 {
        return coeffs;
    }
    let xs = nodes::<N>(NodeKind::Roots);
    for (k, &v) in values.iter().enumerate() {
        let x = xs[k];
        // T_0 = 1
        coeffs[0] = coeffs[0] + v;
        if N > 1 {
            // T_1 = x
            coeffs[1] = coeffs[1] + v * x;
            let two_x = 2.0 * x;
            let mut tkm1 = 1.0_f64;
            let mut tk = x;
            for coeff in coeffs.iter_mut().take(N).skip(2) {
                let tkp1 = two_x * tk - tkm1;
                *coeff = *coeff + v * tkp1;
                tkm1 = tk;
                tk = tkp1;
            }
        }
    }
    let nf = N as f64;
    coeffs[0] = coeffs[0] / nf;
    let scale = 2.0 / nf;
    for coeff in coeffs.iter_mut().skip(1) {
        *coeff = *coeff * scale;
    }
    coeffs
}

/// Fit a Chebyshev series on normalized `[-1, 1]`.
#[inline]
pub fn fit_from_fn<T: ChebyScalar, const N: usize>(f: impl Fn(f64) -> T) -> ChebySeries<T, N> {
    let xs = nodes::<N>(NodeKind::Roots);
    let values = core::array::from_fn(|k| f(xs[k]));
    ChebySeries::new(fit_coeffs(&values))
}

/// Fit a dynamic Chebyshev series of `degree` on normalized `[-1, 1]`.
#[cfg(feature = "alloc")]
pub fn fit_dyn_from_fn(
    degree: usize,
    mut f: impl FnMut(f64) -> f64,
) -> Result<ChebySeriesDyn<f64>, ChebyError> {
    let n = degree.saturating_add(1).max(2);
    let nf = n as f64;
    let mut xs = alloc::vec![0.0_f64; n];
    let mut values = alloc::vec![0.0_f64; n];
    for (k, x_slot) in xs.iter_mut().enumerate() {
        let x = (core::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * nf)).cos();
        *x_slot = x;
        values[k] = f(x);
    }

    let mut coeffs = alloc::vec![0.0_f64; n];
    for (k, &v) in values.iter().enumerate() {
        let x = xs[k];
        coeffs[0] += v;
        if n > 1 {
            coeffs[1] += v * x;
            let two_x = 2.0 * x;
            let mut tkm1 = 1.0_f64;
            let mut tk = x;
            for coeff in coeffs.iter_mut().take(n).skip(2) {
                let tkp1 = two_x * tk - tkm1;
                *coeff += v * tkp1;
                tkm1 = tk;
                tk = tkp1;
            }
        }
    }
    coeffs[0] /= nf;
    let scale = 2.0 / nf;
    for coeff in coeffs.iter_mut().skip(1) {
        *coeff *= scale;
    }

    ChebySeriesDyn::new(coeffs)
}

/// Fit a Chebyshev series on a physical domain.
#[inline]
pub fn fit_from_fn_on<T, X, const N: usize>(
    domain: Domain<X>,
    f: impl Fn(X) -> T,
) -> ChebySeriesOn<T, X, N>
where
    T: ChebyScalar,
    X: crate::core::ChebyTime,
{
    let xs = crate::core::nodes::nodes_mapped::<X, N>(domain, NodeKind::Roots);
    let values = core::array::from_fn(|k| f(xs[k]));
    ChebySeriesOn::new(domain, ChebySeries::new(fit_coeffs(&values)))
}
