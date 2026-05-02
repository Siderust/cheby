//! Coefficient fitting from samples and functions.

use crate::core::{nodes, ChebyScalar, ChebySeries, ChebySeriesOn, Domain, NodeKind};

/// Compute Chebyshev coefficients from values sampled at roots of `T_N`.
#[inline]
pub fn fit_coeffs<T: ChebyScalar, const N: usize>(values: &[T; N]) -> [T; N] {
    let mut coeffs = [T::zero(); N];
    let n = N as f64;
    for (j, coeff) in coeffs.iter_mut().enumerate() {
        let mut sum = T::zero();
        for (k, value) in values.iter().enumerate() {
            let arg = core::f64::consts::PI * j as f64 * (2.0 * k as f64 + 1.0) / (2.0 * n);
            sum = sum + *value * arg.cos();
        }
        *coeff = if j == 0 { sum / n } else { sum * (2.0 / n) };
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
