//! Coefficient-tail error estimates.

use crate::core::ChebyScalar;

/// Sum of magnitudes from `from` to the end of the coefficient slice.
#[inline]
pub fn tail_norm<T: ChebyScalar>(coeffs: &[T], from: usize) -> f64 {
    coeffs.iter().skip(from).map(|c| c.magnitude()).sum()
}

/// Estimate truncation error from the final quarter of coefficients.
#[inline]
pub fn estimated_truncation_error<T: ChebyScalar>(coeffs: &[T]) -> f64 {
    if coeffs.is_empty() {
        0.0
    } else {
        tail_norm(coeffs, (3 * coeffs.len()) / 4)
    }
}

/// Whether coefficients appear converged by absolute and relative tail tests.
#[inline]
pub fn is_converged<T: ChebyScalar>(coeffs: &[T], abs_tol: f64, rel_tol: f64) -> bool {
    let scale = coeffs
        .iter()
        .map(|c| c.magnitude())
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let err = estimated_truncation_error(coeffs);
    err <= abs_tol || err <= rel_tol * scale
}
