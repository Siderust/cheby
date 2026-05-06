//! Coefficient-tail error estimates.

use crate::core::ChebyScalar;

/// Sum of magnitudes from `from` to the end of the coefficient slice.
#[inline]
pub fn tail_norm<T: ChebyScalar>(coeffs: &[T], from: usize) -> f64 {
    coeffs.iter().skip(from).map(|c| c.magnitude()).sum()
}

/// Estimate truncation error from the final quarter of coefficients.
///
/// For smooth analytic targets the magnitude of the dropped tail bounds the
/// truncation error well. For localized features or endpoint stress this is
/// a heuristic, not a guaranteed bound; see [`is_converged`] for a stricter
/// test that also checks the trailing few coefficients individually.
#[inline]
pub fn estimated_truncation_error<T: ChebyScalar>(coeffs: &[T]) -> f64 {
    if coeffs.is_empty() {
        0.0
    } else {
        tail_norm(coeffs, (3 * coeffs.len()) / 4)
    }
}

/// Whether coefficients appear converged.
///
/// Combines the quarter-tail estimate with an explicit check on the last
/// few coefficient magnitudes so that a single oscillating high-order entry
/// cannot be hidden by averaging across the tail.
#[inline]
pub fn is_converged<T: ChebyScalar>(coeffs: &[T], abs_tol: f64, rel_tol: f64) -> bool {
    let scale = coeffs
        .iter()
        .map(|c| c.magnitude())
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let tail_err = estimated_truncation_error(coeffs);
    // Largest magnitude among the last min(3, len) coefficients.
    let trailing = coeffs
        .iter()
        .rev()
        .take(3)
        .map(|c| c.magnitude())
        .fold(0.0_f64, f64::max);
    let err = tail_err.max(trailing);
    err <= abs_tol || err <= rel_tol * scale
}
