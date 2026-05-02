//! Clenshaw evaluation for Chebyshev series.

use super::ChebyScalar;

/// Evaluate a Chebyshev series at normalized coordinate `x`.
#[inline]
pub fn evaluate<T: ChebyScalar>(coeffs: &[T], x: f64) -> T {
    let n = coeffs.len();
    if n == 0 {
        return T::zero();
    }
    if n == 1 {
        return coeffs[0];
    }

    let two_x = 2.0 * x;
    let mut b1 = T::zero();
    let mut b2 = T::zero();
    for k in (1..n).rev() {
        let b0 = b1 * two_x - b2 + coeffs[k];
        b2 = b1;
        b1 = b0;
    }
    coeffs[0] + b1 * x - b2
}

/// Evaluate `df/dx` for a normalized-coordinate Chebyshev series.
#[inline]
pub fn evaluate_derivative<T: ChebyScalar>(coeffs: &[T], x: f64) -> T {
    evaluate_both(coeffs, x).1
}

/// Evaluate both value and derivative with respect to normalized `x`.
#[inline]
pub fn evaluate_both<T: ChebyScalar>(coeffs: &[T], x: f64) -> (T, T) {
    let n = coeffs.len();
    if n == 0 {
        return (T::zero(), T::zero());
    }
    if n == 1 {
        return (coeffs[0], T::zero());
    }

    let two_x = 2.0 * x;
    let mut b1 = T::zero();
    let mut b2 = T::zero();
    let mut db1 = T::zero();
    let mut db2 = T::zero();

    for k in (1..n).rev() {
        let b0 = b1 * two_x - b2 + coeffs[k];
        let db0 = db1 * two_x - db2 + b1 * 2.0;
        b2 = b1;
        b1 = b0;
        db2 = db1;
        db1 = db0;
    }

    let value = coeffs[0] + b1 * x - b2;
    let derivative = b1 + db1 * x - db2;
    (value, derivative)
}
