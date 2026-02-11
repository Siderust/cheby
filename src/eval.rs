// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev polynomial evaluation using the Clenshaw recurrence.
//!
//! Given coefficients `c[0..n]` and normalised argument `tau ∈ [-1, 1]`,
//! evaluates:
//!
//! ```text
//! f(tau) = c[0]*T_0(tau) + c[1]*T_1(tau) + ... + c[n-1]*T_{n-1}(tau)
//! ```
//!
//! where `T_k` are Chebyshev polynomials of the first kind.
//!
//! All functions are generic over [`ChebyScalar`](crate::ChebyScalar).

use crate::ChebyScalar;

/// Evaluate a Chebyshev polynomial using the Clenshaw algorithm.
///
/// - `coeffs`: Chebyshev coefficients `c[0..n]`.
/// - `tau`: normalised time in `[-1, 1]`.
///
/// Returns the polynomial value.
#[inline]
pub fn evaluate<T: ChebyScalar>(coeffs: &[T], tau: f64) -> T {
    let n = coeffs.len();
    if n == 0 {
        return T::zero();
    }
    if n == 1 {
        return coeffs[0];
    }

    // Clenshaw recurrence (backwards):
    //   b_{n+1} = 0
    //   b_n     = 0
    //   b_k     = 2*tau*b_{k+1} - b_{k+2} + c[k]   for k = n-1, ..., 1
    //   result  = c[0] + tau*b_1 - b_2
    let two_tau = 2.0 * tau;
    let mut b_kp1 = T::zero(); // b_{k+1}
    let mut b_kp2 = T::zero(); // b_{k+2}

    for k in (1..n).rev() {
        let b_k = b_kp1 * two_tau - b_kp2 + coeffs[k];
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }

    coeffs[0] + b_kp1 * tau - b_kp2
}

/// Evaluate the derivative of a Chebyshev polynomial.
///
/// Returns `df/dtau` — the derivative with respect to the normalised
/// argument `tau`. To get `df/dt` in physical time units, multiply by
/// `dtau/dt = 1/half_width`.
///
/// Uses the modified Clenshaw recurrence that tracks both value and
/// derivative simultaneously.
#[inline]
pub fn evaluate_derivative<T: ChebyScalar>(coeffs: &[T], tau: f64) -> T {
    let n = coeffs.len();
    if n <= 1 {
        return T::zero();
    }

    let two_tau = 2.0 * tau;
    let mut b_kp1 = T::zero();
    let mut b_kp2 = T::zero();
    let mut db_kp1 = T::zero();
    let mut db_kp2 = T::zero();

    for k in (1..n).rev() {
        let b_k = b_kp1 * two_tau - b_kp2 + coeffs[k];
        let db_k = db_kp1 * two_tau - db_kp2 + b_kp1 * 2.0;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
        db_kp2 = db_kp1;
        db_kp1 = db_k;
    }

    // derivative of (c[0] + tau*b_1 - b_2) w.r.t. tau
    // = b_1 + tau*db_1 - db_2
    b_kp1 + db_kp1 * tau - db_kp2
}

/// Evaluate both the Chebyshev polynomial and its derivative in one pass.
///
/// Returns `(value, d_value_d_tau)`.
#[inline]
pub fn evaluate_both<T: ChebyScalar>(coeffs: &[T], tau: f64) -> (T, T) {
    let n = coeffs.len();
    if n == 0 {
        return (T::zero(), T::zero());
    }
    if n == 1 {
        return (coeffs[0], T::zero());
    }

    let two_tau = 2.0 * tau;
    let mut b_kp1 = T::zero();
    let mut b_kp2 = T::zero();
    let mut db_kp1 = T::zero();
    let mut db_kp2 = T::zero();

    for k in (1..n).rev() {
        let b_k = b_kp1 * two_tau - b_kp2 + coeffs[k];
        let db_k = db_kp1 * two_tau - db_kp2 + b_kp1 * 2.0;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
        db_kp2 = db_kp1;
        db_kp1 = db_k;
    }

    let value = coeffs[0] + b_kp1 * tau - b_kp2;
    let deriv = b_kp1 + db_kp1 * tau - db_kp2;
    (value, deriv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        assert_eq!(evaluate::<f64>(&[], 0.5), 0.0);
        assert_eq!(evaluate_derivative::<f64>(&[], 0.5), 0.0);
        let (v, d) = evaluate_both::<f64>(&[], 0.5);
        assert_eq!(v, 0.0);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_constant() {
        // T_0(x) = 1, so [c] evaluates to c for any tau
        assert!((evaluate(&[3.5], 0.0) - 3.5).abs() < 1e-15);
        assert!((evaluate(&[3.5], 0.7) - 3.5).abs() < 1e-15);
        assert!(evaluate_derivative(&[3.5], 0.7).abs() < 1e-15);
    }

    #[test]
    fn test_linear() {
        // c[0]*T_0 + c[1]*T_1 = 2 + 3*tau
        let coeffs = [2.0, 3.0];
        assert!((evaluate(&coeffs, 0.0) - 2.0).abs() < 1e-15);
        assert!((evaluate(&coeffs, 1.0) - 5.0).abs() < 1e-15);
        assert!((evaluate(&coeffs, -1.0) - (-1.0)).abs() < 1e-15);
        assert!((evaluate_derivative(&coeffs, 0.0) - 3.0).abs() < 1e-15);
        assert!((evaluate_derivative(&coeffs, 1.0) - 3.0).abs() < 1e-15);
    }

    #[test]
    fn test_quadratic() {
        // c[0]*T_0 + c[1]*T_1 + c[2]*T_2 = 1 + 0*tau + 2*(2*tau^2 - 1)
        // = 1 + 4*tau^2 - 2 = -1 + 4*tau^2
        let coeffs = [1.0, 0.0, 2.0];
        assert!((evaluate(&coeffs, 0.0) - (-1.0)).abs() < 1e-14);
        assert!((evaluate(&coeffs, 1.0) - 3.0).abs() < 1e-14);
        // derivative: 8*tau
        assert!((evaluate_derivative(&coeffs, 0.5) - 4.0).abs() < 1e-14);
    }

    #[test]
    fn test_evaluate_both_matches() {
        let coeffs = [1.0, 2.0, 3.0, 4.0, 5.0];
        let tau = 0.37;
        let (val, deriv) = evaluate_both(&coeffs, tau);
        assert!((val - evaluate(&coeffs, tau)).abs() < 1e-14);
        assert!((deriv - evaluate_derivative(&coeffs, tau)).abs() < 1e-14);
    }

    #[test]
    fn test_quantity_type() {
        use qtty::Quantity;
        type Km = qtty::Kilometer;
        type Kilometers = Quantity<Km>;

        let coeffs: [Kilometers; 3] = [
            Kilometers::new(100.0),
            Kilometers::new(50.0),
            Kilometers::new(10.0),
        ];
        let tau = 0.5;
        let val = evaluate(&coeffs, tau);
        // Compare with f64 evaluation
        let f64_coeffs = [100.0, 50.0, 10.0];
        let f64_val = evaluate(&f64_coeffs, tau);
        assert!((val.value() - f64_val).abs() < 1e-12);
    }
}
