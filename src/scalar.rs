// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Scalar trait for Chebyshev operations.
//!
//! [`ChebyScalar`] abstracts over numeric types that can participate in
//! Chebyshev evaluation and fitting. Implemented for `f64` and for all
//! `qtty::Quantity<U>` types.

use std::ops::{Add, Div, Mul, Sub};

/// A scalar type usable as a Chebyshev coefficient or value.
///
/// This trait requires the arithmetic needed by the Clenshaw recurrence
/// and the DCT coefficient computation:
///
/// - Addition and subtraction of two values of the same type.
/// - Multiplication and division by a dimensionless `f64`.
/// - A zero element.
pub trait ChebyScalar:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
    + Sized
{
    /// The additive identity (zero).
    fn zero() -> Self;
}

// ── f64 implementation ──────────────────────────────────────────────────

impl ChebyScalar for f64 {
    #[inline]
    fn zero() -> Self {
        0.0
    }
}

// ── qtty::Quantity blanket implementation ────────────────────────────────

impl<U> ChebyScalar for qtty::Quantity<U>
where
    U: qtty::Unit,
{
    #[inline]
    fn zero() -> Self {
        Self::new(0.0)
    }
}
