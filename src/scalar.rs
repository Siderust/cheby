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

// ── ChebyTime — typed time for segment domain parameters ─────────────────

/// A type usable as the time parameter in Chebyshev segment evaluation.
///
/// Implemented for `f64` (raw seconds or JD) and for any `qtty::Quantity<U>`
/// typed time value.  Using a typed time type removes the need to call
/// `.value()` at call sites.
///
/// The trait requires only the arithmetic needed for τ normalisation:
/// `τ = (t − mid) / half_width`.  For `Quantity<U>`, this division naturally
/// cancels the unit and returns a dimensionless `f64` via qtty's `UnitDiv`
/// machinery — no extraction required.
pub trait ChebyTime:
    Copy
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<f64, Output = Self>
    + std::ops::Div<Self, Output = f64>
    + PartialOrd
    + Sized
{
    /// The additive identity (zero in this time domain).
    fn zero() -> Self;

    /// Returns `1.0 / self` as a raw `f64` scaling factor.
    ///
    /// Used at segment construction time to precompute the derivative
    /// chain-rule factor `dτ/dt = 1/half`. The hot evaluation path uses
    /// the stored `f64`.
    fn recip_f64(self) -> f64;
}

impl ChebyTime for f64 {
    #[inline]
    fn zero() -> Self {
        0.0
    }

    #[inline]
    fn recip_f64(self) -> f64 {
        1.0 / self
    }
}

impl<U> ChebyTime for qtty::Quantity<U>
where
    U: qtty::Unit,
    qtty::Quantity<U>: std::ops::Div<qtty::Quantity<U>, Output = f64>,
{
    #[inline]
    fn zero() -> Self {
        Self::new(0.0)
    }

    #[inline]
    fn recip_f64(self) -> f64 {
        1.0 / self.value()
    }
}
