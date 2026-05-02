//! Scalar and unit-scaling traits.

use core::ops::{Add, Div, Mul, Sub};

/// A coefficient or value type usable in Chebyshev recurrences.
pub trait ChebyScalar:
    Copy + Add<Output = Self> + Sub<Output = Self> + Mul<f64, Output = Self> + Div<f64, Output = Self>
{
    /// Additive identity.
    fn zero() -> Self;

    /// Magnitude as `f64`, used for error estimates and convergence tests.
    fn magnitude(self) -> f64;
}

impl ChebyScalar for f64 {
    #[inline]
    fn zero() -> Self {
        0.0
    }

    #[inline]
    fn magnitude(self) -> f64 {
        self.abs()
    }
}

impl<U> ChebyScalar for qtty::Quantity<U>
where
    U: qtty::Unit,
{
    #[inline]
    fn zero() -> Self {
        Self::new(0.0)
    }

    #[inline]
    fn magnitude(self) -> f64 {
        self.value().abs()
    }
}

/// A coordinate type usable as a Chebyshev domain variable.
///
/// Same-unit subtraction and division normalize the coordinate to a raw
/// dimensionless `f64`, while multiplication by `f64` maps normalized nodes
/// back to the coordinate type.
pub trait ChebyTime:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<f64, Output = Self>
    + Div<Self, Output = f64>
    + PartialOrd
{
    /// Zero value in this coordinate type.
    fn zero() -> Self;
    /// Whether the stored scalar is finite.
    fn is_finite(self) -> bool;
}

impl ChebyTime for f64 {
    #[inline]
    fn zero() -> Self {
        0.0
    }

    #[inline]
    fn is_finite(self) -> bool {
        f64::is_finite(self)
    }
}

impl<U> ChebyTime for qtty::Quantity<U>
where
    U: qtty::Unit,
    qtty::Quantity<U>: Div<qtty::Quantity<U>, Output = f64>,
{
    #[inline]
    fn zero() -> Self {
        Self::new(0.0)
    }

    #[inline]
    fn is_finite(self) -> bool {
        self.value().is_finite()
    }
}

/// Scales a normalized derivative by a physical coordinate width.
///
/// Implemented blanketly for every pair where `Value / Coordinate` is a valid
/// Rust operation. With `qtty`, that makes derivatives dimensionally typed.
pub trait DifferentiateWith<X>: Sized {
    /// Dimensionally correct derivative type.
    type Derivative;

    /// Return `self / x`.
    fn scale_derivative(self, x: X) -> Self::Derivative;
}

impl<T, X> DifferentiateWith<X> for T
where
    T: Div<X>,
{
    type Derivative = <T as Div<X>>::Output;

    #[inline]
    fn scale_derivative(self, x: X) -> Self::Derivative {
        self / x
    }
}

/// Scales a normalized integral by a physical coordinate width.
///
/// Implemented blanketly for every pair where `Value * Coordinate` is valid.
/// With `qtty`, velocity integrated over time returns a length-like quantity.
pub trait IntegrateWith<X>: Sized {
    /// Dimensionally correct integral type.
    type Integral;

    /// Return `self * x`.
    fn scale_integral(self, x: X) -> Self::Integral;
}

impl<T, X> IntegrateWith<X> for T
where
    T: Mul<X>,
{
    type Integral = <T as Mul<X>>::Output;

    #[inline]
    fn scale_integral(self, x: X) -> Self::Integral {
        self * x
    }
}
