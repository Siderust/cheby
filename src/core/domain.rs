//! First-class Chebyshev domains.

use super::{ChebyError, ChebyTime};

/// A closed interval mapped affinely to the Chebyshev coordinate `[-1, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Domain<X> {
    start: X,
    end: X,
    mid: X,
    half: X,
}

impl<X: ChebyTime> Domain<X> {
    /// Construct a domain, returning a typed error for invalid inputs.
    #[inline]
    pub fn try_new(start: X, end: X) -> Result<Self, ChebyError> {
        if !start.is_finite() || !end.is_finite() {
            return Err(ChebyError::NonFiniteInput);
        }
        let width = end - start;
        if width == X::zero() {
            return Err(ChebyError::EmptyDomain);
        }
        if width < X::zero() {
            return Err(ChebyError::InvalidDomain);
        }
        let half = width * 0.5;
        let mid = start + half;
        Ok(Self {
            start,
            end,
            mid,
            half,
        })
    }

    /// Construct a domain and panic if it is invalid.
    #[inline]
    pub fn new(start: X, end: X) -> Self {
        Self::try_new(start, end).expect("invalid Chebyshev domain")
    }

    /// Domain start.
    #[inline]
    pub fn start(&self) -> X {
        self.start
    }

    /// Domain end.
    #[inline]
    pub fn end(&self) -> X {
        self.end
    }

    /// Domain midpoint.
    #[inline]
    pub fn midpoint(&self) -> X {
        self.mid
    }

    /// Domain half-width.
    #[inline]
    pub fn half_width(&self) -> X {
        self.half
    }

    /// Normalize `x` to `tau` where `start -> -1`, midpoint `-> 0`, and `end -> 1`.
    ///
    /// Uses the form `(x - start) / half - 1` rather than `(x - mid) / half` to
    /// avoid catastrophic cancellation when evaluating exactly at the stored
    /// endpoints, where `x - start` is computed without loss of significance.
    #[inline]
    pub fn normalize(&self, x: X) -> f64 {
        (x - self.start) / self.half - 1.0
    }

    /// Whether `x` lies in the closed interval `[start, end]`.
    #[inline]
    pub fn contains(&self, x: X) -> bool {
        x >= self.start && x <= self.end
    }

    /// Map normalized `tau` back into this domain.
    #[inline]
    pub fn denormalize(&self, tau: f64) -> X {
        self.mid + self.half * tau
    }
}
