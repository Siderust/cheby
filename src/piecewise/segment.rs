//! Single Chebyshev segment.

use crate::core::{ChebyError, ChebyScalar, ChebySeries, ChebyTime, DifferentiateWith, Domain};

/// A Chebyshev series valid on one physical domain segment.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChebySegment<T, X, const N: usize> {
    domain: Domain<X>,
    series: ChebySeries<T, N>,
}

impl<T, X, const N: usize> ChebySegment<T, X, N>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    /// Create a segment from a domain and normalized series.
    pub fn try_new(domain: Domain<X>, series: ChebySeries<T, N>) -> Result<Self, ChebyError> {
        if N == 0 {
            return Err(ChebyError::EmptyCoefficientSet);
        }
        Ok(Self { domain, series })
    }

    /// Backward-compatible constructor from coefficients, midpoint, and half-width.
    #[inline]
    pub fn new(coeffs: [T; N], mid: X, half: X) -> Self {
        let domain = Domain::new(mid - half, mid + half);
        Self {
            domain,
            series: ChebySeries::new(coeffs),
        }
    }

    /// Segment domain.
    #[inline]
    pub fn domain(&self) -> Domain<X> {
        self.domain
    }

    /// Coefficients.
    #[inline]
    pub fn coeffs(&self) -> &[T; N] {
        self.series.coeffs()
    }

    /// Midpoint.
    #[inline]
    pub fn mid(&self) -> X {
        self.domain.midpoint()
    }

    /// Half-width.
    #[inline]
    pub fn half(&self) -> X {
        self.domain.half_width()
    }

    /// Normalize a coordinate into this segment.
    #[inline]
    pub fn normalise(&self, x: X) -> f64 {
        self.domain.normalize(x)
    }

    /// American spelling alias.
    #[inline]
    pub fn normalize(&self, x: X) -> f64 {
        self.normalise(x)
    }

    /// Whether `x` is in this segment.
    #[inline]
    pub fn contains(&self, x: X) -> bool {
        self.domain.contains(x)
    }

    /// Fallible evaluation.
    #[inline]
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
        if !self.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        Ok(self.series.evaluate(self.domain.normalize(x)))
    }

    /// Backward-compatible infallible evaluation.
    #[inline]
    pub fn eval(&self, x: X) -> T {
        self.evaluate(x)
            .expect("ChebySegment evaluation outside domain")
    }

    /// Fallible, dimensionally correct derivative.
    #[inline]
    pub fn evaluate_derivative(
        &self,
        x: X,
    ) -> Result<<T as DifferentiateWith<X>>::Derivative, ChebyError>
    where
        T: DifferentiateWith<X>,
    {
        if !self.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        let d_tau = self.series.evaluate_both(self.domain.normalize(x)).1;
        Ok((d_tau * 2.0).scale_derivative(self.domain.end() - self.domain.start()))
    }

    /// Backward-compatible derivative alias.
    #[inline]
    pub fn eval_derivative(&self, x: X) -> <T as DifferentiateWith<X>>::Derivative
    where
        T: DifferentiateWith<X>,
    {
        self.evaluate_derivative(x)
            .expect("ChebySegment derivative outside domain")
    }

    /// Fallible value and derivative evaluation.
    #[inline]
    pub fn evaluate_both(
        &self,
        x: X,
    ) -> Result<(T, <T as DifferentiateWith<X>>::Derivative), ChebyError>
    where
        T: DifferentiateWith<X>,
    {
        if !self.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        let (value, d_tau) = self.series.evaluate_both(self.domain.normalize(x));
        Ok((
            value,
            (d_tau * 2.0).scale_derivative(self.domain.end() - self.domain.start()),
        ))
    }

    /// Backward-compatible value and derivative alias.
    #[inline]
    pub fn eval_both(&self, x: X) -> (T, <T as DifferentiateWith<X>>::Derivative)
    where
        T: DifferentiateWith<X>,
    {
        self.evaluate_both(x)
            .expect("ChebySegment value/derivative outside domain")
    }
}
