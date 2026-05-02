//! Segment tables.

use alloc::vec::Vec;

use crate::approx::fit;
use crate::core::{ChebyError, ChebyScalar, ChebySeries, ChebyTime, DifferentiateWith, Domain};
use crate::piecewise::{lookup, ChebySegment};

/// A uniform table of Chebyshev segments.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChebySegmentTable<T, X, const N: usize> {
    start: X,
    segment_len: X,
    segments: Vec<ChebySegment<T, X, N>>,
}

impl<T, X, const N: usize> ChebySegmentTable<T, X, N>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    /// Build from precomputed uniform segments.
    pub fn try_from_segments(segments: Vec<ChebySegment<T, X, N>>) -> Result<Self, ChebyError> {
        if segments.is_empty() {
            return Err(ChebyError::EmptySegmentTable);
        }
        let start = segments[0].domain().start();
        let segment_len = segments[0].domain().end() - start;
        if segment_len <= X::zero() {
            return Err(ChebyError::NonPositiveSegmentLength);
        }
        for (i, segment) in segments.iter().enumerate() {
            let expected = start + segment_len * i as f64;
            if segment.domain().start() != expected {
                return Err(ChebyError::InvalidDomain);
            }
        }
        Ok(Self {
            start,
            segment_len,
            segments,
        })
    }

    /// Backward-compatible construction from segments and explicit table metadata.
    #[inline]
    pub fn from_segments(segments: Vec<ChebySegment<T, X, N>>, start: X, segment_len: X) -> Self {
        Self {
            start,
            segment_len,
            segments,
        }
    }

    /// Build a uniform table by fitting each segment.
    pub fn from_fn(f: impl Fn(X) -> T, start: X, end: X, segment_len: X) -> Self {
        Self::try_from_fn(f, start, end, segment_len).expect("invalid ChebySegmentTable domain")
    }

    /// Fallible construction from a function.
    pub fn try_from_fn(
        f: impl Fn(X) -> T,
        start: X,
        end: X,
        segment_len: X,
    ) -> Result<Self, ChebyError> {
        if segment_len <= X::zero() {
            return Err(ChebyError::NonPositiveSegmentLength);
        }
        let span = (end - start) / segment_len;
        if span <= 0.0 {
            return Err(ChebyError::InvalidDomain);
        }
        let count = span.ceil() as usize;
        let mut segments = Vec::with_capacity(count);
        for i in 0..count {
            let seg_start = start + segment_len * i as f64;
            let seg_end = if i + 1 == count {
                end
            } else {
                seg_start + segment_len
            };
            let domain = Domain::try_new(seg_start, seg_end)?;
            let coeffs = fit::fit_from_fn_on::<T, X, N>(domain, &f).into_coeffs();
            segments.push(ChebySegment::try_new(domain, ChebySeries::new(coeffs))?);
        }
        Ok(Self {
            start,
            segment_len,
            segments,
        })
    }

    /// Number of segments.
    #[inline]
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Whether the table is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    /// Table start.
    #[inline]
    pub fn start(&self) -> X {
        self.start
    }

    /// Table end.
    #[inline]
    pub fn end(&self) -> X {
        self.start + self.segment_len * self.segments.len() as f64
    }

    /// Uniform segment length used for O(1) lookup.
    #[inline]
    pub fn segment_len(&self) -> X {
        self.segment_len
    }

    /// Segment slice.
    #[inline]
    pub fn segments(&self) -> &[ChebySegment<T, X, N>] {
        &self.segments
    }

    /// Look up a segment.
    #[inline]
    pub fn get_segment(&self, x: X) -> Option<&ChebySegment<T, X, N>> {
        lookup::uniform_index(self.start, self.segment_len, self.segments.len(), x)
            .and_then(|idx| self.segments.get(idx))
            .filter(|segment| segment.contains(x))
    }

    /// Evaluate at `x`.
    #[inline]
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
        self.get_segment(x)
            .ok_or(ChebyError::EvaluationOutOfDomain)?
            .evaluate(x)
    }

    /// Backward-compatible optional evaluation.
    #[inline]
    pub fn eval(&self, x: X) -> Option<T> {
        self.evaluate(x).ok()
    }

    /// Evaluate derivative at `x`.
    #[inline]
    pub fn evaluate_derivative(
        &self,
        x: X,
    ) -> Result<<T as DifferentiateWith<X>>::Derivative, ChebyError>
    where
        T: DifferentiateWith<X>,
    {
        self.get_segment(x)
            .ok_or(ChebyError::EvaluationOutOfDomain)?
            .evaluate_derivative(x)
    }

    /// Backward-compatible optional derivative evaluation.
    #[inline]
    pub fn eval_derivative(&self, x: X) -> Option<<T as DifferentiateWith<X>>::Derivative>
    where
        T: DifferentiateWith<X>,
    {
        self.evaluate_derivative(x).ok()
    }

    /// Evaluate value and derivative at `x`.
    #[inline]
    pub fn eval_both(&self, x: X) -> Option<(T, <T as DifferentiateWith<X>>::Derivative)>
    where
        T: DifferentiateWith<X>,
    {
        self.get_segment(x)
            .and_then(|segment| segment.evaluate_both(x).ok())
    }
}
