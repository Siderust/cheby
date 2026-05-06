//! Segment tables.

use alloc::vec::Vec;

use crate::approx::fit;
use crate::core::{ChebyError, ChebyScalar, ChebySeries, ChebyTime, DifferentiateWith, Domain};
use crate::piecewise::{lookup, ChebySegment};

/// A uniform table of Chebyshev segments.
///
/// The first segments have [`Self::segment_len`] width; the final segment may
/// be shorter when the requested table end is not an exact multiple of that
/// width. Table lookup treats the total range as half-open: `[start, end)`.
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
        let last_idx = segments.len() - 1;
        for (i, segment) in segments.iter().enumerate() {
            let expected = start + segment_len * i as f64;
            let domain = segment.domain();
            if domain.start() != expected {
                return Err(ChebyError::InvalidDomain);
            }
            let width = domain.end() - domain.start();
            if width <= X::zero() {
                return Err(ChebyError::NonPositiveSegmentLength);
            }
            if i == last_idx {
                if width > segment_len {
                    return Err(ChebyError::InvalidDomain);
                }
            } else if width != segment_len {
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
        if !start.is_finite() || !end.is_finite() || !segment_len.is_finite() {
            return Err(ChebyError::NonFiniteInput);
        }
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
        self.segments
            .last()
            .map(|segment| segment.domain().end())
            .unwrap_or(self.start)
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
        if self.segments.is_empty() || x < self.start || x >= self.end() {
            return None;
        }
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
