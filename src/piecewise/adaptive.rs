//! Adaptive segment tables.
//!
//! Segments are produced by recursive bisection: a candidate fit is computed
//! on each subdomain, its tail-error estimate is compared against the
//! tolerance, and on failure the subdomain is split at its midpoint. The
//! resulting segments are stored in ascending-start order with their
//! boundaries cached as a separate sorted `Vec<X>` so that
//! [`AdaptiveSegmentTable::evaluate`] can locate the active segment in
//! `O(log n)` via binary search rather than a linear scan.

use alloc::vec::Vec;

use crate::approx::error_estimation;
use crate::approx::fit;
use crate::core::{ChebyError, ChebyScalar, ChebySeries, ChebyTime, DifferentiateWith, Domain};
use crate::piecewise::ChebySegment;

/// Metadata for an adaptive segment.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SegmentMetadata<T> {
    /// Number of coefficients stored.
    pub degree: usize,
    /// Estimated local truncation error.
    pub estimated_error: Option<T>,
    /// Number of function samples used for the local fit.
    pub samples_used: usize,
    /// Recursion depth at which this segment was emitted.
    pub depth: usize,
}

/// An adaptive piecewise Chebyshev table over a single base domain.
///
/// Internally segments are stored in ascending order, plus a cached
/// `boundaries` vector of length `segments.len() + 1` (start of each segment
/// followed by the end of the last). Lookup uses [`slice::binary_search_by`].
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveSegmentTable<T, X, const N: usize> {
    domain: Domain<X>,
    segments: Vec<ChebySegment<T, X, N>>,
    metadata: Vec<SegmentMetadata<f64>>,
    boundaries: Vec<X>,
}

impl<T, X, const N: usize> AdaptiveSegmentTable<T, X, N>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    /// Build by recursively splitting until a tail estimate meets `tolerance`.
    pub fn from_fn(
        domain: Domain<X>,
        f: impl Fn(X) -> T + Copy,
        tolerance: f64,
        max_depth: usize,
    ) -> Result<Self, ChebyError> {
        if N < 2 {
            return Err(ChebyError::InvalidDegree);
        }
        if !(tolerance.is_finite() && tolerance >= 0.0) {
            return Err(ChebyError::NonFiniteInput);
        }
        let mut segments = Vec::new();
        let mut metadata = Vec::new();
        split::<T, X, N>(
            domain,
            f,
            tolerance,
            max_depth,
            0,
            &mut segments,
            &mut metadata,
        )?;
        if segments.is_empty() {
            return Err(ChebyError::EmptySegmentTable);
        }
        // Verify ordering and contiguity, then cache boundaries.
        let mut boundaries = Vec::with_capacity(segments.len() + 1);
        boundaries.push(segments[0].domain().start());
        for (i, seg) in segments.iter().enumerate() {
            if i > 0 && seg.domain().start() != segments[i - 1].domain().end() {
                return Err(ChebyError::SegmentBoundariesNotContiguous);
            }
            boundaries.push(seg.domain().end());
        }
        Ok(Self {
            domain,
            segments,
            metadata,
            boundaries,
        })
    }

    /// Base domain spanning every segment.
    #[inline]
    pub fn domain(&self) -> Domain<X> {
        self.domain
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

    /// Borrow the segment slice in ascending domain order.
    #[inline]
    pub fn segments(&self) -> &[ChebySegment<T, X, N>] {
        &self.segments
    }

    /// Per-segment metadata, in the same order as [`Self::segments`].
    #[inline]
    pub fn metadata(&self) -> &[SegmentMetadata<f64>] {
        &self.metadata
    }

    /// Cached monotone boundary vector. Length is `len() + 1`.
    #[inline]
    pub fn boundaries(&self) -> &[X] {
        &self.boundaries
    }

    /// Locate the segment containing `x` via binary search.
    pub fn locate(&self, x: X) -> Option<&ChebySegment<T, X, N>> {
        if self.segments.is_empty() {
            return None;
        }
        if x < self.boundaries[0] || x > self.boundaries[self.boundaries.len() - 1] {
            return None;
        }
        // Binary-search the boundary vector for the largest start <= x.
        let probe = self
            .boundaries
            .binary_search_by(|b| b.partial_cmp(&x).unwrap_or(core::cmp::Ordering::Less));
        let idx = match probe {
            Ok(i) => {
                if i == self.segments.len() {
                    self.segments.len() - 1
                } else {
                    i
                }
            }
            Err(i) => {
                if i == 0 {
                    0
                } else {
                    i - 1
                }
            }
        };
        self.segments.get(idx).filter(|s| s.contains(x))
    }

    /// Evaluate at `x`.
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
        self.locate(x)
            .ok_or(ChebyError::EvaluationOutOfDomain)?
            .evaluate(x)
    }

    /// Evaluate the dimensionally-correct derivative at `x`.
    pub fn evaluate_derivative(
        &self,
        x: X,
    ) -> Result<<T as DifferentiateWith<X>>::Derivative, ChebyError>
    where
        T: DifferentiateWith<X>,
    {
        self.locate(x)
            .ok_or(ChebyError::EvaluationOutOfDomain)?
            .evaluate_derivative(x)
    }
}

#[allow(clippy::too_many_arguments)]
fn split<T, X, const N: usize>(
    domain: Domain<X>,
    f: impl Fn(X) -> T + Copy,
    tolerance: f64,
    depth_remaining: usize,
    depth_so_far: usize,
    segments: &mut Vec<ChebySegment<T, X, N>>,
    metadata: &mut Vec<SegmentMetadata<f64>>,
) -> Result<(), ChebyError>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    let fitted = fit::fit_from_fn_on::<T, X, N>(domain, f);
    let coeffs = fitted.into_coeffs();
    let err = error_estimation::estimated_truncation_error(&coeffs);
    if err <= tolerance || depth_remaining == 0 {
        segments.push(ChebySegment::try_new(domain, ChebySeries::new(coeffs))?);
        metadata.push(SegmentMetadata {
            degree: N,
            estimated_error: Some(err),
            samples_used: N,
            depth: depth_so_far,
        });
        return Ok(());
    }
    let mid = domain.midpoint();
    split(
        Domain::try_new(domain.start(), mid)?,
        f,
        tolerance,
        depth_remaining - 1,
        depth_so_far + 1,
        segments,
        metadata,
    )?;
    split(
        Domain::try_new(mid, domain.end())?,
        f,
        tolerance,
        depth_remaining - 1,
        depth_so_far + 1,
        segments,
        metadata,
    )
}
