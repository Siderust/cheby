//! Adaptive segment tables.

use alloc::vec::Vec;

use crate::approx::error_estimation;
use crate::approx::fit;
use crate::core::{ChebyError, ChebyScalar, ChebySeries, ChebyTime, Domain};
use crate::piecewise::ChebySegment;

/// Metadata for an adaptive segment.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SegmentMetadata<T> {
    /// Number of coefficients stored.
    pub degree: usize,
    /// Estimated local error.
    pub estimated_error: Option<T>,
    /// Number of function samples used.
    pub samples_used: usize,
}

/// A simple adaptive segment table.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveSegmentTable<T, X, const N: usize> {
    segments: Vec<ChebySegment<T, X, N>>,
    metadata: Vec<SegmentMetadata<f64>>,
}

impl<T, X, const N: usize> AdaptiveSegmentTable<T, X, N>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    /// Build by recursively splitting until a tail estimate meets tolerance.
    pub fn from_fn(
        domain: Domain<X>,
        f: impl Fn(X) -> T + Copy,
        tolerance: f64,
        max_depth: usize,
    ) -> Result<Self, ChebyError> {
        let mut segments = Vec::new();
        let mut metadata = Vec::new();
        split::<T, X, N>(
            domain,
            f,
            tolerance,
            max_depth,
            &mut segments,
            &mut metadata,
        )?;
        Ok(Self { segments, metadata })
    }

    /// Evaluate by binary-searching segment starts.
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
        self.segments
            .iter()
            .find(|segment| segment.contains(x))
            .ok_or(ChebyError::EvaluationOutOfDomain)?
            .evaluate(x)
    }

    /// Segment metadata.
    pub fn metadata(&self) -> &[SegmentMetadata<f64>] {
        &self.metadata
    }
}

fn split<T, X, const N: usize>(
    domain: Domain<X>,
    f: impl Fn(X) -> T + Copy,
    tolerance: f64,
    depth: usize,
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
    if err <= tolerance || depth == 0 {
        segments.push(ChebySegment::try_new(domain, ChebySeries::new(coeffs))?);
        metadata.push(SegmentMetadata {
            degree: N,
            estimated_error: Some(err),
            samples_used: N,
        });
        return Ok(());
    }
    let mid = domain.midpoint();
    split(
        Domain::try_new(domain.start(), mid)?,
        f,
        tolerance,
        depth - 1,
        segments,
        metadata,
    )?;
    split(
        Domain::try_new(mid, domain.end())?,
        f,
        tolerance,
        depth - 1,
        segments,
        metadata,
    )?;
    Ok(())
}
