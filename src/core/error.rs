//! Error types used by fallible Chebyshev APIs.

/// Errors returned by Chebyshev constructors and fallible evaluations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ChebyError {
    /// The requested domain has no extent.
    EmptyDomain,
    /// The requested domain is reversed, non-finite, or otherwise invalid.
    InvalidDomain,
    /// A segment length was zero or negative.
    NonPositiveSegmentLength,
    /// A requested degree is invalid for the operation.
    InvalidDegree,
    /// No coefficients were supplied where at least one coefficient is required.
    EmptyCoefficientSet,
    /// Evaluation was requested outside the fitted domain.
    EvaluationOutOfDomain,
    /// An input value was NaN or infinite.
    NonFiniteInput,
    /// A segment table was constructed without segments.
    EmptySegmentTable,
}

#[cfg(feature = "std")]
impl std::fmt::Display for ChebyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let message = match self {
            Self::EmptyDomain => "domain has zero width",
            Self::InvalidDomain => "domain is invalid",
            Self::NonPositiveSegmentLength => "segment length must be positive",
            Self::InvalidDegree => "degree is invalid",
            Self::EmptyCoefficientSet => "coefficient set is empty",
            Self::EvaluationOutOfDomain => "evaluation point is outside the domain",
            Self::NonFiniteInput => "input is not finite",
            Self::EmptySegmentTable => "segment table is empty",
        };
        f.write_str(message)
    }
}

#[cfg(feature = "std")]
impl std::error::Error for ChebyError {}
