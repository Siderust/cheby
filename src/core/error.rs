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
    /// Adjacent segments did not meet at a shared boundary.
    SegmentBoundariesNotContiguous,
    /// A binary blob is shorter than the minimum required header/footer.
    BinaryTooShort,
    /// A binary blob carries a format version this build does not support.
    UnsupportedFormatVersion {
        /// Version found in the blob.
        found: u32,
        /// Version expected by this build.
        expected: u32,
    },
    /// A binary blob's stored length does not match its payload size.
    BinaryLengthMismatch,
    /// A binary blob's checksum did not match its payload.
    BinaryChecksumMismatch,
    /// Remez exchange did not converge within the requested iteration budget.
    RemezDidNotConverge,
    /// Remez exchange could not find enough alternation points.
    RemezAlternationFailure,
}

#[cfg(feature = "std")]
impl std::fmt::Display for ChebyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyDomain => f.write_str("domain has zero width"),
            Self::InvalidDomain => f.write_str("domain is invalid"),
            Self::NonPositiveSegmentLength => f.write_str("segment length must be positive"),
            Self::InvalidDegree => f.write_str("degree is invalid"),
            Self::EmptyCoefficientSet => f.write_str("coefficient set is empty"),
            Self::EvaluationOutOfDomain => f.write_str("evaluation point is outside the domain"),
            Self::NonFiniteInput => f.write_str("input is not finite"),
            Self::EmptySegmentTable => f.write_str("segment table is empty"),
            Self::SegmentBoundariesNotContiguous => {
                f.write_str("segment boundaries are not contiguous")
            }
            Self::BinaryTooShort => f.write_str("binary blob is shorter than the format header"),
            Self::UnsupportedFormatVersion { found, expected } => {
                write!(
                    f,
                    "unsupported binary format version {found} (expected {expected})"
                )
            }
            Self::BinaryLengthMismatch => {
                f.write_str("binary blob length does not match the stored coefficient count")
            }
            Self::BinaryChecksumMismatch => f.write_str("binary blob checksum mismatch"),
            Self::RemezDidNotConverge => {
                f.write_str("Remez exchange did not converge within the iteration budget")
            }
            Self::RemezAlternationFailure => {
                f.write_str("Remez exchange could not find enough alternation points")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for ChebyError {}
