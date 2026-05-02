//! Coefficient-table metadata.

/// Current binary table format version.
pub const FORMAT_VERSION: u32 = 1;

/// Metadata describing a coefficient table.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TableMetadata {
    /// Format version.
    pub format_version: u32,
    /// Polynomial degree.
    pub degree: usize,
    /// Segment count.
    pub segment_count: usize,
    /// Coefficient type name.
    pub coefficient_type: &'static str,
    /// Optional unit metadata.
    pub unit: Option<&'static str>,
    /// Simple checksum over payload bytes.
    pub checksum: u64,
}
