//! Piecewise Chebyshev approximations.

pub mod adaptive;
pub mod lookup;
pub mod segment;
pub mod table;

pub use adaptive::{AdaptiveSegmentTable, SegmentMetadata};
pub use segment::ChebySegment;
pub use table::ChebySegmentTable;
