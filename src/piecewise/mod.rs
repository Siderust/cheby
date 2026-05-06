//! Piecewise Chebyshev approximations.
//!
//! # Theory
//!
//! A piecewise table partitions a domain into segments, fits a Chebyshev
//! series on each, and dispatches lookups by segment. [`ChebySegmentTable`]
//! uses a uniform partition; [`AdaptiveSegmentTable`] subdivides by error
//! tolerance and stores the resulting monotone segment boundaries for
//! `O(log S)` lookup.
//!
//! # Features
//!
//! `piecewise` (requires `alloc`).
//!
//! # Performance
//!
//! Uniform tables locate the active segment in `O(1)`; adaptive tables in
//! `O(log S)` via [`core::slice::binary_search_by`]. Per-segment evaluation
//! is `O(N)` Clenshaw.

pub mod adaptive;
pub mod lookup;
pub mod segment;
pub mod table;

pub use adaptive::{AdaptiveSegmentTable, SegmentMetadata};
pub use segment::ChebySegment;
pub use table::ChebySegmentTable;
