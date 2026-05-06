//! Optional serialization helpers.
//!
//! # Features
//!
//! - `binary`: a compact, versioned, checksummed `f64` coefficient layout
//!   (see [`binary`]).
//! - `serde`: derives that integrate with the `serde` ecosystem.
//!
//! Both submodules require `alloc`.

#[cfg(feature = "binary")]
pub mod binary;
pub mod metadata;
#[cfg(feature = "serde")]
pub mod serde;
