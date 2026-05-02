//! Optional serialization helpers.

#[cfg(feature = "binary")]
pub mod binary;
pub mod metadata;
#[cfg(feature = "serde")]
pub mod serde;
