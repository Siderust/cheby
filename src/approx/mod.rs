//! Chebyshev approximation and interpolation.

pub mod error_estimation;
pub mod fit;
pub mod interpolation;

#[cfg(feature = "adaptive")]
pub mod adaptive;
#[cfg(feature = "minimax")]
pub mod minimax;
