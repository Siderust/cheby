//! Chebyshev approximation and interpolation.
//!
//! # Theory
//!
//! Given a function `f` on a [`crate::core::Domain`], `fit` produces a
//! truncated Chebyshev series whose coefficients minimise the L² error
//! at the Chebyshev nodes (interpolation in the polynomial sense).
//! The optional [`adaptive`] module refines the degree until a target
//! tolerance is met, and the optional [`minimax`] module replaces the
//! L² fit with a real Remez exchange that minimises the L∞ error.
//!
//! # Features
//!
//! - `approx` (default): fixed-size and dynamic L² fitting.
//! - `adaptive`: iterative refinement (`alloc` required).
//! - `minimax`: Remez exchange returning explicit convergence diagnostics
//!   (`alloc` required, `f64` only).
//!
//! # Performance
//!
//! Fixed-size fits are `O(N²)` with no allocations. Adaptive fits and
//! Remez allocate per iteration; both are `O(N²)` per pass.

pub mod error_estimation;
pub mod fit;
pub mod interpolation;

#[cfg(feature = "adaptive")]
pub mod adaptive;
#[cfg(feature = "minimax")]
pub mod minimax;
