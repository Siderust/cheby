//! Chebyshev approximation and interpolation.
//!
//! # Theory
//!
//! Given a function `f` on a [`crate::core::Domain`], `fit` produces a
//! Chebyshev interpolating series by sampling `f` at roots of `T_N` and using
//! the discrete cosine orthogonality formula for the coefficients.
//! The optional [`adaptive`] module refines the degree until a target
//! tolerance is met, and the optional [`minimax`] module replaces the
//! root-sample fit with a Remez exchange that minimises the L∞ error.
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
