//! Chebyshev spectral numerics.
//!
//! # Theory
//!
//! [`chebyshev_differentiation_matrix`] builds the dense `N×N` matrix `D`
//! such that `D · f` approximates `f'` at the Lobatto collocation nodes
//! returned by [`collocation_points`]. Convergence is spectral
//! (super-algebraic) for analytic integrands.
//!
//! # Features
//!
//! `spectral` (requires `calculus` and `alloc`).
//!
//! # Performance
//!
//! Matrix construction is `O(N²)`. A subsequent multiply against a sample
//! vector is also `O(N²)` and dominates the cost.

pub mod collocation;
pub mod differentiation;
pub mod matrices;

pub use collocation::collocation_points;
pub use differentiation::chebyshev_differentiation_matrix;
pub use matrices::Matrix;
