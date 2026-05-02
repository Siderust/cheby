//! Chebyshev spectral numerics.

pub mod collocation;
pub mod differentiation;
pub mod matrices;

pub use collocation::collocation_points;
pub use differentiation::chebyshev_differentiation_matrix;
pub use matrices::Matrix;
