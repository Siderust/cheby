//! Chebyshev-related quadrature.
//!
//! # Theory
//!
//! [`clenshaw_curtis`] integrates over `[-1, 1]` with weights derived from
//! the discrete cosine transform on Chebyshev-Lobatto nodes; it is exact
//! on polynomials of degree up to `N - 1`. [`gauss_chebyshev`] integrates
//! against the weight `1/√(1−x²)` using the closed-form Chebyshev nodes
//! and uniform weight `π/N`.
//!
//! # Features
//!
//! `quadrature`, `no_std`, allocation-free.
//!
//! # Performance
//!
//! Both rules are `O(N²)` to build the weight set and `O(N)` per integrand
//! evaluation.

pub mod clenshaw_curtis;
pub mod gauss_chebyshev;
