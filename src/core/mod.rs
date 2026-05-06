//! Core Chebyshev primitives.
//!
//! # Theory
//!
//! Every approximation in `cheby` is a finite Chebyshev expansion
//! `f(x) ≈ Σ aₖ Tₖ(τ)` over a normalised parameter `τ ∈ [-1, 1]`.
//! [`Domain`] handles the affine map between user coordinates and `τ`,
//! [`nodes`] returns Chebyshev abscissae (`Gauss` or `GaussLobatto`),
//! [`evaluate`] uses Clenshaw's recurrence, and [`ChebySeries`] /
//! [`ChebySeriesOn`] hold coefficients with optional unit safety via
//! [`ChebyScalar`].
//!
//! # Features
//!
//! Always available; `no_std` and allocation-free. Enable `alloc` for
//! the dynamic [`ChebySeriesDyn`] type.
//!
//! # Performance
//!
//! All hot paths are `const`-generic on the polynomial degree. Evaluation
//! is `O(N)` per query with no heap traffic; node generation is `O(N)`.

pub mod basis;
pub mod domain;
pub mod error;
pub mod eval;
pub mod nodes;
pub mod scalar;
pub mod series;

pub use domain::Domain;
pub use error::ChebyError;
pub use eval::{evaluate, evaluate_both};
pub use nodes::{nodes, nodes_mapped, NodeKind};
pub use scalar::{ChebyScalar, ChebyTime, DifferentiateWith, IntegrateWith};
pub use series::{ChebySeries, ChebySeriesOn};

#[cfg(feature = "alloc")]
pub use series::ChebySeriesDyn;
