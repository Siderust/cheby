//! Core Chebyshev primitives.

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
