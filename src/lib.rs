// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev polynomial toolkit: node generation, coefficient fitting, and
//! Clenshaw evaluation — generic over scalar type.
//!
//! # Overview
//!
//! This crate provides the full pipeline for piecewise Chebyshev
//! interpolation, as used in JPL DE-series ephemerides and cached
//! lunar/planetary position evaluators:
//!
//! 1. **[`nodes`]** — Chebyshev node generation on `[-1, 1]` or mapped to
//!    an arbitrary interval.
//! 2. **[`fit`]** — DCT-based coefficient computation from function values
//!    at Chebyshev nodes.
//! 3. **[`eval`]** — Clenshaw-recurrence evaluation of a Chebyshev series
//!    (value, derivative, or both in one pass).
//! 4. **[`segment`]** — Piecewise Chebyshev approximation over uniform time
//!    segments, with automatic lookup and `t → τ` normalisation.
//!
//! All core functions are generic over [`ChebyScalar`], so they work with
//! raw `f64` as well as typed quantities (`qtty::Quantity<U>`).

mod eval;
mod fit;
mod nodes;
pub mod scalar;
pub mod segment;

pub use eval::{evaluate, evaluate_both, evaluate_derivative};
pub use fit::{fit_coeffs, fit_from_fn};
pub use nodes::{nodes, nodes_mapped};
pub use scalar::ChebyScalar;
pub use segment::{ChebySegment, ChebySegmentTable};
