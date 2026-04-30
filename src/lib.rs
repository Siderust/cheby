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
//!
//! # Numerical properties
//!
//! Chebyshev expansions are popular precisely because they are well-behaved
//! both in terms of *approximation power* and *floating-point conditioning*.
//! A few facts worth keeping in mind when using this crate:
//!
//! - **Exponential coefficient decay (analytic functions).**
//!   If `f` is analytic on `[-1, 1]` and admits an analytic continuation to
//!   the Bernstein ellipse `E_ρ` (foci ±1, semi-axis sum `ρ > 1`), then its
//!   Chebyshev coefficients satisfy `|a_n| ≤ M · ρ^{-n}` for some constant
//!   `M`. In practice this means coefficients fall off geometrically and
//!   the truncation error after keeping `N` of them is dominated by the
//!   first dropped term:
//!
//!   ```text
//!   |f(x) − Σ_{k=0}^{N-1} a_k T_k(x)|  ≲  Σ_{k≥N} |a_k|  ≈  |a_N|.
//!   ```
//!
//!   So a quick "did I keep enough terms?" check is to look at the magnitude
//!   of the last few coefficients returned by [`fit_coeffs`] / [`fit_from_fn`]:
//!   when they are at the level of `≈ ε · ‖f‖`, the series is converged.
//!
//! - **Clenshaw evaluation.** The [`eval`] module uses the Clenshaw
//!   recurrence rather than evaluating `Σ a_k T_k(x)` directly. This is
//!   numerically more stable for high degrees because it never forms the
//!   ill-conditioned monomial basis and the recurrence is a normwise stable
//!   way to sum a Chebyshev series. See [`evaluate`] for details.
//!
//! - **Domain rescaling.** Working on a general interval `[a, b]` is done
//!   via the standard affine map
//!
//!   ```text
//!   τ(x) = (2x − (a + b)) / (b − a),     τ ∈ [-1, 1].
//!   ```
//!
//!   The [`nodes_mapped`] / [`fit_from_fn`] / [`segment`] APIs apply this
//!   automatically. The caller is responsible for evaluating only at points
//!   `x ∈ [a, b]`: outside this interval `T_k(τ)` grows like
//!   `cosh(k · acosh|τ|)`, so extrapolation diverges *very* fast — this is
//!   not a polite "extrapolation is less accurate" warning, it is essentially
//!   exponential blow-up in the degree.
//!
//! - **Conditioning across `[-1, 1]`.** Evaluation is well-conditioned in
//!   the interior and slightly less so near the endpoints `τ = ±1`, but for
//!   typical engineering degrees (a few tens) the condition number stays
//!   modest and Clenshaw remains accurate to nearly machine precision.

mod eval;
mod fit;
mod nodes;
pub mod scalar;
pub mod segment;

pub use eval::{evaluate, evaluate_both, evaluate_derivative};
pub use fit::{fit_coeffs, fit_from_fn, fit_from_fn_t};
pub use nodes::{nodes, nodes_mapped, nodes_mapped_t};
pub use scalar::{ChebyScalar, ChebyTime};
pub use segment::{ChebySegment, ChebySegmentF, ChebySegmentTable, ChebySegmentTableF};
