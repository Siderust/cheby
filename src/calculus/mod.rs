//! Chebyshev calculus.
//!
//! # Theory
//!
//! Differentiation maps `Σ aₖ Tₖ` to `Σ a'ₖ Tₖ` via the standard upward
//! recurrence; integration runs the same recurrence in reverse and adds
//! a chosen constant. Both operators are exact on the polynomial space
//! and stable for the degrees `cheby` targets.
//!
//! # Features
//!
//! `calculus` (default), `no_std`, allocation-free.
//!
//! # Performance
//!
//! `O(N)` per series, in place, on the const-generic coefficient array.

pub mod derivative;
pub mod integral;
pub mod operators;
