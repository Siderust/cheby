// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

#![forbid(unsafe_code)]
#![cfg_attr(not(feature = "std"), no_std)]
#![doc = include_str!("../README.md")]

#[cfg(feature = "alloc")]
extern crate alloc;

pub mod core;

#[cfg(feature = "approx")]
pub mod approx;
#[cfg(feature = "calculus")]
pub mod calculus;
#[cfg(any(feature = "serde", feature = "binary"))]
pub mod io;
#[cfg(feature = "piecewise")]
pub mod piecewise;
#[cfg(feature = "quadrature")]
pub mod quadrature;
#[cfg(feature = "spectral")]
pub mod spectral;

pub use core::{
    basis, evaluate, evaluate_both, ChebyError, ChebyScalar, ChebySeries, ChebySeriesOn, ChebyTime,
    Domain, NodeKind,
};

#[cfg(feature = "alloc")]
pub use core::ChebySeriesDyn;

#[cfg(feature = "approx")]
pub use approx::fit::fit_coeffs;

#[cfg(feature = "piecewise")]
pub use piecewise::{ChebySegment, ChebySegmentTable};

/// Backward-compatible normalized derivative evaluation.
///
/// This returns `df/dtau`, where `tau` is the normalized coordinate in
/// `[-1, 1]`. Domain-aware derivatives with physical units are available on
/// [`ChebySeriesOn`] and [`piecewise::ChebySegment`].
#[inline]
pub fn evaluate_derivative<T: ChebyScalar>(coeffs: &[T], tau: f64) -> T {
    core::eval::evaluate_derivative(coeffs, tau)
}

/// Backward-compatible root helper for roots of `T_N` on `[-1, 1]`.
#[inline]
pub fn chebyshev_roots<const N: usize>() -> [f64; N] {
    core::nodes::nodes(NodeKind::Roots)
}

/// Backward-compatible root-node helper.
#[inline]
pub fn nodes<const N: usize>() -> [f64; N] {
    core::nodes::nodes(NodeKind::Roots)
}

/// Backward-compatible mapped root-node helper.
#[inline]
pub fn nodes_mapped<const N: usize>(start: f64, end: f64) -> [f64; N] {
    let domain = Domain::new(start, end);
    core::nodes::nodes_mapped(domain, NodeKind::Roots)
}

/// Backward-compatible mapped roots helper.
#[inline]
pub fn nodes_mapped_interval<const N: usize>(start: f64, end: f64) -> [f64; N] {
    let domain = Domain::new(start, end);
    core::nodes::nodes_mapped(domain, NodeKind::Roots)
}

/// Backward-compatible typed mapped roots helper.
#[inline]
pub fn nodes_mapped_t<Tt: ChebyTime, const N: usize>(start: Tt, end: Tt) -> [Tt; N] {
    let domain = Domain::new(start, end);
    core::nodes::nodes_mapped(domain, NodeKind::Roots)
}

/// Backward-compatible root-node name.
#[inline]
pub fn nodes_roots<const N: usize>() -> [f64; N] {
    core::nodes::nodes(NodeKind::Roots)
}

/// Backward-compatible fitting helper returning raw coefficients.
#[cfg(feature = "approx")]
#[inline]
pub fn fit_from_fn<T: ChebyScalar, const N: usize>(
    f: impl Fn(f64) -> T,
    start: f64,
    end: f64,
) -> [T; N] {
    approx::fit::fit_from_fn_on(Domain::new(start, end), f).into_coeffs()
}

/// Backward-compatible typed fitting helper returning raw coefficients.
#[cfg(feature = "approx")]
#[inline]
pub fn fit_from_fn_t<T: ChebyScalar, X: ChebyTime, const N: usize>(
    f: impl Fn(X) -> T,
    start: X,
    end: X,
) -> [T; N] {
    approx::fit::fit_from_fn_on(Domain::new(start, end), f).into_coeffs()
}
