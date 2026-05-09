# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-05-09

### Changed

- The crate has been refactored from a flat source layout into feature-gated
  mathematical modules under `core`, `approx`, `calculus`, `piecewise`,
  `quadrature`, `spectral`, and `io`.
- The public positioning of the crate now centers on unit-safe Chebyshev
  approximation, interpolation, and spectral numerics. Public domain-specific
  APIs for trajectory or ephemeris-style use remain out of scope.
- `Cargo.toml` now exposes modular feature flags:
  `std`, `alloc`, `approx`, `adaptive`, `minimax`, `calculus`, `piecewise`,
  `quadrature`, `spectral`, `serde`, `binary`, `fast-dct`, and `full`.
- The README, examples, benchmarks, and CI configuration were rewritten to
  match the new crate architecture and release posture.

### Added

- `#![forbid(unsafe_code)]` at the crate root.
- Core abstractions:
  `Domain<X>`, `ChebyError`, `ChebySeries<T, N>`, `ChebySeriesDyn<T>`,
  `ChebySeriesOn<T, X, N>`, `ChebyScalar`, `ChebyTime`,
  `DifferentiateWith<X>`, and `IntegrateWith<X>`.
- Chebyshev basis functions in `core::basis`:
  `t(n, x)` and `u(n, x)`.
- First-class node families in `core::nodes`:
  `Roots`, `Extrema`, `Lobatto`, `Gauss`, and `GaussLobatto`.
- Stable Clenshaw evaluation in `core::eval`:
  `evaluate`, `evaluate_both`, and the compatibility `evaluate_derivative`.
- Approximation APIs in `approx`:
  coefficient fitting, function fitting on normalized and typed domains,
  barycentric interpolation, adaptive fitting with `FitReport`, coefficient-tail
  error estimates, and a Remez-style minimax interface.
- Calculus APIs for normalized derivatives and domain-aware definite integrals.
- Piecewise APIs:
  `ChebySegment`, `ChebySegmentTable`, uniform O(1) lookup, and adaptive
  segment-table construction with per-segment metadata.
- Quadrature APIs:
  Clenshaw-Curtis style weights and integration, plus Gauss-Chebyshev weighted
  integration.
- Spectral APIs:
  collocation points, a lightweight dense `Matrix`, and Chebyshev
  differentiation matrices.
- Optional IO support:
  serde-backed public type serialization, binary encoding helpers for dynamic
  `f64` series, and explicit table-format metadata.
- Property tests with `proptest` covering domain normalization, recurrence,
  interpolation reproduction, series derivative/integral consistency,
  piecewise lookup, and adaptive fitting.
- Compile-fail tests with `trybuild` covering mixed-unit domains and invalid
  derivative/integral assignments under `qtty`.
- New examples:
  `basic_series`, `fit_sin`, `adaptive_fit`, `minimax_exp`,
  `derivative_velocity`, `integral_position`, `piecewise_trajectory`,
  `ephemeris_like_table`, `angular_rate`, `star_alt_az_approximation`,
  `clenshaw_curtis_integral`, and `spectral_differentiation`.
- Criterion benchmarks for evaluation, derivative evaluation, fitting,
  piecewise lookup/evaluation, adaptive fitting, and Clenshaw-Curtis
  integration.
- Release-hardening files:
  `.github/workflows/ci.yml` and `deny.toml`.
- Integration test suite (`tests/coverage_boost.rs`) with 55 targeted tests
  raising line coverage from 75% to 91.65%, exceeding the 90% CI quality gate.

### Fixed

- Domain-aware derivatives are now dimensionally correct with `qtty`. A series
  over `Second` returning `Kilometer` now differentiates to
  `Kilometer / Second` instead of returning `Kilometer` with only a
  documentation caveat.
- Domain-aware integrals are now dimensionally correct with `qtty`. A series
  returning velocity-like quantities now integrates over time to a
  position-like quantity.
- Out-of-domain evaluation now consistently reports `ChebyError::EvaluationOutOfDomain`
  through fallible APIs instead of relying on implicit caller discipline.
- Segment-table construction now validates empty tables, segment lengths, and
  uniformity assumptions through typed errors.

- `Domain::normalize` now uses `(x - start) / half - 1` instead of
  `(x - mid) / half`, eliminating catastrophic cancellation when evaluating
  exactly at stored endpoints. `normalize(start)` and `normalize(end)` are
  now exact in IEEE 754 arithmetic.
- The `domain_normalization_properties` proptest now uses tighter input ranges
  (`|start| ≤ 1e3`, `width ≥ 1e-3`) consistent with the achievable floating-point
  precision of the normalization formula.

### Documentation

- The crate README now describes `cheby` as a unit-safe Chebyshev mathematics
  crate, documents the feature flags, and explicitly keeps application-domain
  use cases in examples or downstream crates.
- The changelog, examples, and CI docs now align with the refactored crate and
  the new dimensional-calculus behavior.

## [0.1.0 - 2026/02/12]

### Added

- GitHub Actions CI workflow at `.github/workflows/ci.yml` with:
  - `cargo check`, `cargo fmt --check`, `cargo clippy -- -D warnings`
  - unit/integration/doc tests
  - nightly `cargo-llvm-cov` coverage reporting and PR summary comments
  - coverage quality gate requiring at least 90% line coverage
- Functional integration tests in `tests/functional_pipeline.rs` covering:
  - full fit/evaluate roundtrip on canonical and mapped intervals
  - value/derivative API consistency
  - `ChebySegmentTable` metadata and boundary behavior
  - construction from precomputed segments and empty-table behavior
- New runnable examples:
  - `examples/basic_interpolation.rs`
  - `examples/segment_table.rs`
  - `examples/typed_quantities.rs`
  - plus documentation index at `examples/README.md`
- Expanded project README with:
  - installation and quick-start snippets
  - segment-table usage example
  - testing and coverage commands
  - CI and quality gate documentation
