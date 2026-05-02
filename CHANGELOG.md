# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `evaluate_derivative` now delegates to `evaluate_both`, eliminating a
  duplicate Clenshaw recurrence loop. Behaviour is identical; any future
  fix to the recurrence only needs to be made in one place.
- `nodes_mapped` now delegates to `nodes_mapped_t` (valid since `f64: ChebyTime`),
  removing a manual reimplementation of the same affine-map loop.
- `fit_from_fn` now delegates to `fit_from_fn_t`, removing a manual
  reimplementation of the same sample-and-fit pipeline.
- `ChebySegment` fields (`coeffs`, `mid`, `half`, `half_inv`) are no longer
  `pub`. Read accessors `coeffs()`, `mid()`, `half()`, `half_inv()` are
  provided instead. This prevents external code from writing `seg.half`
  without updating the precomputed `half_inv`, which would silently corrupt
  derivative results.
- `examples/typed_quantities.rs` rewritten to demonstrate the full
  typed-quantity pipeline: `Second` as the time type and `Kilometer` as the
  value type, using `fit_from_fn_t` and `ChebySegment<Kilometer, Second, N>`.
  The derivative output and its implicit `km/s` units are explicitly shown.

### Documentation

- `fit_coeffs`: added `# Complexity` section noting the O(N²) cost of the
  naive DCT and recommending an FFT-based approach for N > ~30.
- `ChebyTime` trait and `recip_f64`: documented the dimension gap — for typed
  time `Quantity<U>`, `recip_f64` returns a raw `f64` that is dimensionally
  `1/[U]`, so the Rust return type of `eval_derivative` / `eval_both` is `T`
  while the physical units are `[T]/[Tt]`. The caller must know the time-unit
  context when interpreting derivative results.
- `ChebySegment::eval_derivative` and `eval_both`: carry the same derivative
  units caveat inline.

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
