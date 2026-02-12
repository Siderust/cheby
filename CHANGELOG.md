# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
