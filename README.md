# cheby

[![Crates.io](https://img.shields.io/crates/v/cheby.svg)](https://crates.io/crates/cheby)
[![Docs](https://docs.rs/cheby/badge.svg)](https://docs.rs/cheby)
[![Code Quality](https://github.com/Siderust/cheby/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Siderust/cheby/actions/workflows/ci.yml)

Chebyshev polynomial toolkit for scientific computing in Rust.

`cheby` provides a full interpolation pipeline suitable for numerical kernels
and ephemeris-style piecewise approximations:

1. Node generation on `[-1, 1]` and mapped intervals.
2. Coefficient fitting via DCT identities.
3. Stable Clenshaw evaluation (value, derivative, or both).
4. Uniform piecewise segment tables with O(1) segment lookup.

All core APIs are generic over `ChebyScalar`, so they work with `f64` and with
typed quantities such as `qtty::Quantity<Kilometer>`.

## Installation

```toml
[dependencies]
cheby = "0.1"
```

## Quick start

```rust
use cheby::{evaluate, fit_coeffs, nodes};

const N: usize = 9;
let xi: [f64; N] = nodes();
let values: [f64; N] = std::array::from_fn(|k| xi[k].sin());
let coeffs = fit_coeffs(&values);

let tau = 0.42;
let approx = evaluate(&coeffs, tau);
println!("sin({tau}) ~= {approx}");
```

## Segment-table workflow

```rust
use cheby::ChebySegmentTable;

let table: ChebySegmentTable<f64, 15> =
    ChebySegmentTable::from_fn(f64::sin, 0.0, std::f64::consts::TAU, std::f64::consts::FRAC_PI_2);

let (value, derivative) = table.eval_both(1.0).unwrap();
println!("f(t)={value}, df/dt={derivative}");
```

## Examples

Runnable examples are provided in `examples/`:

- `basic_interpolation`
- `segment_table`
- `typed_quantities`

Use:

```bash
cargo run --example basic_interpolation
```

See `examples/README.md` for details.

## Tests and coverage

`cheby` includes:

- Unit tests inside modules.
- Functional integration tests in `tests/functional_pipeline.rs`.
- Doctests for public examples.

Run locally:

```bash
cargo test --all-targets
cargo test --doc
cargo +nightly llvm-cov --workspace --all-features --doctests --summary-only
```

Coverage is gated in CI at **>= 90% line coverage**.

## CI

GitHub Actions workflow: `/.github/workflows/ci.yml`

Jobs:

- `check`
- `fmt`
- `clippy` (`-D warnings`)
- `test` (unit/integration + docs)
- `coverage` (`cargo +nightly llvm-cov`, PR summary, 90% gate)

## License

AGPL-3.0-only
