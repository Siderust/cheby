# cheby examples

All examples are runnable with `cargo run --example <name>`.

## `basic_interpolation`

End-to-end interpolation on `[-1, 1]`:

- Generate Chebyshev nodes.
- Sample a function.
- Fit coefficients.
- Evaluate and print approximation error.

Run:

```bash
cargo run --example basic_interpolation
```

## `segment_table`

Piecewise approximation over a physical time range:

- Build `ChebySegmentTable`.
- Inspect metadata (`start`, `end`, segment count, segment length).
- Evaluate value and derivative in one pass.

Run:

```bash
cargo run --example segment_table
```

## `typed_quantities`

Demonstrates typed units with `qtty::Quantity`:

- Fit coefficients where values are `Quantity<Kilometer>`.
- Evaluate while preserving units.

Run:

```bash
cargo run --example typed_quantities
```
