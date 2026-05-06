# cheby examples

All examples are runnable with `cargo run --example <name> [--features ...]`.
Examples whose subsystem requires a Cargo feature have a `required-features`
entry in `Cargo.toml`; the table at the bottom maps each example to the
features it needs.

## Core evaluation and fitting

### `basic_series`
Build a `ChebySeries` by hand and evaluate it.

### `basic_interpolation`
End-to-end interpolation on `[-1, 1]`: nodes, sample, fit, evaluate, error.

### `fit_sin`
Fixed-size fit of `sin` over a typed time domain.

### `typed_quantities`
Fit and evaluate while preserving `qtty::Quantity` units.

## Derivatives and integrals (`calculus`)

### `derivative_velocity`
Differentiate a position series to recover a typed velocity.

### `integral_position`
Integrate a velocity series to recover a typed position.

## Adaptive and minimax fitting (`adaptive`, `minimax`)

### `adaptive_fit`
Refine a fit until a target tolerance is met (`adaptive` feature).

### `minimax_exp`
Minimise the L∞ error for `exp` using Remez exchange (`minimax` feature).

## Piecewise tables (`piecewise`)

### `segment_table`
Uniform piecewise table over a typed domain with metadata.

### `piecewise_trajectory`
A typed multi-segment trajectory with per-segment evaluation.

### `ephemeris_like_table`
Builds a piecewise table reminiscent of an ephemeris layout.

## Quadrature (`quadrature`)

### `clenshaw_curtis_integral`
Clenshaw-Curtis quadrature of `sin` on `[0, π]`.

### `gauss_chebyshev`
Gauss-Chebyshev recovers `π` and `π·J₀(1)`.

## Spectral differentiation (`spectral`)

### `spectral_differentiation`
Build the Chebyshev differentiation matrix and verify it against
`d/dx sin(x) = cos(x)` at the Lobatto nodes.

## Binary I/O (`binary`)

### `binary_roundtrip`
Encode an `f64` series, decode it, and demonstrate that a flipped byte is
rejected by the checksum.

## Application-flavoured demos

These illustrate downstream usage and rely on default features only:

### `angular_rate`
Differentiate an angle series to recover an angular rate.

### `star_alt_az_approximation`
Toy approximation of an alt/az curve using a piecewise table.

## Feature-to-example map

| Feature       | Examples |
|---------------|----------|
| `std` (default) | `basic_series`, `basic_interpolation`, `fit_sin`, `typed_quantities`, `angular_rate` |
| `calculus` (default) | `derivative_velocity`, `integral_position` |
| `piecewise` (default) | `segment_table`, `piecewise_trajectory`, `ephemeris_like_table`, `star_alt_az_approximation` |
| `adaptive`     | `adaptive_fit` |
| `minimax`      | `minimax_exp` |
| `quadrature`   | `clenshaw_curtis_integral`, `gauss_chebyshev` |
| `spectral`     | `spectral_differentiation` |
| `binary`       | `binary_roundtrip` |
