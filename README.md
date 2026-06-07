# cheby

`cheby`: unit-safe Chebyshev approximation and spectral numerics for Rust.

This crate is a reusable mathematical toolkit for Chebyshev basis functions,
nodes, series evaluation, interpolation, approximation, calculus, piecewise
tables, quadrature, and spectral numerics. It intentionally does not expose
public modules for physics, astronomy, trajectories, ephemerides, mechanics,
robotics, signal filters, or other application domains. Domain-specific uses
belong in `examples/` or downstream crates.

## What is this crate?

`cheby` provides numerical primitives centered on Chebyshev polynomials:

- stable Clenshaw evaluation,
- root and Lobatto node families,
- fixed and dynamic coefficient series,
- typed domains that map physical coordinates to `[-1, 1]`,
- coefficient fitting and barycentric interpolation,
- adaptive fitting and tail-based error estimates,
- dimensionally correct derivatives and integrals with `qtty`,
- piecewise segment tables,
- Chebyshev-related quadrature,
- collocation points and differentiation matrices.

## Core concepts

The core API is always available:

```rust
use cheby::{ChebySeries, Domain};

let domain = Domain::new(0.0, 2.0);
let series = ChebySeries::new([1.0, 2.0, 0.5]);
let value = series.evaluate(domain.normalize(1.25));
```

`Domain<X>` maps `start -> -1`, midpoint `-> 0`, and `end -> 1`.
`ChebySeries<T, N>` stores coefficients of `T_0..T_{N-1}`.

## Mathematical model

A Chebyshev series approximates a function on a finite interval `[a, b]` by
first mapping the physical coordinate `x` into a dimensionless coordinate
`tau` on `[-1, 1]`:

```text
mid  = (a + b) / 2
half = (b - a) / 2
tau  = (x - mid) / half
x    = mid + half * tau
```

The approximating polynomial is stored in the first-kind Chebyshev basis:

```text
p(x) = sum_j a_j T_j(tau)
T_0(tau) = 1
T_1(tau) = tau
T_{j+1}(tau) = 2 tau T_j(tau) - T_{j-1}(tau)
```

`fit_coeffs` expects samples at the roots of `T_N`:

```text
tau_k = cos(pi * (2k + 1) / (2N)),  k = 0..N-1
a_0   = (1/N) sum_k f(tau_k)
a_j   = (2/N) sum_k f(tau_k) cos(j pi (2k + 1) / (2N)), j > 0
```

Evaluation uses Clenshaw recurrence, which evaluates the same series without
expanding it into monomials. This is usually more stable and faster than
building powers `tau^j` explicitly.

For a domain-aware series, derivatives and integrals apply the affine scale:

```text
dp/dx = (1 / half) dp/dtau = (2 / (b - a)) dp/dtau
Integral from a to x = half * Integral from -1 to tau of p(s) ds
```

Use this crate when a smooth function on a bounded interval will be evaluated
many times: ephemeris segments, trajectory tables, calibration curves,
spectral collocation, quadrature rules, and compact tabulation of expensive
functions are typical applications.

## Why Chebyshev?

Chebyshev expansions are well conditioned on finite intervals, avoid the worst
Runge behavior of uniform interpolation, and often converge rapidly for smooth
or analytic functions. `cheby` evaluates them with Clenshaw recurrence instead
of unstable monomial expansion.

## Unit safety with qtty

Domain variables and values can be `qtty` quantities. Mixed-unit domains fail
at compile time, and domain-aware calculus returns dimensionally correct types.

```rust
use qtty::{Kilometer, Second};

let domain = cheby::Domain::new(Second::new(0.0), Second::new(10.0));
let height = cheby::approx::fit::fit_from_fn_on::<Kilometer, Second, 16>(
    domain,
    |t| Kilometer::new(t.value() * t.value()),
);
let velocity = height.evaluate_derivative(Second::new(4.0)).unwrap();
```

The derivative type is `Kilometer / Second`, not `Kilometer`.

## Interpolation

`approx::interpolation::BarycentricInterpolator` evaluates interpolants from
sampled nodes and values while preserving value units.

## Approximation

`approx::fit` fits coefficients from samples or functions. `adaptive` adds a
builder that tries increasing coefficient counts up to the caller's
`max_degree` until coefficient-tail tolerances are met.
`minimax` exposes a Remez-style interface with convergence diagnostics.

## Piecewise tables

`piecewise::ChebySegment` ties a series to a domain. `ChebySegmentTable`
stores uniform segments with O(1) lookup and explicit out-of-domain errors.

## Derivatives and integrals

Normalized derivatives are available from `core::eval`. Domain-aware series
and segments apply the chain rule and use `qtty` multiplication/division so:

- position over time differentiates to velocity,
- velocity over time integrates to position,
- angular rate over time integrates to angle.

## Root finding

With the `alloc` feature, dynamic series support real root finding on the
normalized interval `[-1, 1]`:

```rust
use cheby::{fit_dyn_from_fn, RootOptions};

let series = fit_dyn_from_fn(12, |tau| tau * tau - 0.25).unwrap();
let roots = series.roots(); // normalized tau in [-1, 1]

let opts = RootOptions {
    zero_tol: 1e-10,
    ..RootOptions::default()
};
let refined = series.roots_with(opts);
```

The solver converts the Chebyshev expansion to a monic power polynomial and finds
roots with a **Durand–Kerner** iteration, then verifies residuals on `[-1, 1]`.
When the primary pass misses roots (tangency, repeated roots), a fallback
Chebyshev-node scan uses Brent refinement and minimum search. Constant or
identically-zero series return an empty root list.

[`RootOptions`] controls bracket width (`unit_tol`), residual tolerance
(`zero_tol`), and deduplication (`dedupe_eps`). Root finding is intended for
moderate polynomial degrees and can be less stable than Clenshaw evaluation at
high degree.

Map normalized roots to a physical domain with [`ChebySeriesDynOn::roots`],
or shift the constant term without refitting via
[`ChebySeriesDyn::shifted_constant`] when searching multiple level sets of the
same fitted signal.

## Quadrature

`quadrature` contains Chebyshev-related rules only: Clenshaw-Curtis style
integration and Gauss-Chebyshev weighted integration.

## Spectral methods

`spectral` contains collocation points, differentiation matrices, and matrix
storage. It does not include full PDE or physics solvers.

## Examples

Application-shaped examples live under `examples/`, including typed derivative,
typed integral, piecewise trajectory-style data, angular-rate integration,
Clenshaw-Curtis integration, and spectral differentiation.

## Feature flags

Default features are `std`, `approx`, `calculus`, and `piecewise`.

- `alloc`: dynamic series, adaptive methods, and piecewise tables.
- `adaptive`: adaptive fitting.
- `minimax`: Remez-style approximation.
- `quadrature`: Chebyshev quadrature.
- `spectral`: collocation and differentiation matrices.
- `serde`: serde support for supported table/metadata types.
- `binary`: compact binary helpers for coefficient tables.
- `full`: all crate features.

No domain-specific feature flags are provided.

## Safety guarantees

The crate uses `#![forbid(unsafe_code)]`. Fallible constructors return
`ChebyError`; panicking compatibility helpers are named or documented as such.

## Validation

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-targets --all-features
cargo test --doc --all-features
```

Additional release checks include `cargo audit`, `cargo deny`, `cargo llvm-cov`,
and Miri for core functionality when available.
