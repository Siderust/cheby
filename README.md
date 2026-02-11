# cheby

Chebyshev polynomial toolkit for scientific computing in Rust.

Provides the full Chebyshev interpolation pipeline:

1. **Node generation** — Chebyshev nodes on `[-1, 1]` or mapped to an arbitrary interval.
2. **Coefficient fitting** — DCT-based coefficient computation from function values at Chebyshev nodes.
3. **Clenshaw evaluation** — Numerically stable evaluation of Chebyshev series (value, derivative, or both in one pass).
4. **Segment management** — Piecewise Chebyshev approximation over uniform time segments, with automatic lookup and `t → τ` normalization.

All core functions are generic over a `ChebyScalar` trait, so they work with
raw `f64` as well as typed quantities (e.g., `qtty::Quantity<Kilometer>`).

## Quick start

```rust
use cheby::{nodes, fit_coeffs, evaluate, evaluate_both};

// 1. Generate 9 Chebyshev nodes on [-1, 1]
let xi: [f64; 9] = nodes();

// 2. Sample a function at those nodes
let values: [f64; 9] = std::array::from_fn(|k| xi[k].sin());

// 3. Fit Chebyshev coefficients
let coeffs = fit_coeffs(&values);

// 4. Evaluate at an arbitrary point
let approx = evaluate(&coeffs, 0.42);
let (val, deriv) = evaluate_both(&coeffs, 0.42);
```

## License

AGPL-3.0
