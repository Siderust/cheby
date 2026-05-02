//! Using `cheby` with typed quantities via `qtty` — typed time *and* typed value.
//!
//! This example demonstrates the full typed-quantity pipeline:
//! - `Second` as the time domain type (`Tt`).
//! - `Kilometer` as the value type (`T`).
//! - `fit_from_fn_t` to sample a function with a typed time argument.
//! - `ChebySegment<Kilometer, Second, N>` for evaluation and derivative.
//!
//! The derivative result has Rust type `Kilometer` but physically represents
//! `Kilometer / Second` — the caller must know the time-unit context.
//!
//! Run with:
//! `cargo run --example typed_quantities`

use cheby::{fit_from_fn_t, ChebySegment};
use qtty::{Kilometer, Second};

fn main() {
    const N: usize = 15;

    let start = Second::new(0.0);
    let end = Second::new(std::f64::consts::TAU); // one full period
    let half = (end - start) * 0.5;
    let mid = start + half;

    // Fit an altitude-like profile over one second-typed period.
    // The closure receives a typed `Second` argument, so no `.value()` is
    // needed at the call site — but sin/cos require a raw f64, hence the
    // single extraction inside the closure.
    let coeffs: [Kilometer; N] = fit_from_fn_t(
        |t: Second| Kilometer::new(7000.0 + 250.0 * t.value().sin()),
        start,
        end,
    );

    let seg: ChebySegment<Kilometer, Second, N> = ChebySegment::new(coeffs, mid, half);

    // Evaluate value and derivative at a few typed time points.
    for raw_t in [0.5, 1.5, 3.0, 5.5] {
        let t = Second::new(raw_t);
        let (val, dval_dt) = seg.eval_both(t);
        // dval_dt has Rust type Kilometer, but physically is km/s.
        let exact_val = Kilometer::new(7000.0 + 250.0 * raw_t.sin());
        let exact_dval = 250.0 * raw_t.cos(); // km/s
        println!(
            "t={raw_t:.2} s: val={val:.6}  (exact {exact_val:.6})  \
             dval/dt={:.6} km/s  (exact {exact_dval:.6} km/s)",
            dval_dt.value() // extract raw f64: type is Kilometer but units are km/s
        );
    }
}
