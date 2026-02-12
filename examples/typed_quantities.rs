//! Using `cheby` with typed quantities via `qtty`.
//!
//! Run with:
//! `cargo run --example typed_quantities`

use cheby::{evaluate, fit_from_fn};
use qtty::Quantity;

type Km = qtty::Kilometer;
type Kilometers = Quantity<Km>;

fn main() {
    const N: usize = 15;
    let start = -1.0;
    let end = 1.0;

    // Fit an altitude-like profile, but keep explicit length units.
    let coeffs: [Kilometers; N] = fit_from_fn(
        |t| Kilometers::new(7000.0 + 250.0 * (2.0 * std::f64::consts::PI * t).sin()),
        start,
        end,
    );

    // Evaluate at tau = 0.0 (the midpoint).
    let mid_value = evaluate(&coeffs, 0.0);
    println!("midpoint value = {:.6} km", mid_value.value());
}
