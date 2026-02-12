//! Basic Chebyshev interpolation pipeline on `[-1, 1]`.
//!
//! Run with:
//! `cargo run --example basic_interpolation`

use cheby::{evaluate, fit_coeffs, nodes};

fn main() {
    const N: usize = 13;

    // 1) Generate canonical Chebyshev nodes.
    let xi: [f64; N] = nodes();

    // 2) Sample a target function at those nodes.
    let values: [f64; N] = std::array::from_fn(|k| (2.0 * xi[k]).sin() + 0.25 * xi[k]);

    // 3) Fit Chebyshev coefficients.
    let coeffs = fit_coeffs(&values);

    // 4) Evaluate at arbitrary tau values in [-1, 1].
    for &tau in &[-0.8, -0.2, 0.3, 0.9] {
        let approx = evaluate(&coeffs, tau);
        let exact = (2.0 * tau).sin() + 0.25 * tau;
        println!(
            "tau={tau:+.2}: approx={approx:.12}, exact={exact:.12}, abs_err={:.3e}",
            (approx - exact).abs()
        );
    }
}
