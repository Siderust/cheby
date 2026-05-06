//! Compare the Chebyshev differentiation matrix against the analytic
//! derivative of `sin(x)` at the Lobatto collocation nodes.
//!
//! Run with `cargo run --example spectral_differentiation --features spectral`.

use cheby::spectral::{chebyshev_differentiation_matrix, collocation_points};

fn main() {
    const N: usize = 17;
    let nodes: [f64; N] = collocation_points::<N>();
    let f: [f64; N] = std::array::from_fn(|k| nodes[k].sin());
    let d = chebyshev_differentiation_matrix(N);

    let mut max_err = 0.0_f64;
    for (i, &node) in nodes.iter().enumerate() {
        let mut acc = 0.0;
        for (j, &fj) in f.iter().enumerate() {
            acc += d.get(i, j) * fj;
        }
        let err = (acc - node.cos()).abs();
        if err > max_err {
            max_err = err;
        }
    }
    println!("D·sin vs cos: max |error| = {max_err:.3e} on {N} Lobatto nodes");
    println!("matrix shape: {}x{}", d.rows(), d.cols());
}
