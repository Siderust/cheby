//! Gauss-Chebyshev quadrature recovers two known weighted integrals:
//!
//! - `∫_{-1}^{1} 1/√(1−x²) dx = π`
//! - `∫_{-1}^{1} cos(x)/√(1−x²) dx = π J₀(1)`
//!
//! Run with `cargo run --example gauss_chebyshev --features quadrature`.

use cheby::quadrature::gauss_chebyshev;

fn main() {
    const N: usize = 32;
    let pi_estimate: f64 = gauss_chebyshev::integrate_weighted::<_, N>(|_| 1.0);
    println!("∫ 1/√(1−x²) dx ≈ {pi_estimate}");

    let j0_estimate: f64 = gauss_chebyshev::integrate_weighted::<_, N>(f64::cos);
    println!("∫ cos(x)/√(1−x²) dx ≈ {j0_estimate}");
    println!("(reference: π·J₀(1) ≈ {})", core::f64::consts::PI * 0.7651976865579666);
}
