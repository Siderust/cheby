//! Accuracy checks for quadrature and spectral differentiation.
#![allow(clippy::needless_range_loop)]
#![cfg(all(feature = "quadrature", feature = "spectral"))]

use approx::assert_abs_diff_eq;
use cheby::core::{Domain, NodeKind};
use cheby::quadrature::{clenshaw_curtis, gauss_chebyshev};
use cheby::spectral::{chebyshev_differentiation_matrix, collocation_points};

#[test]
fn clenshaw_curtis_integrates_polynomial() {
    // ∫_{0}^{2} x^3 dx = 4
    let v = clenshaw_curtis::integrate::<f64, f64, 9>(Domain::new(0.0, 2.0), |x| x * x * x);
    assert_abs_diff_eq!(v, 4.0, epsilon = 1e-9);
}

#[test]
fn clenshaw_curtis_integrates_smooth_function() {
    // ∫_{0}^{π} sin x dx = 2
    let v = clenshaw_curtis::integrate::<f64, f64, 33>(
        Domain::new(0.0, core::f64::consts::PI),
        f64::sin,
    );
    assert_abs_diff_eq!(v, 2.0, epsilon = 1e-3);
}

#[test]
fn gauss_chebyshev_recovers_known_weighted_integrals() {
    // ∫_{-1}^{1} 1 / sqrt(1 - x^2) dx = π
    let v = gauss_chebyshev::integrate_weighted::<f64, 32>(|_| 1.0);
    assert_abs_diff_eq!(v, core::f64::consts::PI, epsilon = 1e-12);

    // ∫_{-1}^{1} cos(x) / sqrt(1 - x^2) dx = π * J_0(1)
    let j0_1 = 0.7651976865579665_f64;
    let v = gauss_chebyshev::integrate_weighted::<f64, 64>(f64::cos);
    assert_abs_diff_eq!(v, core::f64::consts::PI * j0_1, epsilon = 1e-10);

    // Polynomial: x^2; closed form is π/2.
    let v = gauss_chebyshev::integrate_weighted::<f64, 8>(|x| x * x);
    assert_abs_diff_eq!(v, core::f64::consts::PI / 2.0, epsilon = 1e-12);
}

#[test]
fn spectral_differentiation_matches_analytic_derivative() {
    // For f(x) = sin(x), differentiate at Lobatto nodes: D * f ≈ cos(x).
    const N: usize = 17;
    let nodes: [f64; N] = collocation_points::<N>();
    let f: [f64; N] = std::array::from_fn(|k| nodes[k].sin());
    let d = chebyshev_differentiation_matrix(N);

    for i in 0..N {
        let mut acc = 0.0;
        for j in 0..N {
            acc += d.get(i, j) * f[j];
        }
        let exact = nodes[i].cos();
        assert_abs_diff_eq!(acc, exact, epsilon = 1e-9);
    }
}

#[test]
fn spectral_differentiation_zero_on_constant() {
    const N: usize = 9;
    let f = [3.5_f64; N];
    let d = chebyshev_differentiation_matrix(N);
    for i in 0..N {
        let mut acc = 0.0;
        for j in 0..N {
            acc += d.get(i, j) * f[j];
        }
        assert_abs_diff_eq!(acc, 0.0, epsilon = 1e-12);
    }
}

#[test]
fn lobatto_nodes_are_endpoints_inclusive() {
    let nodes: [f64; 5] = cheby::core::nodes::nodes(NodeKind::GaussLobatto);
    assert_abs_diff_eq!(nodes[0], 1.0, epsilon = 1e-15);
    assert_abs_diff_eq!(nodes[4], -1.0, epsilon = 1e-15);
}
