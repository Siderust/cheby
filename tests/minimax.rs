//! Real Remez exchange — convergence and failure modes.
#![cfg(feature = "minimax")]

use approx::assert_abs_diff_eq;
use cheby::approx::minimax::{remez, RemezOptions};
use cheby::core::{ChebyError, Domain};

#[test]
fn remez_converges_for_smooth_target() {
    let result = remez(Domain::new(-1.0, 1.0), 7, f64::exp, RemezOptions::default()).unwrap();
    assert!(result.converged);
    assert_eq!(result.alternations.len(), 9); // degree + 2
    assert!(result.max_error < 1e-6);
    // Equioscillation: dense max error should be very close to leveled error.
    assert_abs_diff_eq!(result.max_error, result.leveled_error, epsilon = 1e-6);
}

#[test]
fn remez_alternations_are_monotone() {
    let result = remez(Domain::new(-1.0, 1.0), 5, f64::cos, RemezOptions::default()).unwrap();
    for w in result.alternations.windows(2) {
        assert!(w[0] < w[1]);
    }
}

#[test]
fn remez_rejects_zero_degree() {
    let err = remez(Domain::new(-1.0, 1.0), 0, f64::sin, RemezOptions::default()).unwrap_err();
    assert_eq!(err, ChebyError::InvalidDegree);
}

#[test]
fn remez_typed_domain_evaluates() {
    use qtty::Second;
    let result = remez(
        Domain::new(Second::new(0.0), Second::new(1.0)),
        5,
        |t: Second| t.value().exp(),
        RemezOptions::default(),
    )
    .unwrap();
    assert!(result.converged);
    let v = result.series_on.evaluate(Second::new(0.5)).unwrap();
    assert_abs_diff_eq!(v, 0.5_f64.exp(), epsilon = 1e-4);
}

#[test]
fn remez_does_not_converge_with_tiny_budget() {
    let opts = RemezOptions {
        max_iterations: 1,
        tolerance: 1e-30,
        grid_size: 64,
    };
    // exp on [-1,1] with degree 3 won't reach 1e-30 in a single iteration.
    let res = remez(Domain::new(-1.0, 1.0), 3, f64::exp, opts);
    assert!(matches!(res, Err(ChebyError::RemezDidNotConverge)));
}

#[test]
fn remez_converges_better_than_root_fit() {
    let opts = RemezOptions::default();
    let degree = 6;
    let domain = Domain::new(-1.0, 1.0);
    let f = |x: f64| x.exp();

    let remez_res = remez(domain, degree, f, opts).unwrap();

    // Root-sample fit at degree+1 nodes (interpolating polynomial).
    let root_series = cheby::approx::fit::fit_from_fn::<f64, 7>(f);
    let root_max_err = (0..1024)
        .map(|i| {
            let tau = -1.0 + 2.0 * i as f64 / 1023.0;
            (f(tau) - root_series.evaluate(tau)).abs()
        })
        .fold(0.0_f64, f64::max);

    // Minimax must be at least as good; almost always strictly better.
    assert!(
        remez_res.max_error <= root_max_err + 1e-12,
        "remez={} root={}",
        remez_res.max_error,
        root_max_err
    );
}
