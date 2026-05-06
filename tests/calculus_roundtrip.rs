//! Calculus round-trip and integral correctness.

use approx::assert_abs_diff_eq;
use cheby::approx::fit::fit_from_fn_on;
use cheby::core::{ChebySeries, Domain};

#[test]
fn integral_then_derivative_recovers_original() {
    // Pad with trailing zeros so the antiderivative stays representable.
    let series = ChebySeries::new([1.5_f64, -0.4, 0.2, 0.0, 0.0]);
    let recovered = series.integral(0.0).derivative();
    for tau in [-0.9, -0.3, 0.0, 0.4, 0.85] {
        assert_abs_diff_eq!(
            recovered.evaluate(tau),
            series.evaluate(tau),
            epsilon = 1e-12
        );
    }
}

#[test]
fn definite_integral_matches_analytic() {
    // ∫_{0}^{x} sin(t) dt = 1 - cos(x), and we test for x in (0, π/2).
    let domain = Domain::new(0.0, core::f64::consts::PI);
    let series = fit_from_fn_on::<f64, f64, 33>(domain, f64::sin);
    for &x in &[0.5_f64, 1.0, 1.5, 2.5, 3.0] {
        let v = series.evaluate_integral_from_start(x).unwrap();
        let exact = 1.0 - x.cos();
        assert_abs_diff_eq!(v, exact, epsilon = 1e-9);
    }
}

#[test]
fn definite_integral_rejects_out_of_domain_point() {
    let domain = Domain::new(0.0, 1.0);
    let series = ChebySeries::new([1.0_f64, 0.0]);
    let series = cheby::ChebySeriesOn::new(domain, series);
    assert_eq!(
        series.evaluate_integral_from_start(1.5).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain
    );
}

#[test]
fn typed_integral_velocity_to_position() {
    use qtty::{Kilometer, Second};
    type Velocity = qtty::velocity::Velocity<qtty::unit::Kilometer, qtty::unit::Second>;

    let domain = Domain::new(Second::new(0.0), Second::new(10.0));
    // Constant velocity = 2 km/s ⇒ position(t) = 2t.
    let series = cheby::ChebySeriesOn::new(
        domain,
        ChebySeries::new([Velocity::new(2.0), Velocity::new(0.0)]),
    );
    let pos: Kilometer = series
        .evaluate_integral_from_start(Second::new(3.0))
        .unwrap();
    assert_abs_diff_eq!(pos.value(), 6.0, epsilon = 1e-12);
}

#[test]
fn cached_definite_integral_matches_one_shot() {
    let domain = Domain::new(0.0, core::f64::consts::PI);
    let series = fit_from_fn_on::<f64, f64, 17>(domain, f64::sin);
    let cached = series.definite_integral_from_start();
    for &x in &[0.0_f64, 0.3, 0.7, 1.5, 2.4, 3.0, core::f64::consts::PI] {
        let one_shot = series.evaluate_integral_from_start(x).unwrap();
        let from_cache = cached.evaluate(x).unwrap();
        assert_abs_diff_eq!(from_cache, one_shot, epsilon = 1e-14);
        let exact = 1.0 - x.cos();
        assert_abs_diff_eq!(from_cache, exact, epsilon = 1e-9);
    }
}

#[test]
fn cached_definite_integral_rejects_out_of_domain() {
    let domain = Domain::new(0.0, 1.0);
    let series = cheby::ChebySeriesOn::new(domain, ChebySeries::new([1.0_f64, 0.0]));
    let cached = series.definite_integral_from_start();
    assert_eq!(
        cached.evaluate(1.5).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain
    );
}
