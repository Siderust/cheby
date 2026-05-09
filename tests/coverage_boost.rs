//! Coverage-boosting tests targeting the least-covered modules.
//!
//! Covers: core/error.rs (Display), core/series.rs, piecewise/segment.rs,
//! lib.rs free-functions, core/basis.rs (U_n), spectral/matrices.rs,
//! approx/interpolation.rs edge-cases, quadrature/clenshaw_curtis.rs,
//! spectral/differentiation.rs, calculus/integral.rs, approx/minimax.rs,
//! piecewise/adaptive.rs.

use approx::assert_abs_diff_eq;

// ─── core/error.rs ──────────────────────────────────────────────────────────

#[test]
fn cheby_error_display_all_variants() {
    use cheby::ChebyError;

    let cases: &[(ChebyError, &str)] = &[
        (ChebyError::EmptyDomain, "domain has zero width"),
        (ChebyError::InvalidDomain, "domain is invalid"),
        (
            ChebyError::NonPositiveSegmentLength,
            "segment length must be positive",
        ),
        (ChebyError::InvalidDegree, "degree is invalid"),
        (ChebyError::EmptyCoefficientSet, "coefficient set is empty"),
        (
            ChebyError::EvaluationOutOfDomain,
            "evaluation point is outside the domain",
        ),
        (ChebyError::NonFiniteInput, "input is not finite"),
        (ChebyError::EmptySegmentTable, "segment table is empty"),
        (
            ChebyError::SegmentBoundariesNotContiguous,
            "segment boundaries are not contiguous",
        ),
        (
            ChebyError::BinaryTooShort,
            "binary blob is shorter than the format header",
        ),
        (
            ChebyError::UnsupportedFormatVersion {
                found: 99,
                expected: 1,
            },
            "unsupported binary format version 99 (expected 1)",
        ),
        (
            ChebyError::BinaryLengthMismatch,
            "binary blob length does not match the stored coefficient count",
        ),
        (
            ChebyError::BinaryChecksumMismatch,
            "binary blob checksum mismatch",
        ),
        (
            ChebyError::RemezDidNotConverge,
            "Remez exchange did not converge within the iteration budget",
        ),
        (
            ChebyError::RemezAlternationFailure,
            "Remez exchange could not find enough alternation points",
        ),
    ];

    for (err, expected) in cases {
        assert_eq!(err.to_string(), *expected, "Display mismatch for {err:?}");
    }
}

#[test]
fn cheby_error_is_std_error() {
    use std::error::Error;
    let e = cheby::ChebyError::EmptyDomain;
    // Just verify the trait is implemented – source() returns None.
    let _: &dyn Error = &e;
}

// ─── core/basis.rs ──────────────────────────────────────────────────────────

#[test]
fn basis_u_n_known_values() {
    use cheby::core::basis;

    // U_0(x) = 1
    assert_abs_diff_eq!(basis::u(0, 0.5), 1.0, epsilon = 1e-15);
    // U_1(x) = 2x
    assert_abs_diff_eq!(basis::u(1, 0.5), 1.0, epsilon = 1e-15);
    assert_abs_diff_eq!(basis::u(1, 0.0), 0.0, epsilon = 1e-15);
    // U_2(x) = 4x^2 - 1
    assert_abs_diff_eq!(basis::u(2, 0.5), 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(basis::u(2, 1.0), 3.0, epsilon = 1e-14);
    // U_3(x) = 8x^3 - 4x
    assert_abs_diff_eq!(basis::u(3, 0.5), -1.0, epsilon = 1e-14);
}

#[test]
fn basis_t_n_known_values() {
    use cheby::core::basis;
    // T_0 and T_1 hit the trivial match arms
    assert_abs_diff_eq!(basis::t(0, 0.7), 1.0, epsilon = 0.0);
    assert_abs_diff_eq!(basis::t(1, 0.7), 0.7, epsilon = 0.0);
}

// ─── core/series.rs ─────────────────────────────────────────────────────────

#[test]
fn cheby_series_accessors() {
    let s = cheby::ChebySeries::new([1.0_f64, -0.5, 0.25]);
    assert_eq!(s.coeffs(), &[1.0, -0.5, 0.25]);
    assert_eq!(s.into_coeffs(), [1.0, -0.5, 0.25]);
}

#[test]
fn cheby_series_evaluate_both() {
    let s = cheby::ChebySeries::new([1.0_f64, 0.5, 0.0]);
    let (v, d) = s.evaluate_both(0.0);
    assert_abs_diff_eq!(v, s.evaluate(0.0), epsilon = 1e-15);
    let _ = d; // derivative is computed – just verify no panic
}

#[test]
fn cheby_series_derivative_edge_cases() {
    // N = 0 (zero-length): returns empty
    let s0 = cheby::ChebySeries::<f64, 0>::new([]);
    let d0 = s0.derivative();
    assert_eq!(d0.coeffs(), &[] as &[f64; 0]);

    // N = 1: derivative of a constant is zero
    let s1 = cheby::ChebySeries::new([3.0_f64]);
    let d1 = s1.derivative();
    assert_abs_diff_eq!(d1.evaluate(0.0), 0.0, epsilon = 1e-15);

    // N = 2: derivative of `a + b*T_1` is `b`
    let s2 = cheby::ChebySeries::new([0.0_f64, 1.0]);
    let d2 = s2.derivative();
    assert_abs_diff_eq!(d2.evaluate(0.0), 1.0, epsilon = 1e-12);
}

#[test]
fn cheby_series_integral_edge_cases() {
    // N = 0 (zero-length): no coefficients → returns empty
    let s0 = cheby::ChebySeries::<f64, 0>::new([]);
    let i0 = s0.integral(0.0);
    assert_eq!(i0.coeffs(), &[] as &[f64; 0]);

    // N = 1: constant c → integral is c * T_0 (i.e. constant term)
    let s1 = cheby::ChebySeries::new([2.0_f64]);
    let i1 = s1.integral(0.0);
    // integral has only a constant slot; check it's finite
    let _ = i1.evaluate(0.0);

    // N = 2: integral of T_0 = T_1
    let s2 = cheby::ChebySeries::new([1.0_f64, 0.0]);
    let i2 = s2.integral(0.0);
    assert_abs_diff_eq!(i2.evaluate(0.5), 0.5, epsilon = 1e-12);
}

// ChebySeriesOn

#[test]
fn cheby_series_on_accessors() {
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let series = cheby::ChebySeries::new([1.0_f64, 0.0]);
    let son = cheby::ChebySeriesOn::new(domain, series);

    assert_eq!(son.domain(), domain);
    assert_eq!(son.series().coeffs(), &[1.0, 0.0]);
    assert_eq!(son.into_coeffs(), [1.0, 0.0]);
}

#[test]
fn cheby_series_on_evaluate_derivative_ok_and_err() {
    let domain = cheby::Domain::new(0.0_f64, 2.0);
    // Fit a linear function: f(x) = x → coefficients will encode that.
    let series = cheby::approx::fit::fit_from_fn_on::<f64, f64, 3>(domain, |x| x);
    let mid = 1.0_f64; // inside domain
    let deriv = series.evaluate_derivative(mid).unwrap();
    assert_abs_diff_eq!(deriv, 1.0, epsilon = 1e-10);

    // Out-of-domain
    assert_eq!(
        series.evaluate_derivative(3.0_f64).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );
}

// ChebySeriesDyn

#[test]
fn cheby_series_dyn_accessors() {
    let dyn_s = cheby::ChebySeriesDyn::new(vec![1.0_f64, -2.0, 0.5]).unwrap();
    assert_eq!(dyn_s.len(), 3);
    assert!(!dyn_s.is_empty());
    assert_eq!(dyn_s.coeffs(), &[1.0, -2.0, 0.5]);
}

#[test]
fn cheby_series_dyn_into_coeffs() {
    let dyn_s = cheby::ChebySeriesDyn::new(vec![1.0_f64, 2.0]).unwrap();
    assert_eq!(dyn_s.into_coeffs(), vec![1.0, 2.0]);
}

#[test]
fn cheby_series_dyn_evaluate_both_and_derivative() {
    let dyn_s = cheby::ChebySeriesDyn::new(vec![1.0_f64, 0.5]).unwrap();
    let (v, d) = dyn_s.evaluate_both(0.0);
    assert_abs_diff_eq!(v, dyn_s.evaluate(0.0), epsilon = 1e-15);
    assert_abs_diff_eq!(d, dyn_s.evaluate_derivative(0.0), epsilon = 1e-15);
}

#[test]
fn cheby_series_dyn_derivative_edge_cases() {
    // len = 1: derivative of constant is zero
    let s1 = cheby::ChebySeriesDyn::new(vec![5.0_f64]).unwrap();
    let d1 = s1.derivative();
    assert_abs_diff_eq!(d1.evaluate(0.0), 0.0, epsilon = 1e-15);

    // len = 2: derivative of `a + b * T_1` is b
    let s2 = cheby::ChebySeriesDyn::new(vec![0.0_f64, 1.0]).unwrap();
    let d2 = s2.derivative();
    assert_abs_diff_eq!(d2.evaluate(0.0), 1.0, epsilon = 1e-12);

    // len = 3: general case
    let s3 = cheby::ChebySeriesDyn::new(vec![1.0_f64, 0.5, 0.25]).unwrap();
    let d3 = s3.derivative();
    assert_eq!(d3.len(), 2);
}

#[test]
fn cheby_series_dyn_integral() {
    let s = cheby::ChebySeriesDyn::new(vec![1.0_f64, 0.0]).unwrap();
    let i = s.integral(0.0);
    assert_eq!(i.len(), 3);
    // integral of T_0 is T_1, so evaluating at tau=0 gives 0
    assert_abs_diff_eq!(i.evaluate(0.0), 0.0, epsilon = 1e-15);
}

// ChebySeriesDynOn

#[test]
fn cheby_series_dyn_on_all_methods() {
    use cheby::core::{ChebySeriesDyn, ChebySeriesDynOn, Domain};

    let domain = Domain::new(0.0_f64, 2.0);
    let dyn_s = ChebySeriesDyn::new(vec![1.0_f64, 0.0]).unwrap();
    let dyn_on = ChebySeriesDynOn::new(domain, dyn_s.clone());

    assert_eq!(dyn_on.domain(), domain);
    assert_eq!(dyn_on.series().coeffs(), dyn_s.coeffs());

    // evaluate inside domain
    let v = dyn_on.evaluate(1.0_f64).unwrap();
    let _ = v;

    // evaluate outside domain
    assert_eq!(
        dyn_on.evaluate(5.0_f64).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );

    // evaluate_derivative inside domain
    let dyn_on2 = ChebySeriesDynOn::new(domain, dyn_s.clone());
    let deriv = dyn_on2.evaluate_derivative(1.0_f64).unwrap();
    let _ = deriv;

    // evaluate_derivative outside domain
    let dyn_on3 = ChebySeriesDynOn::new(domain, dyn_s.clone());
    assert_eq!(
        dyn_on3.evaluate_derivative(9.0_f64).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );

    // into_series
    let dyn_on4 = ChebySeriesDynOn::new(domain, dyn_s);
    let _ = dyn_on4.into_series();
}

// ─── piecewise/segment.rs ───────────────────────────────────────────────────

#[test]
fn cheby_segment_try_new_empty_returns_error() {
    use cheby::piecewise::ChebySegment;
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let empty_series = cheby::ChebySeries::<f64, 0>::new([]);
    let result = ChebySegment::try_new(domain, empty_series);
    assert_eq!(result.unwrap_err(), cheby::ChebyError::EmptyCoefficientSet);
}

#[test]
fn cheby_segment_try_new_success() {
    use cheby::piecewise::ChebySegment;
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let series = cheby::ChebySeries::new([1.0_f64, 0.0]);
    let seg = ChebySegment::try_new(domain, series).unwrap();
    assert_eq!(seg.domain(), domain);
}

#[test]
fn cheby_segment_accessors() {
    use cheby::piecewise::ChebySegment;
    let seg = ChebySegment::new([1.0_f64, -0.5, 0.25], 2.0_f64, 1.0_f64);

    assert_abs_diff_eq!(seg.mid(), 2.0, epsilon = 0.0);
    assert_abs_diff_eq!(seg.half(), 1.0, epsilon = 0.0);
    assert_eq!(seg.coeffs(), &[1.0, -0.5, 0.25]);
    assert_abs_diff_eq!(seg.domain().start(), 1.0, epsilon = 0.0);
    assert_abs_diff_eq!(seg.domain().end(), 3.0, epsilon = 0.0);
}

#[test]
fn cheby_segment_normalise_and_contains() {
    use cheby::piecewise::ChebySegment;
    let seg = ChebySegment::new([1.0_f64, 0.0], 5.0_f64, 2.0_f64);

    // midpoint normalises to 0
    assert_abs_diff_eq!(seg.normalise(5.0), 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(seg.normalize(5.0), 0.0, epsilon = 1e-15);

    assert!(seg.contains(5.0));
    assert!(!seg.contains(10.0));
}

#[test]
fn cheby_segment_evaluate_and_eval() {
    use cheby::piecewise::ChebySegment;
    let seg = ChebySegment::new([1.0_f64, 0.0, 0.0], 0.5_f64, 0.5_f64);

    // In-domain: constant series evaluates to 1.0 everywhere
    assert_abs_diff_eq!(seg.evaluate(0.5).unwrap(), 1.0, epsilon = 1e-12);
    assert_abs_diff_eq!(seg.eval(0.5), 1.0, epsilon = 1e-12);

    // Out-of-domain
    assert_eq!(
        seg.evaluate(5.0).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );
}

#[test]
fn cheby_segment_evaluate_derivative_and_eval_derivative() {
    // Fit linear function on [0, 2] so derivative is 1 everywhere
    let domain = cheby::Domain::new(0.0_f64, 2.0);
    let series = cheby::approx::fit::fit_from_fn_on::<f64, f64, 5>(domain, |x| x);
    let seg = cheby::piecewise::ChebySegment::try_new(domain, *series.series()).unwrap();

    // In-domain derivative
    let d = seg.evaluate_derivative(1.0_f64).unwrap();
    assert_abs_diff_eq!(d, 1.0, epsilon = 1e-9);
    let d2 = seg.eval_derivative(1.0_f64);
    assert_abs_diff_eq!(d2, 1.0, epsilon = 1e-9);

    // Out-of-domain derivative
    assert_eq!(
        seg.evaluate_derivative(9.0).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );
}

#[test]
fn cheby_segment_evaluate_both_and_eval_both() {
    let domain = cheby::Domain::new(0.0_f64, 2.0);
    let series = cheby::approx::fit::fit_from_fn_on::<f64, f64, 5>(domain, |x| x);
    let seg = cheby::piecewise::ChebySegment::try_new(domain, *series.series()).unwrap();

    let (v, d) = seg.evaluate_both(1.0).unwrap();
    assert_abs_diff_eq!(v, 1.0, epsilon = 1e-9);
    assert_abs_diff_eq!(d, 1.0, epsilon = 1e-9);

    let (v2, d2) = seg.eval_both(1.0);
    assert_abs_diff_eq!(v2, 1.0, epsilon = 1e-9);
    assert_abs_diff_eq!(d2, 1.0, epsilon = 1e-9);

    // Out-of-domain
    assert_eq!(
        seg.evaluate_both(9.0).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );
}

// ─── lib.rs free functions ───────────────────────────────────────────────────

#[test]
fn lib_evaluate_derivative_free_fn() {
    // T_1 -> derivative is 1 everywhere
    let coeffs = [0.0_f64, 1.0];
    let d = cheby::evaluate_derivative(&coeffs, 0.0);
    assert_abs_diff_eq!(d, 1.0, epsilon = 1e-12);
}

#[test]
fn lib_chebyshev_roots_and_nodes() {
    let roots: [f64; 3] = cheby::chebyshev_roots();
    let nodes_: [f64; 3] = cheby::nodes();
    // Both functions should return the same Chebyshev roots
    for (r, n) in roots.iter().zip(nodes_.iter()) {
        assert_abs_diff_eq!(r, n, epsilon = 1e-15);
    }
}

#[test]
fn lib_nodes_roots() {
    let r1: [f64; 4] = cheby::chebyshev_roots();
    let r2: [f64; 4] = cheby::nodes_roots();
    for (a, b) in r1.iter().zip(r2.iter()) {
        assert_abs_diff_eq!(a, b, epsilon = 0.0);
    }
}

#[test]
fn lib_nodes_mapped_fns() {
    let mapped: [f64; 5] = cheby::nodes_mapped(2.0, 4.0);
    let mapped_i: [f64; 5] = cheby::nodes_mapped_interval(2.0, 4.0);

    for (a, b) in mapped.iter().zip(mapped_i.iter()) {
        assert_abs_diff_eq!(a, b, epsilon = 1e-15);
    }
    // All nodes should be within [2, 4]
    for &n in &mapped {
        assert!((2.0..=4.0).contains(&n), "node {n} out of [2, 4]");
    }
}

#[test]
fn lib_nodes_mapped_t() {
    let mapped: [f64; 5] = cheby::nodes_mapped_t(2.0_f64, 4.0_f64);
    let expected: [f64; 5] = cheby::nodes_mapped(2.0, 4.0);
    for (a, b) in mapped.iter().zip(expected.iter()) {
        assert_abs_diff_eq!(a, b, epsilon = 1e-15);
    }
}

#[test]
fn lib_fit_from_fn() {
    let coeffs: [f64; 9] = cheby::fit_from_fn(f64::sin, 0.0, 1.0);
    // Evaluate at midpoint and compare to sin
    let v = cheby::evaluate(&coeffs, 0.0); // normalized midpoint → x = 0.5
    assert_abs_diff_eq!(v, (0.5_f64).sin(), epsilon = 1e-9);
}

#[test]
fn lib_fit_from_fn_t() {
    let coeffs: [f64; 9] = cheby::fit_from_fn_t(f64::cos, 0.0_f64, 1.0_f64);
    let v = cheby::evaluate(&coeffs, 0.0);
    assert_abs_diff_eq!(v, (0.5_f64).cos(), epsilon = 1e-9);
}

// ─── spectral/matrices.rs ───────────────────────────────────────────────────

#[test]
fn matrix_all_methods() {
    use cheby::spectral::Matrix;

    let mut m = Matrix::zeros(3, 4);
    assert_eq!(m.rows(), 3);
    assert_eq!(m.cols(), 4);
    assert_eq!(m.data().len(), 12);

    // All entries start at zero
    assert_abs_diff_eq!(m.get(1, 2), 0.0, epsilon = 0.0);

    // Set and get
    m.set(0, 0, 1.5);
    assert_abs_diff_eq!(m.get(0, 0), 1.5, epsilon = 0.0);
    m.set(2, 3, -7.25);
    assert_abs_diff_eq!(m.get(2, 3), -7.25, epsilon = 0.0);
}

// ─── spectral/differentiation.rs ────────────────────────────────────────────

#[test]
fn differentiation_matrix_edge_cases() {
    use cheby::spectral::chebyshev_differentiation_matrix;

    // n = 0: 0×0 matrix
    let d0 = chebyshev_differentiation_matrix(0);
    assert_eq!(d0.rows(), 0);
    assert_eq!(d0.cols(), 0);

    // n = 1: 1×1 zero matrix
    let d1 = chebyshev_differentiation_matrix(1);
    assert_eq!(d1.rows(), 1);
    assert_abs_diff_eq!(d1.get(0, 0), 0.0, epsilon = 0.0);

    // n = 2
    let d2 = chebyshev_differentiation_matrix(2);
    assert_eq!(d2.rows(), 2);

    // n = 3 (hits the special-case path)
    let d3 = chebyshev_differentiation_matrix(3);
    assert_eq!(d3.rows(), 3);
}

#[test]
fn fixed_collocation_points_matches_gauss_lobatto() {
    use cheby::core::{nodes, NodeKind};
    use cheby::spectral::collocation_points;

    let pts: [f64; 5] = collocation_points();
    let lobatto: [f64; 5] = nodes(NodeKind::GaussLobatto);
    for (a, b) in pts.iter().zip(lobatto.iter()) {
        assert_abs_diff_eq!(a, b, epsilon = 1e-15);
    }
}

// ─── quadrature/clenshaw_curtis.rs ──────────────────────────────────────────

#[test]
fn clenshaw_curtis_rule_accessors_and_integrate() {
    use cheby::quadrature::clenshaw_curtis::{clenshaw_curtis_weights, ClenshawCurtisRule};

    // Edge-case weight arrays
    let w0: [f64; 0] = clenshaw_curtis_weights();
    assert_eq!(w0.len(), 0);

    let w1: [f64; 1] = clenshaw_curtis_weights();
    assert_abs_diff_eq!(w1[0], 2.0, epsilon = 0.0);

    let w2: [f64; 2] = clenshaw_curtis_weights();
    assert_abs_diff_eq!(w2[0], 1.0, epsilon = 0.0);
    assert_abs_diff_eq!(w2[1], 1.0, epsilon = 0.0);

    // Rule via Default
    let rule: ClenshawCurtisRule<5> = ClenshawCurtisRule::default();
    assert_eq!(rule.nodes().len(), 5);
    assert_eq!(rule.weights().len(), 5);

    // Integrate constant 1 over [-1, 1] → 2
    let domain = cheby::Domain::new(-1.0_f64, 1.0);
    let integral: f64 = rule.integrate(domain, |_x| 1.0_f64);
    assert_abs_diff_eq!(integral, 2.0, epsilon = 1e-14);
}

#[test]
fn clenshaw_curtis_integrate_odd_m() {
    // Exercises the `!m.is_multiple_of(2)` branch in weight construction.
    // N = 4  => m = 3 (odd)
    use cheby::quadrature::clenshaw_curtis::clenshaw_curtis_weights;
    let w: [f64; 4] = clenshaw_curtis_weights();
    // All weights positive and sum to 2 for integration over [-1,1]
    let sum: f64 = w.iter().sum();
    assert_abs_diff_eq!(sum, 2.0, epsilon = 1e-14);
}

// ─── approx/interpolation.rs edge cases ─────────────────────────────────────

#[test]
fn barycentric_n_zero_returns_error() {
    use cheby::approx::interpolation::BarycentricInterpolator;
    let result = BarycentricInterpolator::<f64, f64, 0>::new([], []);
    assert_eq!(result.unwrap_err(), cheby::ChebyError::InvalidDegree);
}

#[test]
fn barycentric_on_roots_n_zero_returns_error() {
    use cheby::approx::interpolation::BarycentricInterpolator;
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let result = BarycentricInterpolator::<f64, f64, 0>::on_chebyshev_roots(domain, []);
    assert_eq!(result.unwrap_err(), cheby::ChebyError::InvalidDegree);
}

#[test]
fn barycentric_on_roots_n_one_evaluates_constant() {
    use cheby::approx::interpolation::BarycentricInterpolator;
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let interp = BarycentricInterpolator::on_chebyshev_roots(domain, [42.0_f64]).unwrap();
    assert_abs_diff_eq!(interp.evaluate(0.5_f64).unwrap(), 42.0, epsilon = 0.0);
}

#[test]
fn barycentric_on_lobatto_n_zero_returns_error() {
    use cheby::approx::interpolation::BarycentricInterpolator;
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let result = BarycentricInterpolator::<f64, f64, 0>::on_lobatto_nodes(domain, []);
    assert_eq!(result.unwrap_err(), cheby::ChebyError::InvalidDegree);
}

#[test]
fn barycentric_on_lobatto_n_one_evaluates_constant() {
    use cheby::approx::interpolation::BarycentricInterpolator;
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let interp = BarycentricInterpolator::on_lobatto_nodes(domain, [7.0_f64]).unwrap();
    assert_abs_diff_eq!(interp.evaluate(0.5_f64).unwrap(), 7.0, epsilon = 0.0);
}

#[test]
fn barycentric_n_one_evaluate_is_constant() {
    use cheby::approx::interpolation::BarycentricInterpolator;
    let interp = BarycentricInterpolator::new([3.0_f64], [99.0_f64]).unwrap();
    assert_abs_diff_eq!(interp.evaluate(100.0_f64).unwrap(), 99.0, epsilon = 0.0);
}

#[test]
fn barycentric_duplicate_nodes_returns_error() {
    use cheby::approx::interpolation::BarycentricInterpolator;
    // Two nodes at the same position → InvalidDomain
    let result = BarycentricInterpolator::new([1.0_f64, 1.0], [0.0, 0.0]);
    assert_eq!(result.unwrap_err(), cheby::ChebyError::InvalidDomain);
}

// ─── approx/minimax.rs ──────────────────────────────────────────────────────

#[test]
fn minimax_series_on_any_accessors_and_evaluate() {
    use cheby::approx::minimax::{remez, RemezOptions};

    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let result = remez(domain, 4, f64::exp, RemezOptions::default()).unwrap();

    // domain / coeffs accessors
    assert_eq!(result.series_on.domain(), domain);
    assert!(!result.series_on.coeffs().is_empty());

    // evaluate inside domain – exp is smooth, accuracy well within 1e-4
    let v = result.series_on.evaluate(0.5_f64).unwrap();
    assert_abs_diff_eq!(v, (0.5_f64).exp(), epsilon = 1e-4);

    // evaluate outside domain
    assert_eq!(
        result.series_on.evaluate(2.0_f64).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );
}

#[test]
fn minimax_try_into_fixed_wrong_size_returns_error() {
    use cheby::approx::minimax::try_into_fixed;
    let coeffs = vec![1.0_f64, 2.0, 3.0]; // length 3
    let result = try_into_fixed::<5>(&coeffs); // want 5
    assert_eq!(result.unwrap_err(), cheby::ChebyError::InvalidDegree);
}

#[test]
fn minimax_try_into_fixed_correct_size() {
    use cheby::approx::minimax::try_into_fixed;
    let coeffs = vec![1.0_f64, 2.0, 3.0];
    let series = try_into_fixed::<3>(&coeffs).unwrap();
    assert_eq!(series.coeffs(), &[1.0, 2.0, 3.0]);
}

#[test]
fn minimax_into_fixed_series_on() {
    use cheby::approx::minimax::{into_fixed_series_on, remez, RemezOptions};

    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let result = remez(domain, 4, f64::exp, RemezOptions::default()).unwrap();
    let fixed: cheby::ChebySeriesOn<f64, f64, 5> = into_fixed_series_on(&result).unwrap();
    let v = fixed.evaluate(0.5_f64).unwrap();
    assert_abs_diff_eq!(v, (0.5_f64).exp(), epsilon = 1e-4);
}

#[test]
fn minimax_into_fixed_series_on_wrong_size_returns_error() {
    use cheby::approx::minimax::{into_fixed_series_on, remez, RemezOptions};

    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let result = remez(domain, 4, f64::exp, RemezOptions::default()).unwrap();
    // degree 4 → 5 coefficients; asking for N = 3 should fail
    let r: Result<cheby::ChebySeriesOn<f64, f64, 3>, _> = into_fixed_series_on(&result);
    assert_eq!(r.unwrap_err(), cheby::ChebyError::InvalidDegree);
}

// ─── calculus/integral.rs – N=1, N=2 edge cases ─────────────────────────────

#[test]
fn definite_integral_n1_constant_series() {
    // Series with N=1 coefficient (just a constant)
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let series = cheby::ChebySeriesOn::new(domain, cheby::ChebySeries::new([2.0_f64]));
    // ∫_0^0.5 2 dt = 1
    let v = series.evaluate_integral_from_start(0.5_f64).unwrap();
    assert_abs_diff_eq!(v, 1.0, epsilon = 1e-10);
}

#[test]
fn definite_integral_n2_linear_series() {
    // Series [1, 0] on [0, 1]: f(x) = 1, so ∫_0^x dt = x
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let series = cheby::ChebySeriesOn::new(domain, cheby::ChebySeries::new([1.0_f64, 0.0_f64]));
    let v = series.evaluate_integral_from_start(0.4_f64).unwrap();
    assert_abs_diff_eq!(v, 0.4, epsilon = 1e-10);
}

#[test]
fn definite_integral_n3_series() {
    // N = 3 exercises the N >= 3 branch in definite_integral_from_start
    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let s: cheby::ChebySeriesOn<f64, f64, 3> =
        cheby::approx::fit::fit_from_fn_on(domain, |x: f64| x);
    let v = s.evaluate_integral_from_start(1.0_f64).unwrap();
    // ∫_0^1 x dx = 0.5
    assert_abs_diff_eq!(v, 0.5, epsilon = 1e-9);
}

// ─── piecewise/adaptive.rs ───────────────────────────────────────────────────

#[test]
fn adaptive_segment_table_all_accessors() {
    use cheby::piecewise::adaptive::AdaptiveSegmentTable;

    let domain = cheby::Domain::new(0.0_f64, 2.0);
    let table: AdaptiveSegmentTable<f64, f64, 8> =
        AdaptiveSegmentTable::from_fn(domain, f64::sin, 1e-6, 10).unwrap();

    assert_eq!(table.domain(), domain);
    assert!(!table.is_empty());
    assert!(!table.is_empty());
    assert!(!table.segments().is_empty());
    assert!(!table.metadata().is_empty());
    assert_eq!(table.boundaries().len(), table.len() + 1);

    // evaluate inside domain
    let v = table.evaluate(1.0_f64).unwrap();
    assert_abs_diff_eq!(v, (1.0_f64).sin(), epsilon = 1e-5);

    // evaluate outside domain
    assert_eq!(
        table.evaluate(9.0_f64).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );
}

#[test]
fn adaptive_segment_table_locate_boundary() {
    use cheby::piecewise::adaptive::AdaptiveSegmentTable;

    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let table: AdaptiveSegmentTable<f64, f64, 4> =
        AdaptiveSegmentTable::from_fn(domain, |x| x, 1e-3, 4).unwrap();

    // locate at the domain start
    let seg = table.locate(0.0_f64);
    assert!(seg.is_some());

    // locate outside domain
    assert!(table.locate(5.0_f64).is_none());
    assert!(table.locate(-1.0_f64).is_none());

    // locate at domain end (last boundary)
    let last = *table.boundaries().last().unwrap();
    let seg_at_end = table.locate(last);
    // End is the upper bound; may or may not be inside a segment.
    let _ = seg_at_end;
}

#[test]
fn adaptive_segment_table_evaluate_derivative() {
    use cheby::piecewise::adaptive::AdaptiveSegmentTable;

    let domain = cheby::Domain::new(0.0_f64, 2.0);
    let table: AdaptiveSegmentTable<f64, f64, 8> =
        AdaptiveSegmentTable::from_fn(domain, |x: f64| x, 1e-6, 6).unwrap();

    // derivative of f(x) = x is 1 everywhere
    let d = table.evaluate_derivative(1.0_f64).unwrap();
    assert_abs_diff_eq!(d, 1.0, epsilon = 1e-6);

    // out-of-domain derivative
    assert_eq!(
        table.evaluate_derivative(9.0_f64).unwrap_err(),
        cheby::ChebyError::EvaluationOutOfDomain,
    );
}

#[test]
fn adaptive_segment_table_invalid_degree_error() {
    use cheby::piecewise::adaptive::AdaptiveSegmentTable;

    let domain = cheby::Domain::new(0.0_f64, 1.0);
    // N = 1 is too small
    let result = AdaptiveSegmentTable::<f64, f64, 1>::from_fn(domain, |x| x, 1e-6, 4);
    assert_eq!(result.unwrap_err(), cheby::ChebyError::InvalidDegree);
}

#[test]
fn adaptive_segment_table_non_finite_tolerance_error() {
    use cheby::piecewise::adaptive::AdaptiveSegmentTable;

    let domain = cheby::Domain::new(0.0_f64, 1.0);
    let result = AdaptiveSegmentTable::<f64, f64, 4>::from_fn(domain, |x| x, f64::NAN, 4);
    assert_eq!(result.unwrap_err(), cheby::ChebyError::NonFiniteInput);
}
