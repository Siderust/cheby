use approx::assert_abs_diff_eq;
use cheby::approx::interpolation::BarycentricInterpolator;
use cheby::core::{basis, nodes, ChebySeries, Domain, NodeKind};
use cheby::piecewise::ChebySegmentTable;
use proptest::prelude::*;

proptest! {
    #[test]
    fn domain_normalization_properties(start in -1e3_f64..1e3, width in 1e-3_f64..1e3) {
        let domain = Domain::try_new(start, start + width).unwrap();
        // Maximum floating-point error in normalize is O(eps * |start| / half).
        // With these ranges (|start| ≤ 1e3, half ≥ 5e-4), that is < 1e-9.
        assert_abs_diff_eq!(domain.normalize(domain.start()), -1.0, epsilon = 1e-9);
        assert_abs_diff_eq!(domain.normalize(domain.midpoint()), 0.0, epsilon = 1e-9);
        assert_abs_diff_eq!(domain.normalize(domain.end()), 1.0, epsilon = 1e-9);
    }

    #[test]
    fn chebyshev_recurrence_holds(x in -1.0_f64..1.0, n in 1_usize..20) {
        let lhs = basis::t(n + 1, x);
        let rhs = 2.0 * x * basis::t(n, x) - basis::t(n - 1, x);
        assert!((lhs - rhs).abs() < 1e-12);
    }
}

#[test]
fn interpolation_reproduces_samples() {
    const N: usize = 9;
    let xs = nodes::<N>(NodeKind::Roots);
    let ys = std::array::from_fn(|k| xs[k].sin());
    let interp = BarycentricInterpolator::new(xs, ys).unwrap();
    for k in 0..N {
        assert_abs_diff_eq!(interp.evaluate(xs[k]).unwrap(), ys[k], epsilon = 1e-12);
    }
}

#[test]
fn interpolation_single_node_is_constant_even_at_zero() {
    // Regression: a single node at zero used to make `scale == 0`, producing
    // a NaN division during evaluation instead of the constant interpolant.
    let interp = BarycentricInterpolator::new([0.0_f64], [42.0_f64]).unwrap();
    for &x in &[-2.5_f64, -0.1, 0.0, 0.7, 17.0] {
        assert_abs_diff_eq!(interp.evaluate(x).unwrap(), 42.0, epsilon = 0.0);
    }
}

#[test]
fn interpolation_chebyshev_roots_constructor_matches_generic() {
    use cheby::core::{nodes_mapped, Domain, NodeKind};
    const N: usize = 12;
    let domain = Domain::new(-1.0_f64, 3.0);
    let xs = nodes_mapped::<f64, N>(domain, NodeKind::Roots);
    let ys: [f64; N] = std::array::from_fn(|k| (xs[k] + 0.3).cos());
    let generic = BarycentricInterpolator::new(xs, ys).unwrap();
    let fast = BarycentricInterpolator::on_chebyshev_roots(domain, ys).unwrap();
    for tau in [-0.95_f64, -0.4, 0.0, 0.2, 0.85] {
        let x = domain.denormalize(tau);
        assert_abs_diff_eq!(
            fast.evaluate(x).unwrap(),
            generic.evaluate(x).unwrap(),
            epsilon = 1e-12
        );
    }
}

#[test]
fn interpolation_lobatto_constructor_reproduces_samples() {
    use cheby::core::{nodes_mapped, Domain, NodeKind};
    const N: usize = 9;
    let domain = Domain::new(0.0_f64, 2.0);
    let xs = nodes_mapped::<f64, N>(domain, NodeKind::Lobatto);
    let ys: [f64; N] = std::array::from_fn(|k| xs[k].sin());
    let interp = BarycentricInterpolator::on_lobatto_nodes(domain, ys).unwrap();
    for k in 0..N {
        assert_abs_diff_eq!(interp.evaluate(xs[k]).unwrap(), ys[k], epsilon = 1e-12);
    }
}

#[test]
fn derivative_of_integral_recovers_series_approximately() {
    let series = ChebySeries::new([1.0, -0.5, 0.25, 0.0, 0.0]);
    let recovered = series.integral(0.0).derivative();
    for tau in [-0.75, -0.25, 0.25, 0.75] {
        assert_abs_diff_eq!(
            recovered.evaluate(tau),
            series.evaluate(tau),
            epsilon = 1e-12
        );
    }
}

#[test]
fn piecewise_lookup_chooses_correct_segment() {
    let table: ChebySegmentTable<f64, f64, 5> = ChebySegmentTable::from_fn(|x| x, 0.0, 4.0, 1.0);
    assert_abs_diff_eq!(table.get_segment(0.25).unwrap().domain().start(), 0.0);
    assert_abs_diff_eq!(table.get_segment(1.25).unwrap().domain().start(), 1.0);
    assert_abs_diff_eq!(table.get_segment(3.25).unwrap().domain().start(), 3.0);
}

#[test]
#[cfg(feature = "adaptive")]
fn adaptive_fitting_converges_for_smooth_function() {
    let result =
        cheby::approx::adaptive::AdaptiveFit::new(Domain::new(-1.0, 1.0), |x: f64| x.sin())
            .max_degree(64)
            .absolute_tolerance(1e-9)
            .relative_tolerance(1e-9)
            .build::<f64>()
            .unwrap();
    assert!(result.report.converged);
}

#[test]
#[cfg(feature = "adaptive")]
fn adaptive_fitting_honors_non_power_of_two_limit() {
    let result =
        cheby::approx::adaptive::AdaptiveFit::new(Domain::new(-1.0, 1.0), |x: f64| (3.0 * x).sin())
            .max_degree(20)
            .absolute_tolerance(0.0)
            .relative_tolerance(0.0)
            .build::<f64>()
            .unwrap();
    assert_eq!(result.report.degree, 20);
}

#[test]
#[cfg(feature = "adaptive")]
fn adaptive_fitting_allows_small_degree_limit() {
    let result =
        cheby::approx::adaptive::AdaptiveFit::new(Domain::new(-1.0, 1.0), |x: f64| x.sin())
            .max_degree(7)
            .absolute_tolerance(0.0)
            .relative_tolerance(0.0)
            .build::<f64>()
            .unwrap();
    assert_eq!(result.report.degree, 7);
}
