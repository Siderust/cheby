use approx::assert_abs_diff_eq;
use cheby::approx::interpolation::BarycentricInterpolator;
use cheby::core::{basis, nodes, ChebySeries, Domain, NodeKind};
use cheby::piecewise::ChebySegmentTable;
use proptest::prelude::*;

proptest! {
    #[test]
    fn domain_normalization_properties(start in -1e6_f64..1e6, width in 1e-9_f64..1e6) {
        let domain = Domain::try_new(start, start + width).unwrap();
        assert_abs_diff_eq!(domain.normalize(domain.start()), -1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(domain.normalize(domain.midpoint()), 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(domain.normalize(domain.end()), 1.0, epsilon = 1e-12);
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
