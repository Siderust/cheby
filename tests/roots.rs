//! Comprehensive tests for dynamic Chebyshev root finding.

use cheby::approx::fit::fit_dyn_from_fn;
use cheby::{ChebyError, ChebySeriesDyn, Domain, RootOptions};
use proptest::prelude::*;

fn opts() -> RootOptions {
    RootOptions {
        zero_tol: 1e-9,
        unit_tol: 1e-11,
        dedupe_eps: 1e-10,
    }
}

fn poly(coeffs: &[f64]) -> ChebySeriesDyn<f64> {
    ChebySeriesDyn::new(coeffs.to_vec()).unwrap()
}

fn assert_residuals(poly: &ChebySeriesDyn<f64>, roots: &[f64], zero_tol: f64) {
    for &root in roots {
        assert!((-1.0..=1.0).contains(&root), "root out of range: {root}");
        let residual = poly.evaluate(root).abs();
        assert!(residual <= zero_tol, "residual {residual} at root {root}");
    }
}

#[test]
fn linear_root() {
    // p(tau) = tau - 0.25  =>  c0 = -0.25, c1 = 1
    let p = poly(&[-0.25, 1.0]);
    let roots = p.roots_with(opts());
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 0.25).abs() < 1e-10);
    assert_residuals(&p, &roots, opts().zero_tol);
}

#[test]
fn root_at_minus_one() {
    // p(tau) = tau + 1
    let p = poly(&[1.0, 1.0]);
    let roots = p.roots_with(opts());
    assert_eq!(roots.len(), 1);
    assert!((roots[0] + 1.0).abs() < 1e-10);
    assert_residuals(&p, &roots, opts().zero_tol);
}

#[test]
fn root_at_plus_one() {
    // p(tau) = 1 - tau
    let p = poly(&[1.0, -1.0]);
    let roots = p.roots_with(opts());
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 1.0).abs() < 1e-10);
    assert_residuals(&p, &roots, opts().zero_tol);
}

#[test]
fn no_real_roots_on_interval() {
    // p(tau) = tau^2 + 1  (via direct monomial-ish fit is poor; use evaluation form)
    let p = fit_dyn_from_fn(8, |x| x * x + 1.0).unwrap();
    let roots = p.roots_with(opts());
    assert!(roots.is_empty(), "{roots:?}");
}

#[test]
fn double_root_tangency() {
    let p = fit_dyn_from_fn(12, |x| (x - 0.3).powi(2)).unwrap();
    let roots = p.roots_with(RootOptions {
        zero_tol: 1e-6,
        ..opts()
    });
    assert!(!roots.is_empty());
    assert!(roots.iter().any(|r| (*r - 0.3).abs() < 1e-3));
    assert_residuals(&p, &roots, 1e-6);
}

#[test]
fn triple_root() {
    let p = fit_dyn_from_fn(16, |x| (x + 0.2).powi(3)).unwrap();
    let roots = p.roots_with(RootOptions {
        zero_tol: 1e-5,
        ..opts()
    });
    assert!(!roots.is_empty());
    assert!(roots.iter().any(|r| (*r + 0.2).abs() < 5e-3));
}

#[test]
fn multiple_separated_roots() {
    let p = fit_dyn_from_fn(16, |x| x * x * x - x).unwrap();
    let roots = p.roots_with(opts());
    assert_eq!(roots.len(), 3, "{roots:?}");
    assert_residuals(&p, &roots, opts().zero_tol);
}

#[test]
fn very_close_roots() {
    let p = fit_dyn_from_fn(24, |x| (x - 0.40) * (x - 0.41)).unwrap();
    let roots = p.roots_with(RootOptions {
        zero_tol: 1e-7,
        dedupe_eps: 1e-4,
        ..opts()
    });
    assert!(roots.len() >= 2, "{roots:?}");
    assert_residuals(&p, &roots, 1e-7);
}

#[test]
fn chebyshev_polynomial_roots() {
    for n in [4usize, 8, 12] {
        let mut coeffs = vec![0.0_f64; n + 1];
        coeffs[n] = 1.0;
        let p = ChebySeriesDyn::new(coeffs).unwrap();
        let options = opts();
        let roots = p.roots_with(options);
        assert_eq!(roots.len(), n, "T_{n} should have {n} roots, got {roots:?}");
        assert_residuals(&p, &roots, options.zero_tol);
    }
}

#[test]
fn non_finite_coefficients_yield_no_roots() {
    let p = poly(&[1.0, f64::NAN]);
    assert!(p.roots_with(opts()).is_empty());
    let p = poly(&[f64::INFINITY, 1.0]);
    assert!(p.roots_with(opts()).is_empty());
}

#[test]
fn invalid_root_options_are_sanitized() {
    let p = poly(&[-0.25, 1.0]);
    let bad = RootOptions {
        unit_tol: f64::NAN,
        zero_tol: -1.0,
        dedupe_eps: 0.0,
    };
    assert_eq!(bad.validated(), Err(ChebyError::NonFiniteInput));
    let roots = p.roots_with(bad);
    assert_eq!(roots.len(), 1);
}

#[test]
fn deduplication_merges_nearby_roots() {
    let p = fit_dyn_from_fn(10, |x| x - 0.5).unwrap();
    let roots = p.roots_with(RootOptions {
        dedupe_eps: 1.0,
        ..opts()
    });
    assert_eq!(roots.len(), 1);
}

#[test]
fn shifted_constant_finds_offset_roots() {
    let base = fit_dyn_from_fn(10, |x| x).unwrap();
    let shifted = base.shifted_constant(-0.5);
    let roots = shifted.roots_with(opts());
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 0.5).abs() < 1e-8);
}

#[test]
fn domain_mapped_roots() {
    let domain = Domain::new(0.0_f64, 10.0);
    let series = fit_dyn_from_fn(10, |tau| tau - 0.0).unwrap();
    let on = cheby::ChebySeriesDynOn::new(domain, series);
    let roots = on.roots_with(opts());
    assert_eq!(roots.len(), 1);
    assert!((roots[0] - 5.0).abs() < 1e-8);
}

#[test]
fn constant_zero_series_has_no_roots() {
    assert!(poly(&[0.0]).roots_with(opts()).is_empty());
    assert!(poly(&[1e-15]).roots_with(opts()).is_empty());
}

#[test]
fn constant_nonzero_series_has_no_roots() {
    assert!(poly(&[1.0]).roots_with(opts()).is_empty());
    assert!(poly(&[-2.5]).roots_with(opts()).is_empty());
}

#[test]
fn degenerate_linear_has_no_roots() {
    let p = poly(&[1.0, 1e-15]);
    assert!(p.roots_with(opts()).is_empty());
}

#[test]
fn fit_degree_zero_produces_one_coefficient() {
    let p = fit_dyn_from_fn(0, |_| 3.5).unwrap();
    assert_eq!(p.coeffs().len(), 1);
    assert!((p.evaluate(0.0) - 3.5).abs() < 1e-12);
    assert!(p.roots_with(opts()).is_empty());
}

#[test]
fn fit_degree_one_produces_two_coefficients() {
    let p = fit_dyn_from_fn(1, |x| 2.0 * x + 1.0).unwrap();
    assert_eq!(p.coeffs().len(), 2);
    let roots = p.roots_with(opts());
    assert_eq!(roots.len(), 1);
    assert!((roots[0] + 0.5).abs() < 1e-8);
    assert_residuals(&p, &roots, opts().zero_tol);
}

#[test]
fn chebyshev_t16_high_degree_roots() {
    let n = 16usize;
    let mut coeffs = vec![0.0_f64; n + 1];
    coeffs[n] = 1.0;
    let p = ChebySeriesDyn::new(coeffs).unwrap();
    let options = opts();
    let roots = p.roots_with(options);
    assert_eq!(roots.len(), n, "{roots:?}");
    assert_residuals(&p, &roots, options.zero_tol);
}

#[test]
fn high_degree_t24_may_miss_roots_but_residuals_hold() {
    let n = 24usize;
    let mut coeffs = vec![0.0_f64; n + 1];
    coeffs[n] = 1.0;
    let p = ChebySeriesDyn::new(coeffs).unwrap();
    let options = opts();
    let roots = p.roots_with(options);
    assert!(
        roots.len() >= n / 2,
        "expected most T_24 roots, got {roots:?}"
    );
    assert!(roots.len() <= n);
    assert_residuals(&p, &roots, options.zero_tol);
}

#[test]
fn strict_zero_tol_can_exclude_approximate_roots() {
    let p = fit_dyn_from_fn(12, |x| (x - 0.3).powi(2)).unwrap();
    let loose = p.roots_with(RootOptions {
        zero_tol: 1e-6,
        ..opts()
    });
    let strict = p.roots_with(RootOptions {
        zero_tol: 1e-14,
        ..opts()
    });
    assert!(!loose.is_empty());
    assert!(strict.is_empty() || strict.len() <= loose.len());
    assert_residuals(&p, &loose, 1e-6);
}

proptest! {
    #[test]
    fn fitted_quadratic_roots_satisfy_residual(a in -0.9f64..0.9, b in -0.9f64..0.9) {
        let p = fit_dyn_from_fn(12, move |x| (x - a) * (x - b)).unwrap();
        let options = opts();
        let roots = p.roots_with(options);
        prop_assert!(!roots.is_empty());
        for root in roots {
            prop_assert!((-1.0..=1.0).contains(&root));
            prop_assert!(p.evaluate(root).abs() <= options.zero_tol);
        }
    }
}

#[test]
fn fixed_and_dynamic_fit_agree_on_samples() {
    use cheby::approx::fit::{fit_coeffs, fit_coeffs_dyn, fit_from_fn};
    let values: [f64; 8] = core::array::from_fn(|k| {
        let n = 8.0;
        let x = (std::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * n)).cos();
        x.sin()
    });
    let xs = cheby::nodes::<8>();
    let fixed = fit_coeffs(&values);
    let dynamic = fit_coeffs_dyn(&values, &xs);
    for (a, b) in fixed.iter().zip(dynamic.iter()) {
        assert!((a - b).abs() < 1e-14, "fixed={a} dynamic={b}");
    }
    let series = fit_from_fn::<f64, 8>(|x| x.sin());
    for (a, b) in series.coeffs().iter().zip(dynamic.iter()) {
        assert!((a - b).abs() < 1e-14);
    }
}
