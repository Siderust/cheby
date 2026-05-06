//! Coverage for adaptive piecewise tables and binary-search lookup.
#![cfg(feature = "piecewise")]

use approx::assert_abs_diff_eq;
use cheby::core::{ChebyError, Domain};
use cheby::piecewise::AdaptiveSegmentTable;

#[test]
fn adaptive_table_splits_for_steep_region() {
    let domain = Domain::new(-1.0, 1.0);
    // tanh(20 x) is steep near zero and forces splitting.
    let table: AdaptiveSegmentTable<f64, f64, 16> =
        AdaptiveSegmentTable::from_fn(domain, |x: f64| (20.0 * x).tanh(), 1e-8, 8).unwrap();
    assert!(table.len() > 1, "expected >1 segments, got {}", table.len());
    // Boundaries are sorted and contiguous.
    let b = table.boundaries();
    for w in b.windows(2) {
        assert!(w[0] <= w[1]);
    }
    assert_eq!(b.first().copied().unwrap(), -1.0);
    assert_eq!(b.last().copied().unwrap(), 1.0);
    // Evaluation accuracy.
    for i in 0..50 {
        let x = -1.0 + 2.0 * i as f64 / 49.0;
        let approx = table.evaluate(x).unwrap();
        let exact = (20.0 * x).tanh();
        assert_abs_diff_eq!(approx, exact, epsilon = 1e-6);
    }
}

#[test]
fn adaptive_table_lookup_is_o_log_n_correct() {
    let domain = Domain::new(0.0, 16.0);
    let table: AdaptiveSegmentTable<f64, f64, 12> =
        AdaptiveSegmentTable::from_fn(domain, |x: f64| (5.0 * x).sin(), 1e-6, 10).unwrap();

    // Check every probe lands on a segment whose domain contains it.
    for i in 0..200 {
        let x = 16.0 * i as f64 / 199.0;
        let seg = table.locate(x).expect("must find a segment");
        assert!(seg.contains(x), "segment must contain probed x");
    }
}

#[test]
fn adaptive_table_out_of_domain_errors() {
    let domain = Domain::new(0.0, 1.0);
    let table: AdaptiveSegmentTable<f64, f64, 8> =
        AdaptiveSegmentTable::from_fn(domain, |x: f64| x * x, 1e-12, 4).unwrap();
    assert_eq!(
        table.evaluate(-0.1).unwrap_err(),
        ChebyError::EvaluationOutOfDomain
    );
    assert_eq!(
        table.evaluate(1.5).unwrap_err(),
        ChebyError::EvaluationOutOfDomain
    );
}

#[test]
fn adaptive_table_metadata_records_depth() {
    let domain = Domain::new(-1.0, 1.0);
    let table: AdaptiveSegmentTable<f64, f64, 8> =
        AdaptiveSegmentTable::from_fn(domain, |x: f64| (50.0 * x).tanh(), 1e-6, 8).unwrap();
    let meta = table.metadata();
    assert_eq!(meta.len(), table.len());
    assert!(meta.iter().any(|m| m.depth > 0));
}

#[test]
fn adaptive_table_rejects_invalid_inputs() {
    let domain = Domain::new(-1.0, 1.0);
    let r: Result<AdaptiveSegmentTable<f64, f64, 8>, _> =
        AdaptiveSegmentTable::from_fn(domain, |_| 0.0, f64::NAN, 4);
    assert_eq!(r.unwrap_err(), ChebyError::NonFiniteInput);
    let r: Result<AdaptiveSegmentTable<f64, f64, 8>, _> =
        AdaptiveSegmentTable::from_fn(domain, |_| 0.0, -1.0, 4);
    assert_eq!(r.unwrap_err(), ChebyError::NonFiniteInput);
}
