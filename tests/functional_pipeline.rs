use approx::assert_abs_diff_eq;
use cheby::{
    evaluate, evaluate_both, evaluate_derivative, fit_coeffs, fit_from_fn, nodes, nodes_mapped,
    ChebySegment, ChebySegmentTable,
};

#[test]
fn functional_roundtrip_on_unit_interval() {
    const N: usize = 17;
    let xi: [f64; N] = nodes();
    let values: [f64; N] = std::array::from_fn(|k| (2.5 * xi[k]).sin() + 0.2 * xi[k] * xi[k]);
    let coeffs = fit_coeffs(&values);

    for &tau in &[-0.95, -0.5, 0.0, 0.37, 0.9] {
        let approx = evaluate(&coeffs, tau);
        let exact = (2.5 * tau).sin() + 0.2 * tau * tau;
        assert_abs_diff_eq!(approx, exact, epsilon = 1e-12);
    }
}

#[test]
fn functional_roundtrip_on_physical_interval() {
    const N: usize = 21;
    let start = -2.0;
    let end = 4.0;
    let coeffs: [f64; N] = fit_from_fn(|t| t.exp() / 10.0, start, end);

    // t -> tau conversion for a mapped segment.
    let mid = 0.5 * (start + end);
    let half = 0.5 * (end - start);

    for &t in &[-1.8, -0.3, 1.0, 2.4, 3.7] {
        let tau = (t - mid) / half;
        let approx = evaluate(&coeffs, tau);
        let exact = t.exp() / 10.0;
        assert_abs_diff_eq!(approx, exact, epsilon = 1e-8);
    }
}

#[test]
fn functional_value_and_derivative_consistency() {
    let coeffs: [f64; 15] = fit_from_fn(f64::sin, -1.0, 1.0);

    for &tau in &[-0.8, -0.2, 0.0, 0.6] {
        let (v, d) = evaluate_both(&coeffs, tau);
        assert_abs_diff_eq!(v, evaluate(&coeffs, tau), epsilon = 1e-13);
        assert_abs_diff_eq!(d, evaluate_derivative(&coeffs, tau), epsilon = 1e-13);
    }
}

#[test]
fn nodes_mapping_matches_affine_formula() {
    const N: usize = 11;
    let start = 10.0;
    let end = 22.0;
    let mapped = nodes_mapped::<N>(start, end);
    let unit = nodes::<N>();
    let mid = 0.5 * (start + end);
    let half = 0.5 * (end - start);

    for k in 0..N {
        let expected = mid + half * unit[k];
        assert_abs_diff_eq!(mapped[k], expected, epsilon = 1e-15);
    }
}

#[test]
fn segment_table_end_to_end() {
    const N: usize = 13;
    let start = 0.0;
    let end = 6.0;
    let segment_len = 1.5;
    let table: ChebySegmentTable<f64, N> =
        ChebySegmentTable::from_fn(|t| 0.5 * t.cos() + 2.0, start, end, segment_len);

    assert_eq!(table.len(), 4);
    assert!(!table.is_empty());
    assert_abs_diff_eq!(table.start(), start, epsilon = 0.0);
    assert_abs_diff_eq!(table.end(), end, epsilon = 0.0);
    assert_abs_diff_eq!(table.segment_len(), segment_len, epsilon = 0.0);
    assert_eq!(table.segments().len(), table.len());

    for &t in &[0.1, 1.0, 2.1, 3.2, 4.7, 5.9] {
        let approx = table.eval(t).expect("in range");
        let exact = 0.5 * t.cos() + 2.0;
        assert_abs_diff_eq!(approx, exact, epsilon = 1e-9);
    }

    // Table uses half-open indexing: [start, end), so exactly at end is out of range.
    assert!(table.eval(start - 1e-6).is_none());
    assert!(table.eval(end).is_none());
}

#[test]
fn segment_table_from_precomputed_segments_and_empty_case() {
    const N: usize = 9;
    let start = 2.0;
    let segment_len = 2.0;
    let half = segment_len / 2.0;

    let mk = |seg_start: f64| -> ChebySegment<f64, N> {
        let coeffs = fit_from_fn::<f64, N>(|t| t * t + 1.0, seg_start, seg_start + segment_len);
        ChebySegment::new(coeffs, seg_start + half, half)
    };

    let table = ChebySegmentTable::from_segments(vec![mk(2.0), mk(4.0)], start, segment_len);
    assert_eq!(table.len(), 2);
    assert!(!table.is_empty());

    let t = 5.25;
    let approx = table.eval(t).expect("inside second segment");
    let exact = t * t + 1.0;
    assert_abs_diff_eq!(approx, exact, epsilon = 1e-9);

    let empty: ChebySegmentTable<f64, N> = ChebySegmentTable::from_segments(vec![], 10.0, 1.0);
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert!(empty.get_segment(10.0).is_none());
    assert!(empty.eval(10.0).is_none());
    assert!(empty.eval_derivative(10.0).is_none());
    assert!(empty.eval_both(10.0).is_none());
}
