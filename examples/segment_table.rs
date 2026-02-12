//! Piecewise Chebyshev approximation over a time range.
//!
//! Run with:
//! `cargo run --example segment_table`

use cheby::ChebySegmentTable;

fn main() {
    // Model function for demo purposes.
    let f = |t: f64| (0.8 * t).sin();

    // Build 4 uniform segments across [0, 8).
    let table: ChebySegmentTable<f64, 11> = ChebySegmentTable::from_fn(f, 0.0, 8.0, 2.0);

    println!(
        "segments={}, start={}, end={}, segment_len={}",
        table.len(),
        table.start(),
        table.end(),
        table.segment_len()
    );

    for &t in &[0.25, 1.50, 3.10, 6.75] {
        let (value, deriv) = table
            .eval_both(t)
            .expect("point must be inside the table range");
        println!(
            "t={t:.2}: value={value:.10}, derivative={deriv:.10}, exact={:.10}",
            f(t)
        );
    }
}
