fn main() {
    let table: cheby::ChebySegmentTable<f64, f64, 12> =
        cheby::ChebySegmentTable::from_fn(|t| t.sin(), 0.0, 4.0, 1.0);
    println!("{}", table.evaluate(1.25).unwrap());
}
