fn main() {
    let table: cheby::ChebySegmentTable<f64, f64, 10> =
        cheby::ChebySegmentTable::from_fn(|t| (0.1 * t).cos(), 0.0, 32.0, 4.0);
    println!("{}", table.evaluate(12.0).unwrap());
}
