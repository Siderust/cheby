use cheby::ChebySeries;

fn main() {
    let series = ChebySeries::new([1.0, 2.0, 0.5]);
    println!("{}", series.evaluate(0.25));
}
