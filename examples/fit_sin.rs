fn main() {
    let series = cheby::approx::fit::fit_from_fn::<f64, 16>(f64::sin);
    println!("{}", series.evaluate(0.5));
}
