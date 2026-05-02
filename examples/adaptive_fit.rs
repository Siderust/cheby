fn main() {
    let result = cheby::approx::adaptive::AdaptiveFit::new(cheby::Domain::new(-1.0, 1.0), f64::sin)
        .max_degree(64)
        .absolute_tolerance(1e-10)
        .build::<f64>()
        .unwrap();
    println!("{} {}", result.report.degree, result.report.estimated_error);
}
