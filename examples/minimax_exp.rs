fn main() {
    let result = cheby::approx::minimax::remez(
        cheby::Domain::new(-1.0, 1.0),
        7,
        f64::exp,
        cheby::approx::minimax::RemezOptions::default(),
    )
    .unwrap();
    println!("{} {}", result.iterations, result.max_error);
}
