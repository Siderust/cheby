fn main() {
    let altitude =
        cheby::approx::fit::fit_from_fn::<f64, 16>(|hour_angle| (0.6 * hour_angle).sin());
    println!("{}", altitude.evaluate(0.2));
}
