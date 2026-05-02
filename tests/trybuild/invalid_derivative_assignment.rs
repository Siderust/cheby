fn main() {
    let domain = cheby::Domain::new(qtty::Second::new(0.0), qtty::Second::new(1.0));
    let series = cheby::ChebySeriesOn::new(
        domain,
        cheby::ChebySeries::new([qtty::Kilometer::new(0.0), qtty::Kilometer::new(1.0)]),
    );
    let velocity: qtty::Kilometer = series
        .evaluate_derivative(qtty::Second::new(0.5))
        .unwrap();
    let _ = velocity;
}
