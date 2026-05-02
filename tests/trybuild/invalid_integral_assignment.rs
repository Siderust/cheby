fn main() {
    type Velocity = qtty::velocity::Velocity<qtty::unit::Kilometer, qtty::unit::Second>;
    let domain = cheby::Domain::new(qtty::Second::new(0.0), qtty::Second::new(1.0));
    let series = cheby::ChebySeriesOn::new(
        domain,
        cheby::ChebySeries::new([Velocity::new(1.0), Velocity::new(0.0)]),
    );
    let position: Velocity = series
        .evaluate_integral_from_start(qtty::Second::new(0.5))
        .unwrap();
    let _ = position;
}
