use qtty::{Kilometer, Second};

type Velocity = qtty::velocity::Velocity<qtty::unit::Kilometer, qtty::unit::Second>;

fn main() {
    let domain = cheby::Domain::new(Second::new(0.0), Second::new(10.0));
    let series = cheby::ChebySeriesOn::new(
        domain,
        cheby::ChebySeries::new([Velocity::new(2.0), Velocity::new(0.0)]),
    );
    let position: Kilometer = series
        .evaluate_integral_from_start(Second::new(3.0))
        .unwrap();
    println!("{}", position.value());
}
