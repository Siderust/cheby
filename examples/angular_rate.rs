use qtty::{Radian, Second};

type AngularRate = qtty::angular_rate::AngularRate<qtty::unit::Radian, qtty::unit::Second>;

fn main() {
    let domain = cheby::Domain::new(Second::new(0.0), Second::new(2.0));
    let series = cheby::ChebySeriesOn::new(
        domain,
        cheby::ChebySeries::new([AngularRate::new(1.0), AngularRate::new(0.0)]),
    );
    let angle: Radian = series
        .evaluate_integral_from_start(Second::new(1.0))
        .unwrap();
    println!("{}", angle.value());
}
