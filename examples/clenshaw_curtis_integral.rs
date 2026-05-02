fn main() {
    let integral = cheby::quadrature::clenshaw_curtis::integrate::<f64, f64, 33>(
        cheby::Domain::new(0.0, core::f64::consts::PI),
        f64::sin,
    );
    println!("{}", integral);
}
