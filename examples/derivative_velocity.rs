use qtty::{Kilometer, Second};

fn main() {
    let domain = cheby::Domain::new(Second::new(0.0), Second::new(10.0));
    let series = cheby::approx::fit::fit_from_fn_on::<Kilometer, Second, 16>(domain, |t| {
        Kilometer::new(t.value() * t.value())
    });
    let velocity = series.evaluate_derivative(Second::new(4.0)).unwrap();
    println!("{}", velocity.value());
}
