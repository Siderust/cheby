//! Hot-path benchmarks for `cheby`.
//!
//! Covers, in order: scalar evaluation, derivative evaluation (raw + typed),
//! fixed-size fitting, adaptive fitting, minimax (Remez), piecewise lookup
//! and evaluation, Clenshaw-Curtis quadrature, Gauss-Chebyshev quadrature,
//! and spectral differentiation matrix construction.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use qtty::{Kilometer, Second};

fn evaluate_f64(c: &mut Criterion) {
    let coeffs = [1.0, 0.5, -0.25, 0.125, -0.0625, 0.03125, 0.0, 0.0];
    c.bench_function("evaluate_f64", |b| {
        b.iter(|| cheby::evaluate(black_box(&coeffs), black_box(0.25)))
    });
}

fn evaluate_quantity(c: &mut Criterion) {
    let coeffs = [
        Kilometer::new(1.0),
        Kilometer::new(0.5),
        Kilometer::new(-0.25),
        Kilometer::new(0.125),
    ];
    c.bench_function("evaluate_quantity", |b| {
        b.iter(|| cheby::evaluate(black_box(&coeffs), black_box(0.25)))
    });
}

fn evaluate_derivative_f64(c: &mut Criterion) {
    let coeffs = [1.0, 0.5, -0.25, 0.125, -0.0625, 0.03125, 0.0, 0.0];
    c.bench_function("evaluate_derivative_f64", |b| {
        b.iter(|| cheby::evaluate_derivative(black_box(&coeffs), black_box(0.25)))
    });
}

fn evaluate_derivative_quantity(c: &mut Criterion) {
    let domain = cheby::Domain::new(Second::new(0.0), Second::new(1.0));
    let series = cheby::ChebySeriesOn::new(
        domain,
        cheby::ChebySeries::new([Kilometer::new(1.0), Kilometer::new(0.5)]),
    );
    c.bench_function("evaluate_derivative_quantity", |b| {
        b.iter(|| {
            series
                .evaluate_derivative(black_box(Second::new(0.25)))
                .unwrap()
        })
    });
}

fn fit_fixed(c: &mut Criterion) {
    c.bench_function("fit_fixed_n16_sin", |b| {
        b.iter(|| cheby::approx::fit::fit_from_fn::<f64, 16>(|x| black_box(x).sin()))
    });
}

fn fit_quantity(c: &mut Criterion) {
    c.bench_function("fit_quantity_n16", |b| {
        b.iter(|| {
            cheby::approx::fit::fit_from_fn::<Kilometer, 16>(|x| Kilometer::new(black_box(x).sin()))
        })
    });
}

fn adaptive_fit(c: &mut Criterion) {
    c.bench_function("adaptive_fit_sin", |b| {
        b.iter(|| {
            cheby::approx::adaptive::AdaptiveFit::new(cheby::Domain::new(-1.0, 1.0), |x| {
                f64::sin(x)
            })
            .max_degree(32)
            .build::<f64>()
            .unwrap()
        })
    });
}

fn minimax_remez(c: &mut Criterion) {
    c.bench_function("minimax_remez_exp_deg7", |b| {
        b.iter(|| {
            cheby::approx::minimax::remez(
                cheby::Domain::new(-1.0, 1.0),
                7,
                |x: f64| black_box(x).exp(),
                cheby::approx::minimax::RemezOptions::default(),
            )
            .unwrap()
        })
    });
}

fn piecewise_lookup(c: &mut Criterion) {
    let table: cheby::ChebySegmentTable<f64, f64, 8> =
        cheby::ChebySegmentTable::from_fn(|x| x.sin(), 0.0, 8.0, 1.0);
    c.bench_function("piecewise_lookup", |b| {
        b.iter(|| table.get_segment(black_box(4.25)))
    });
}

fn piecewise_evaluate(c: &mut Criterion) {
    let table: cheby::ChebySegmentTable<f64, f64, 8> =
        cheby::ChebySegmentTable::from_fn(|x| x.sin(), 0.0, 8.0, 1.0);
    c.bench_function("piecewise_evaluate", |b| {
        b.iter(|| table.evaluate(black_box(4.25)).unwrap())
    });
}

fn adaptive_piecewise_lookup(c: &mut Criterion) {
    let table: cheby::piecewise::AdaptiveSegmentTable<f64, f64, 16> =
        cheby::piecewise::AdaptiveSegmentTable::from_fn(
            cheby::Domain::new(0.0, 8.0),
            |x: f64| x.sin(),
            1e-10,
            8,
        )
        .unwrap();
    c.bench_function("adaptive_piecewise_lookup", |b| {
        b.iter(|| table.evaluate(black_box(4.25)).unwrap())
    });
}

fn clenshaw_curtis_integrate(c: &mut Criterion) {
    c.bench_function("clenshaw_curtis_integrate_n32", |b| {
        b.iter(|| {
            cheby::quadrature::clenshaw_curtis::integrate::<f64, f64, 32>(
                cheby::Domain::new(0.0, 1.0),
                |x| black_box(x).exp(),
            )
        })
    });
}

fn gauss_chebyshev(c: &mut Criterion) {
    c.bench_function("gauss_chebyshev_weighted_n32", |b| {
        b.iter(|| {
            cheby::quadrature::gauss_chebyshev::integrate_weighted::<f64, 32>(|x| {
                black_box(x).cos()
            })
        })
    });
}

fn spectral_diff_matrix(c: &mut Criterion) {
    c.bench_function("spectral_diff_matrix_n16", |b| {
        b.iter(|| cheby::spectral::chebyshev_differentiation_matrix(black_box(16)))
    });
}

criterion_group!(
    benches,
    evaluate_f64,
    evaluate_quantity,
    evaluate_derivative_f64,
    evaluate_derivative_quantity,
    fit_fixed,
    fit_quantity,
    adaptive_fit,
    minimax_remez,
    piecewise_lookup,
    piecewise_evaluate,
    adaptive_piecewise_lookup,
    clenshaw_curtis_integrate,
    gauss_chebyshev,
    spectral_diff_matrix
);
criterion_main!(benches);
