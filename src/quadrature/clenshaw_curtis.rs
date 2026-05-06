//! Clenshaw-Curtis quadrature on Chebyshev-Lobatto nodes.
//!
//! Computes weights `w_k` such that `sum_k w_k f(x_k)` integrates polynomials
//! of degree up to `N - 1` exactly on `[-1, 1]`, where `N` is the node count
//! and `x_k = cos(k π / (N - 1))` are the Chebyshev extrema (Lobatto nodes).
//!
//! The implementation follows Trefethen, *Spectral Methods in MATLAB*,
//! Program 30 (`clencurt`). Cost is `O(N^2)`.

use crate::core::{ChebyScalar, ChebyTime, Domain, IntegrateWith, NodeKind};

/// Clenshaw-Curtis weights for `N` Lobatto nodes on `[-1, 1]`.
pub fn clenshaw_curtis_weights<const N: usize>() -> [f64; N] {
    let mut weights = [0.0; N];
    if N == 0 {
        return weights;
    }
    if N == 1 {
        weights[0] = 2.0;
        return weights;
    }
    if N == 2 {
        weights[0] = 1.0;
        weights[1] = 1.0;
        return weights;
    }

    let m = N - 1; // Trefethen's "N"
    let theta: [f64; N] =
        core::array::from_fn(|k| core::f64::consts::PI * k as f64 / m as f64);

    // Endpoint weights.
    if m.is_multiple_of(2) {
        weights[0] = 1.0 / ((m as f64) * (m as f64) - 1.0);
        weights[N - 1] = weights[0];
    } else {
        weights[0] = 1.0 / (m as f64 * m as f64);
        weights[N - 1] = weights[0];
    }

    // Interior weights.
    for (k, w) in weights.iter_mut().enumerate().take(N - 1).skip(1) {
        let mut v = 1.0_f64;
        let half = m / 2;
        let upper = if m.is_multiple_of(2) { half - 1 } else { half };
        for j in 1..=upper {
            v -= 2.0 * (2.0 * j as f64 * theta[k]).cos() / (4.0 * (j * j) as f64 - 1.0);
        }
        if m.is_multiple_of(2) {
            v -= (m as f64 * theta[k]).cos() / ((m as f64) * (m as f64) - 1.0);
        }
        *w = 2.0 * v / m as f64;
    }

    weights
}

/// Integrate a function over `domain` using `N`-point Clenshaw-Curtis.
pub fn integrate<T, X, const N: usize>(
    domain: Domain<X>,
    f: impl Fn(X) -> T,
) -> <T as IntegrateWith<X>>::Integral
where
    T: ChebyScalar + IntegrateWith<X>,
    X: ChebyTime,
{
    let nodes = crate::core::nodes::nodes::<N>(NodeKind::GaussLobatto);
    let weights = clenshaw_curtis_weights::<N>();
    let mut sum = T::zero();
    for k in 0..N {
        sum = sum + f(domain.denormalize(nodes[k])) * weights[k];
    }
    sum.scale_integral(domain.half_width())
}
