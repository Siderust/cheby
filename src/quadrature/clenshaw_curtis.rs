//! Clenshaw-Curtis quadrature.

use crate::core::{ChebyScalar, ChebyTime, Domain, IntegrateWith, NodeKind};

/// Return simple Clenshaw-Curtis-like weights on `[-1, 1]`.
///
/// The implementation uses trapezoidal weights over Chebyshev-Lobatto
/// abscissae. It is stable for examples and small rules, and keeps the API
/// lightweight.
pub fn clenshaw_curtis_weights<const N: usize>() -> [f64; N] {
    if N == 0 {
        return [0.0; N];
    }
    if N == 1 {
        return [2.0; N];
    }
    let xs = crate::core::nodes::nodes::<N>(NodeKind::GaussLobatto);
    let mut weights = [0.0; N];
    for k in 0..N {
        let left = if k == 0 {
            xs[0]
        } else {
            0.5 * (xs[k - 1] + xs[k])
        };
        let right = if k + 1 == N {
            xs[N - 1]
        } else {
            0.5 * (xs[k] + xs[k + 1])
        };
        weights[k] = (left - right).abs();
    }
    weights
}

/// Integrate a function over `domain`.
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
