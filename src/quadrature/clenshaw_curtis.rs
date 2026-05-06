//! Clenshaw-Curtis quadrature on Chebyshev-Lobatto nodes.
//!
//! Computes weights `w_k` such that `sum_k w_k f(x_k)` integrates polynomials
//! of degree up to `N - 1` exactly on `[-1, 1]`, where `N` is the node count
//! and `x_k = cos(k π / (N - 1))` are the Chebyshev extrema (Lobatto nodes).
//!
//! The implementation follows Trefethen, *Spectral Methods in MATLAB*,
//! Program 30 (`clencurt`). Cost is `O(N^2)`.
//!
//! For repeated quadrature on the same `N`, build a [`ClenshawCurtisRule`]
//! once and reuse it; the per-call [`integrate`] convenience rebuilds nodes
//! and weights every time.

use crate::core::{ChebyScalar, ChebyTime, Domain, IntegrateWith, NodeKind};

/// A precomputed `N`-point Clenshaw-Curtis rule on `[-1, 1]`.
///
/// Holds the Chebyshev-Lobatto nodes and the matching weights so that
/// repeated calls do not pay the `O(N²)` weight construction every time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClenshawCurtisRule<const N: usize> {
    nodes: [f64; N],
    weights: [f64; N],
}

impl<const N: usize> Default for ClenshawCurtisRule<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> ClenshawCurtisRule<N> {
    /// Build the cached rule.
    #[inline]
    pub fn new() -> Self {
        Self {
            nodes: crate::core::nodes::nodes::<N>(NodeKind::GaussLobatto),
            weights: clenshaw_curtis_weights::<N>(),
        }
    }

    /// Borrow the cached Lobatto nodes on `[-1, 1]`.
    #[inline]
    pub fn nodes(&self) -> &[f64; N] {
        &self.nodes
    }

    /// Borrow the cached quadrature weights.
    #[inline]
    pub fn weights(&self) -> &[f64; N] {
        &self.weights
    }

    /// Integrate `f` over `domain` using this cached rule.
    #[inline]
    pub fn integrate<T, X>(
        &self,
        domain: Domain<X>,
        f: impl Fn(X) -> T,
    ) -> <T as IntegrateWith<X>>::Integral
    where
        T: ChebyScalar + IntegrateWith<X>,
        X: ChebyTime,
    {
        let mut sum = T::zero();
        for k in 0..N {
            sum = sum + f(domain.denormalize(self.nodes[k])) * self.weights[k];
        }
        sum.scale_integral(domain.half_width())
    }
}

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
    let theta: [f64; N] = core::array::from_fn(|k| core::f64::consts::PI * k as f64 / m as f64);

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
///
/// Convenience wrapper that builds a fresh [`ClenshawCurtisRule`] each call.
/// For repeated quadrature on the same `N`, construct the rule once and call
/// [`ClenshawCurtisRule::integrate`].
pub fn integrate<T, X, const N: usize>(
    domain: Domain<X>,
    f: impl Fn(X) -> T,
) -> <T as IntegrateWith<X>>::Integral
where
    T: ChebyScalar + IntegrateWith<X>,
    X: ChebyTime,
{
    ClenshawCurtisRule::<N>::new().integrate(domain, f)
}
