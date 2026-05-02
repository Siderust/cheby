//! Gauss-Chebyshev quadrature.

use crate::core::{ChebyScalar, NodeKind};

/// Integrate `f(x) / sqrt(1 - x^2)` on `[-1, 1]` with `N` roots.
pub fn integrate_weighted<T: ChebyScalar, const N: usize>(f: impl Fn(f64) -> T) -> T {
    let xs = crate::core::nodes::nodes::<N>(NodeKind::Gauss);
    let mut sum = T::zero();
    for x in xs {
        sum = sum + f(x);
    }
    sum * (core::f64::consts::PI / N as f64)
}
