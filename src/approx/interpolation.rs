//! Barycentric interpolation.

use crate::core::{nodes, ChebyError, ChebyScalar, ChebyTime, Domain, NodeKind};

/// Barycentric interpolator over supplied nodes and values.
///
/// For `N == 1`, the interpolant is the constant `values[0]` and `scale` is
/// unused; for `N > 1` `scale = nodes[0] - nodes[N-1]` and must be non-zero.
#[derive(Debug, Clone, PartialEq)]
pub struct BarycentricInterpolator<T, X, const N: usize> {
    nodes: [X; N],
    values: [T; N],
    weights: [f64; N],
    scale: X,
}

impl<T, X, const N: usize> BarycentricInterpolator<T, X, N>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    /// Create an interpolator. Nodes must be distinct.
    ///
    /// Generic O(N²) weight construction. For Chebyshev roots or Lobatto
    /// nodes prefer [`Self::on_chebyshev_roots`] or [`Self::on_lobatto_nodes`],
    /// which use the closed-form weights and are both faster and numerically
    /// more stable.
    pub fn new(nodes: [X; N], values: [T; N]) -> Result<Self, ChebyError> {
        if N == 0 {
            return Err(ChebyError::InvalidDegree);
        }
        if N == 1 {
            return Ok(Self {
                nodes,
                values,
                weights: [1.0; N],
                scale: X::zero(),
            });
        }
        let scale = nodes[0] - nodes[N - 1];
        if scale == X::zero() {
            return Err(ChebyError::InvalidDomain);
        }
        let mut weights = [0.0; N];
        for j in 0..N {
            let mut product = 1.0;
            for k in 0..N {
                if j != k {
                    let d = (nodes[j] - nodes[k]) / scale;
                    if d == 0.0 {
                        return Err(ChebyError::InvalidDomain);
                    }
                    product *= d;
                }
            }
            weights[j] = 1.0 / product;
        }
        Ok(Self {
            nodes,
            values,
            weights,
            scale,
        })
    }

    /// Build a barycentric interpolant on Chebyshev roots of `T_N` mapped
    /// into `domain`, using the closed-form weights
    /// `w_k = (-1)^k sin((2k+1)π / (2N))`.
    ///
    /// Setup is `O(N)` and avoids the `O(N²)` weight product chain. The
    /// scaling cancels in the barycentric formula, so absolute magnitudes do
    /// not matter.
    pub fn on_chebyshev_roots(domain: Domain<X>, values: [T; N]) -> Result<Self, ChebyError> {
        if N == 0 {
            return Err(ChebyError::InvalidDegree);
        }
        let xs = nodes::<N>(NodeKind::Roots);
        let mapped: [X; N] = core::array::from_fn(|k| domain.denormalize(xs[k]));
        if N == 1 {
            return Ok(Self {
                nodes: mapped,
                values,
                weights: [1.0; N],
                scale: X::zero(),
            });
        }
        let mut weights = [0.0; N];
        let mut sign = 1.0_f64;
        for (k, w) in weights.iter_mut().enumerate() {
            let theta = core::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * N as f64);
            *w = sign * theta.sin();
            sign = -sign;
        }
        let scale = mapped[0] - mapped[N - 1];
        Ok(Self {
            nodes: mapped,
            values,
            weights,
            scale,
        })
    }

    /// Build a barycentric interpolant on Chebyshev-Lobatto nodes of order
    /// `N - 1` mapped into `domain`, using the closed-form weights
    /// `w_0 = 1/2, w_{N-1} = (-1)^{N-1}/2`, interior `w_k = (-1)^k`.
    ///
    /// Setup is `O(N)`.
    pub fn on_lobatto_nodes(domain: Domain<X>, values: [T; N]) -> Result<Self, ChebyError> {
        if N == 0 {
            return Err(ChebyError::InvalidDegree);
        }
        let xs = nodes::<N>(NodeKind::Lobatto);
        let mapped: [X; N] = core::array::from_fn(|k| domain.denormalize(xs[k]));
        if N == 1 {
            return Ok(Self {
                nodes: mapped,
                values,
                weights: [1.0; N],
                scale: X::zero(),
            });
        }
        let mut weights = [0.0; N];
        let mut sign = 1.0_f64;
        for w in weights.iter_mut() {
            *w = sign;
            sign = -sign;
        }
        weights[0] *= 0.5;
        weights[N - 1] *= 0.5;
        let scale = mapped[0] - mapped[N - 1];
        Ok(Self {
            nodes: mapped,
            values,
            weights,
            scale,
        })
    }

    /// Evaluate the interpolant at `x`.
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
        if N == 1 {
            return Ok(self.values[0]);
        }
        let mut numerator = T::zero();
        let mut denominator = 0.0;
        for k in 0..N {
            let d = (x - self.nodes[k]) / self.scale;
            if d.abs() <= 8.0 * f64::EPSILON {
                return Ok(self.values[k]);
            }
            let term = self.weights[k] / d;
            numerator = numerator + self.values[k] * term;
            denominator += term;
        }
        Ok(numerator / denominator)
    }
}
