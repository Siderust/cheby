//! Barycentric interpolation.

use crate::core::{ChebyError, ChebyScalar, ChebyTime};

/// Barycentric interpolator over supplied nodes and values.
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
    pub fn new(nodes: [X; N], values: [T; N]) -> Result<Self, ChebyError> {
        if N == 0 {
            return Err(ChebyError::InvalidDegree);
        }
        let scale = if N > 1 {
            nodes[0] - nodes[N - 1]
        } else {
            X::zero() + nodes[0] * 0.0 + nodes[0] * 1.0
        };
        if N > 1 && scale == X::zero() {
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

    /// Evaluate the interpolant at `x`.
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
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
