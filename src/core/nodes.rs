//! Chebyshev node families.

use super::{ChebyTime, Domain};

/// Supported Chebyshev node families.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum NodeKind {
    /// Roots of `T_N`.
    Roots,
    /// Extrema of `T_{N-1}`.
    Extrema,
    /// Alias for [`NodeKind::Extrema`].
    Lobatto,
    /// Alias for [`NodeKind::Roots`].
    Gauss,
    /// Alias for [`NodeKind::Extrema`].
    GaussLobatto,
}

/// Return `N` nodes in `[-1, 1]`.
#[inline]
pub fn nodes<const N: usize>(kind: NodeKind) -> [f64; N] {
    match kind {
        NodeKind::Roots | NodeKind::Gauss => {
            let n = N as f64;
            core::array::from_fn(|k| {
                (core::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * n)).cos()
            })
        }
        NodeKind::Extrema | NodeKind::Lobatto | NodeKind::GaussLobatto => {
            if N == 1 {
                [0.0; N]
            } else {
                let denom = (N - 1) as f64;
                core::array::from_fn(|k| (core::f64::consts::PI * k as f64 / denom).cos())
            }
        }
    }
}

/// Return `N` nodes mapped from `[-1, 1]` into `domain`.
#[inline]
pub fn nodes_mapped<X: ChebyTime, const N: usize>(domain: Domain<X>, kind: NodeKind) -> [X; N] {
    let unit = nodes::<N>(kind);
    core::array::from_fn(|k| domain.denormalize(unit[k]))
}
