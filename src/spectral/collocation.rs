//! Collocation points.

use crate::core::NodeKind;

/// Chebyshev-Lobatto collocation points.
#[inline]
pub fn collocation_points<const N: usize>() -> [f64; N] {
    crate::core::nodes::nodes(NodeKind::GaussLobatto)
}
