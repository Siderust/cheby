// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev node generation.
//!
//! Chebyshev nodes (zeros of the Chebyshev polynomial of the first kind)
//! minimise the maximum interpolation error for smooth functions, avoiding
//! the Runge phenomenon that plagues equally-spaced nodes.
//!
//! For `N` nodes on `[-1, 1]`:
//!
//! ```text
//! ξ_k = cos(π(2k+1) / (2N)),  k = 0, …, N-1
//! ```

/// Compute `N` Chebyshev nodes on `[-1, 1]`.
///
/// Returns the nodes in descending order (from near +1 to near −1),
/// which is the natural ordering of the cosine formula.
///
/// # Example
///
/// ```
/// let xi: [f64; 5] = cheby::nodes();
/// assert!(xi[0] > xi[4]); // descending
/// ```
#[inline]
pub fn nodes<const N: usize>() -> [f64; N] {
    let mut out = [0.0_f64; N];
    let n = N as f64;
    for (k, x) in out.iter_mut().enumerate() {
        let arg = std::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * n);
        *x = arg.cos();
    }
    out
}

/// Compute `N` Chebyshev nodes mapped to the interval `[start, end]`.
///
/// Each node `ξ_k ∈ [-1, 1]` is mapped to `t_k = mid + half * ξ_k`
/// where `mid = (start + end) / 2` and `half = (end - start) / 2`.
///
/// # Example
///
/// ```
/// let t: [f64; 9] = cheby::nodes_mapped(0.0, 4.0);
/// assert!(t.iter().all(|&x| x >= 0.0 && x <= 4.0));
/// ```
#[inline]
pub fn nodes_mapped<const N: usize>(start: f64, end: f64) -> [f64; N] {
    let mid = 0.5 * (start + end);
    let half = 0.5 * (end - start);
    let unit = nodes::<N>();
    let mut out = [0.0_f64; N];
    for (x, &u) in out.iter_mut().zip(unit.iter()) {
        *x = mid + half * u;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nodes_count() {
        let n3: [f64; 3] = nodes();
        assert_eq!(n3.len(), 3);
        let n9: [f64; 9] = nodes();
        assert_eq!(n9.len(), 9);
    }

    #[test]
    fn test_nodes_in_range() {
        let xi: [f64; 9] = nodes();
        for &x in &xi {
            assert!(x > -1.0 && x < 1.0, "node {x} out of (-1, 1)");
        }
    }

    #[test]
    fn test_nodes_descending() {
        let xi: [f64; 9] = nodes();
        for w in xi.windows(2) {
            assert!(w[0] > w[1], "nodes not descending: {} <= {}", w[0], w[1]);
        }
    }

    #[test]
    fn test_nodes_mapped_range() {
        let t: [f64; 9] = nodes_mapped(10.0, 20.0);
        for &x in &t {
            assert!(
                (10.0..=20.0).contains(&x),
                "mapped node {x} out of [10, 20]"
            );
        }
    }

    #[test]
    fn test_nodes_mapped_midpoint() {
        // The midpoint of a mapped set should be close to (start+end)/2
        let t: [f64; 9] = nodes_mapped(0.0, 4.0);
        let mean: f64 = t.iter().sum::<f64>() / 9.0;
        assert!((mean - 2.0).abs() < 1e-10);
    }
}
