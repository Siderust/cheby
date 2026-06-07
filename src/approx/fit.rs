//! Coefficient fitting from samples and functions.
//!
//! The fixed-size fitter samples at roots of `T_N` and applies the direct
//! discrete-cosine coefficient formula. With `T_j(x_k) = cos(j θ_k)` evaluated
//! by Chebyshev recurrence in `j`, the work is `O(N²)` multiplies-and-adds
//! and only `O(N)` calls to `cos` (used to build the node grid once).

use crate::core::{nodes, ChebyError, ChebyScalar, ChebySeries, ChebySeriesOn, Domain, NodeKind};

#[cfg(feature = "alloc")]
use crate::core::ChebySeriesDyn;

/// Compute Chebyshev coefficients from values sampled at roots of `T_N`.
///
/// Uses the Chebyshev three-term recurrence to generate `T_j(x_k)` for each
/// sample `k` in a single pass over `j`, avoiding `N²` `cos()` calls.
#[inline]
pub fn fit_coeffs<T: ChebyScalar, const N: usize>(values: &[T; N]) -> [T; N] {
    let xs = nodes::<N>(NodeKind::Roots);
    let mut coeffs = [T::zero(); N];
    accumulate_coeffs(values.as_slice(), &xs, &mut coeffs);
    coeffs
}

/// Compute dynamic Chebyshev coefficients from samples at supplied unit nodes.
#[cfg(feature = "alloc")]
pub fn fit_coeffs_dyn(values: &[f64], unit_nodes: &[f64]) -> Vec<f64> {
    let mut coeffs = alloc::vec![0.0_f64; values.len()];
    accumulate_coeffs(values, unit_nodes, &mut coeffs);
    coeffs
}

fn accumulate_coeffs<T: ChebyScalar>(values: &[T], unit_nodes: &[f64], coeffs: &mut [T]) {
    let n = values.len();
    if n == 0 || coeffs.len() != n {
        return;
    }
    debug_assert_eq!(n, unit_nodes.len());

    for coeff in coeffs.iter_mut() {
        *coeff = T::zero();
    }

    for (k, &v) in values.iter().enumerate() {
        let x = unit_nodes[k];
        coeffs[0] = coeffs[0] + v;
        if n > 1 {
            coeffs[1] = coeffs[1] + v * x;
            let two_x = 2.0 * x;
            let mut tkm1 = 1.0_f64;
            let mut tk = x;
            for coeff in coeffs.iter_mut().take(n).skip(2) {
                let tkp1 = two_x * tk - tkm1;
                *coeff = *coeff + v * tkp1;
                tkm1 = tk;
                tk = tkp1;
            }
        }
    }

    let nf = n as f64;
    coeffs[0] = coeffs[0] / nf;
    let scale = 2.0 / nf;
    for coeff in coeffs.iter_mut().skip(1) {
        *coeff = *coeff * scale;
    }
}

/// Fit a Chebyshev series on normalized `[-1, 1]`.
#[inline]
pub fn fit_from_fn<T: ChebyScalar, const N: usize>(f: impl Fn(f64) -> T) -> ChebySeries<T, N> {
    let xs = nodes::<N>(NodeKind::Roots);
    let values = core::array::from_fn(|k| f(xs[k]));
    ChebySeries::new(fit_coeffs(&values))
}

/// Fit a dynamic Chebyshev series of `degree` on normalized `[-1, 1]`.
#[cfg(feature = "alloc")]
pub fn fit_dyn_from_fn(
    degree: usize,
    mut f: impl FnMut(f64) -> f64,
) -> Result<ChebySeriesDyn<f64>, ChebyError> {
    let n = degree.saturating_add(1).max(2);
    let unit_nodes = unit_root_nodes(n);
    let values: Vec<f64> = unit_nodes.iter().map(|&x| f(x)).collect();
    let coeffs = fit_coeffs_dyn(&values, &unit_nodes);
    ChebySeriesDyn::new(coeffs)
}

#[cfg(feature = "alloc")]
fn unit_root_nodes(n: usize) -> Vec<f64> {
    let nf = n as f64;
    (0..n)
        .map(|k| (core::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * nf)).cos())
        .collect()
}

/// Fit a Chebyshev series on a physical domain.
#[inline]
pub fn fit_from_fn_on<T, X, const N: usize>(
    domain: Domain<X>,
    f: impl Fn(X) -> T,
) -> ChebySeriesOn<T, X, N>
where
    T: ChebyScalar,
    X: crate::core::ChebyTime,
{
    let xs = crate::core::nodes::nodes_mapped::<X, N>(domain, NodeKind::Roots);
    let values = core::array::from_fn(|k| f(xs[k]));
    ChebySeriesOn::new(domain, ChebySeries::new(fit_coeffs(&values)))
}
