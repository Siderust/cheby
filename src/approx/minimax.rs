//! Minimax (Remez exchange) approximation on `f64`.
//!
//! # Algorithm
//!
//! For a continuous `f: [a, b] -> R` and a target degree `N`, this module
//! computes a polynomial `p(x) = sum_{k=0}^{N} c_k T_k(tau)` (with
//! `tau = domain.normalize(x)`) that approximately minimizes
//! `max_{x in [a,b]} |f(x) - p(x)|`.
//!
//! Each iteration:
//!
//! 1. Solve the `(N+2) x (N+2)` linear system
//!    `sum_k c_k T_k(tau_i) + (-1)^i E = f(x_i)` over the current
//!    alternation set `x_0 < x_1 < ... < x_{N+1}` (Gaussian elimination
//!    with partial pivoting).
//! 2. Sample the residual `r(x) = f(x) - p(x)` on a dense grid,
//!    partition the grid by sign of `r`, and locate the per-sign-block
//!    maxima of `|r|`.
//! 3. If fewer than `N+2` alternations are available the iteration
//!    halts with [`ChebyError::RemezAlternationFailure`]. Otherwise the
//!    best `N+2` contiguous alternation points become the next set.
//! 4. Convergence is declared when `(max|r| - |E|) / max|r|` falls below
//!    `RemezOptions::tolerance`. The iteration budget is
//!    `RemezOptions::max_iterations`; exhausting it returns
//!    [`ChebyError::RemezDidNotConverge`].
//!
//! Root-sample (interpolating) coefficient fitting is **not** minimax;
//! it lives in [`crate::approx::fit`].
//!
//! # Performance notes
//!
//! Cost per iteration is `O((N+2)^3 + grid_size * N)`. Pick `degree` and
//! `grid_size` accordingly; `grid_size >= 16 * (degree + 2)` is usually
//! enough for smooth targets.

use alloc::vec;
use alloc::vec::Vec;

use crate::core::{
    basis, ChebyError, ChebySeries, ChebySeriesDyn, ChebySeriesOn, ChebyTime, Domain,
};

/// Options for the Remez exchange algorithm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemezOptions {
    /// Maximum exchange iterations.
    pub max_iterations: usize,
    /// Relative convergence tolerance on the residual equioscillation gap.
    pub tolerance: f64,
    /// Dense residual sampling grid size used to locate alternations.
    pub grid_size: usize,
}

impl Default for RemezOptions {
    fn default() -> Self {
        Self {
            max_iterations: 32,
            tolerance: 1e-12,
            grid_size: 1024,
        }
    }
}

/// Result of a successful Remez exchange.
#[derive(Debug, Clone, PartialEq)]
pub struct RemezResult<X> {
    /// Domain-aware fitted series.
    pub series_on: ChebySeriesOnAny<X>,
    /// Underlying coefficient view.
    pub series: ChebySeriesDyn<f64>,
    /// Final equioscillation level returned by the linear solve (`|E|`).
    pub leveled_error: f64,
    /// Final dense-grid sup-norm error `max |f - p|`.
    pub max_error: f64,
    /// Iterations actually performed.
    pub iterations: usize,
    /// Whether convergence tolerance was met.
    pub converged: bool,
    /// Final alternation points in domain coordinates.
    pub alternations: Vec<X>,
    /// Domain used for fitting.
    pub domain: Domain<X>,
}

/// Domain-aware Remez series wrapper exposed for convenience evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct ChebySeriesOnAny<X> {
    domain: Domain<X>,
    coeffs: Vec<f64>,
}

impl<X: ChebyTime> ChebySeriesOnAny<X> {
    /// Domain of validity.
    #[inline]
    pub fn domain(&self) -> Domain<X> {
        self.domain
    }
    /// Coefficient slice in the first-kind Chebyshev basis.
    #[inline]
    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }
    /// Evaluate at a domain point.
    #[inline]
    pub fn evaluate(&self, x: X) -> Result<f64, ChebyError> {
        if !self.domain.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        Ok(crate::core::eval::evaluate(
            &self.coeffs,
            self.domain.normalize(x),
        ))
    }
}

/// Backward-compatible fixed-degree convenience wrapper.
///
/// `Convert` a [`RemezResult`] coefficient slice into a fixed-size
/// [`ChebySeries`] when the degree is known statically.
pub fn try_into_fixed<const N: usize>(coeffs: &[f64]) -> Result<ChebySeries<f64, N>, ChebyError> {
    if coeffs.len() != N {
        return Err(ChebyError::InvalidDegree);
    }
    let mut arr = [0.0; N];
    arr.copy_from_slice(coeffs);
    Ok(ChebySeries::new(arr))
}

/// Run the Remez exchange algorithm for an `f64`-valued target.
///
/// The polynomial degree is `degree`; the resulting series has
/// `degree + 1` coefficients in the first-kind Chebyshev basis.
pub fn remez<X>(
    domain: Domain<X>,
    degree: usize,
    f: impl Fn(X) -> f64,
    options: RemezOptions,
) -> Result<RemezResult<X>, ChebyError>
where
    X: ChebyTime,
{
    if degree == 0 {
        return Err(ChebyError::InvalidDegree);
    }
    let max_iter = options.max_iterations.max(1);
    let m = degree + 2; // alternation count
    let grid = options.grid_size.max(8 * m);
    if !(options.tolerance.is_finite() && options.tolerance >= 0.0) {
        return Err(ChebyError::NonFiniteInput);
    }

    // Initial alternation set: Chebyshev extrema (Lobatto) of degree N+1
    // mapped into the domain, ascending.
    let mut taus: Vec<f64> = (0..m)
        .map(|k| -(core::f64::consts::PI * k as f64 / (m - 1) as f64).cos())
        .collect();

    // Pre-allocated work buffers.
    let n = degree + 1;
    let mut matrix = vec![0.0; m * (m + 1)];
    let mut coeffs = vec![0.0; n];
    let mut prev_max_err = f64::INFINITY;
    let mut max_err = 0.0;
    let mut leveled = 0.0;
    let mut iterations = 0;
    let mut converged = false;

    for iter in 1..=max_iter {
        iterations = iter;

        // Build linear system: row i = [T_0(tau_i) ... T_N(tau_i) | (-1)^i | f(x_i)]
        for (i, &tau) in taus.iter().enumerate() {
            let row = i * (m + 1);
            for k in 0..n {
                matrix[row + k] = basis::t(k, tau);
            }
            matrix[row + n] = if i % 2 == 0 { 1.0 } else { -1.0 };
            matrix[row + m] = f(domain.denormalize(tau));
        }

        if !solve_in_place(&mut matrix, m) {
            return Err(ChebyError::RemezDidNotConverge);
        }

        for k in 0..n {
            coeffs[k] = matrix[k * (m + 1) + m];
        }
        leveled = matrix[n * (m + 1) + m].abs();

        // Sample residual on dense grid.
        let mut residual = vec![0.0_f64; grid];
        for (i, r) in residual.iter_mut().enumerate() {
            let tau = -1.0 + 2.0 * i as f64 / (grid - 1) as f64;
            let x = domain.denormalize(tau);
            *r = f(x) - crate::core::eval::evaluate(&coeffs, tau);
        }

        // Locate per-sign-block extrema.
        let mut blocks: Vec<(f64, f64)> = Vec::new(); // (tau_max, residual_max_signed)
        let mut cur_sign = 0i8;
        let mut cur_best_idx = 0usize;
        let mut cur_best_abs = -1.0_f64;
        for (i, &r) in residual.iter().enumerate() {
            let s = if r > 0.0 {
                1
            } else if r < 0.0 {
                -1
            } else {
                0
            };
            if s == 0 {
                continue;
            }
            if cur_sign == 0 {
                cur_sign = s;
                cur_best_idx = i;
                cur_best_abs = r.abs();
            } else if s == cur_sign {
                if r.abs() > cur_best_abs {
                    cur_best_abs = r.abs();
                    cur_best_idx = i;
                }
            } else {
                let tau = -1.0 + 2.0 * cur_best_idx as f64 / (grid - 1) as f64;
                blocks.push((tau, residual[cur_best_idx]));
                cur_sign = s;
                cur_best_idx = i;
                cur_best_abs = r.abs();
            }
        }
        if cur_sign != 0 {
            let tau = -1.0 + 2.0 * cur_best_idx as f64 / (grid - 1) as f64;
            blocks.push((tau, residual[cur_best_idx]));
        }

        if blocks.len() < m {
            return Err(ChebyError::RemezAlternationFailure);
        }

        // Pick the contiguous window of m blocks maximizing min |r|.
        let mut best_start = 0usize;
        let mut best_min = -1.0_f64;
        for s in 0..=blocks.len() - m {
            let window_min = blocks
                .iter()
                .skip(s)
                .take(m)
                .map(|b| b.1.abs())
                .fold(f64::INFINITY, f64::min);
            if window_min > best_min {
                best_min = window_min;
                best_start = s;
            }
        }

        let mut new_taus: Vec<f64> = blocks[best_start..best_start + m]
            .iter()
            .map(|&(t, _)| t)
            .collect();
        new_taus.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));

        max_err = blocks[best_start..best_start + m]
            .iter()
            .map(|&(_, r)| r.abs())
            .fold(0.0_f64, f64::max);

        let gap = (max_err - leveled).abs();
        let denom = max_err.max(f64::EPSILON);
        if gap / denom <= options.tolerance {
            converged = true;
            taus = new_taus;
            break;
        }

        // Stagnation guard.
        if (prev_max_err - max_err).abs() <= options.tolerance * denom && iter > 1 {
            taus = new_taus;
            converged = true;
            break;
        }
        prev_max_err = max_err;
        taus = new_taus;
    }

    if !converged {
        return Err(ChebyError::RemezDidNotConverge);
    }

    let series = ChebySeriesDyn::new(coeffs.clone())?;
    Ok(RemezResult {
        series_on: ChebySeriesOnAny {
            domain,
            coeffs: coeffs.clone(),
        },
        series,
        leveled_error: leveled,
        max_error: max_err,
        iterations,
        converged,
        alternations: taus.into_iter().map(|t| domain.denormalize(t)).collect(),
        domain,
    })
}

/// Augmented Gaussian elimination with partial pivoting.
///
/// Matrix layout: `m` rows, `m + 1` columns (last column is RHS).
/// On success the RHS column holds the solution.
fn solve_in_place(a: &mut [f64], m: usize) -> bool {
    let stride = m + 1;
    for k in 0..m {
        // Pivot.
        let mut piv_row = k;
        let mut piv_abs = a[k * stride + k].abs();
        for r in (k + 1)..m {
            let v = a[r * stride + k].abs();
            if v > piv_abs {
                piv_abs = v;
                piv_row = r;
            }
        }
        if piv_abs < 1e-300 {
            return false;
        }
        if piv_row != k {
            for c in 0..stride {
                a.swap(k * stride + c, piv_row * stride + c);
            }
        }
        // Eliminate.
        let pivot = a[k * stride + k];
        for r in (k + 1)..m {
            let factor = a[r * stride + k] / pivot;
            if factor == 0.0 {
                continue;
            }
            for c in k..stride {
                a[r * stride + c] -= factor * a[k * stride + c];
            }
        }
    }
    // Back-substitute, store solution in RHS column (`stride - 1`).
    for r in (0..m).rev() {
        let mut sum = a[r * stride + (stride - 1)];
        for c in (r + 1)..m {
            sum -= a[r * stride + c] * a[c * stride + (stride - 1)];
        }
        a[r * stride + (stride - 1)] = sum / a[r * stride + r];
    }
    true
}

/// Convert a Remez result into a fixed-size [`ChebySeriesOn`] when the
/// caller knows the degree at compile time.
pub fn into_fixed_series_on<X: ChebyTime, const N: usize>(
    result: &RemezResult<X>,
) -> Result<ChebySeriesOn<f64, X, N>, ChebyError> {
    let series = try_into_fixed::<N>(result.series.coeffs())?;
    Ok(ChebySeriesOn::new(result.domain, series))
}
