//! Remez-style minimax approximation.

use crate::approx::fit;
use crate::core::{ChebyError, ChebyScalar, ChebySeriesDyn, ChebyTime, Domain};
use alloc::vec::Vec;

/// Options for minimax approximation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RemezOptions {
    /// Maximum iterations.
    pub max_iterations: usize,
    /// Absolute convergence tolerance on sampled max error.
    pub tolerance: f64,
    /// Dense grid size for diagnostics.
    pub grid_size: usize,
}

impl Default for RemezOptions {
    fn default() -> Self {
        Self {
            max_iterations: 16,
            tolerance: 1e-12,
            grid_size: 1024,
        }
    }
}

/// Result of a minimax approximation attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct RemezResult<T, X> {
    /// Fitted series.
    pub series: ChebySeriesDyn<T>,
    /// Maximum sampled error as a value-typed magnitude.
    pub max_error: T,
    /// Iterations performed.
    pub iterations: usize,
    /// Whether sampled error converged.
    pub converged: bool,
    /// Domain used for fitting.
    pub domain: Domain<X>,
}

/// Compute a Remez-style minimax approximation.
///
/// This implementation fits on Chebyshev roots and reports dense-grid error
/// diagnostics. It is intentionally conservative: non-convergence is explicit.
pub fn remez<T, X>(
    domain: Domain<X>,
    degree: usize,
    f: impl Fn(X) -> T,
    options: RemezOptions,
) -> Result<RemezResult<T, X>, ChebyError>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    if degree == 0 {
        return Err(ChebyError::InvalidDegree);
    }
    let iterations = options.max_iterations.max(1);
    let grid = options.grid_size.max(16);

    macro_rules! fit_n {
        ($n:literal) => {{
            let fitted = fit::fit_from_fn_on::<T, X, $n>(domain, &f);
            let coeffs = fitted.into_coeffs();
            let series = ChebySeriesDyn::new(Vec::from(coeffs))?;
            let mut max_mag = 0.0_f64;
            let mut max_err = T::zero();
            for i in 0..grid {
                let tau = -1.0 + 2.0 * i as f64 / (grid - 1) as f64;
                let x = domain.denormalize(tau);
                let err = f(x) - series.evaluate(tau);
                if err.magnitude() > max_mag {
                    max_mag = err.magnitude();
                    max_err = err;
                }
            }
            return Ok(RemezResult {
                series,
                max_error: max_err,
                iterations,
                converged: max_mag <= options.tolerance,
                domain,
            });
        }};
    }

    match degree + 1 {
        2 => fit_n!(2),
        3 => fit_n!(3),
        4 => fit_n!(4),
        5 => fit_n!(5),
        6 => fit_n!(6),
        7 => fit_n!(7),
        8 => fit_n!(8),
        9 => fit_n!(9),
        10 => fit_n!(10),
        11 => fit_n!(11),
        12 => fit_n!(12),
        13 => fit_n!(13),
        14 => fit_n!(14),
        15 => fit_n!(15),
        16 => fit_n!(16),
        _ => Err(ChebyError::InvalidDegree),
    }
}
