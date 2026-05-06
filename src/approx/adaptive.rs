//! Adaptive Chebyshev fitting.

use alloc::vec;
use alloc::vec::Vec;

use crate::approx::error_estimation;
use crate::core::{ChebyError, ChebyScalar, ChebySeriesDyn, ChebyTime, Domain};

/// Metadata returned by adaptive fitting.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FitReport {
    /// Number of coefficients in the returned series.
    pub degree: usize,
    /// Tail-based error estimate.
    pub estimated_error: f64,
    /// Whether tolerance criteria were satisfied.
    pub converged: bool,
}

/// Result of an adaptive fit.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveFitResult<T> {
    /// Fitted dynamic series.
    pub series: ChebySeriesDyn<T>,
    /// Fit metadata.
    pub report: FitReport,
}

/// Builder for adaptive fitting.
pub struct AdaptiveFit<F, X> {
    domain: Domain<X>,
    f: F,
    max_degree: usize,
    abs_tol: f64,
    rel_tol: f64,
}

impl<F, X> AdaptiveFit<F, X> {
    /// Start an adaptive fit builder.
    pub fn new(domain: Domain<X>, f: F) -> Self {
        Self {
            domain,
            f,
            max_degree: 64,
            abs_tol: 1e-12,
            rel_tol: 1e-12,
        }
    }

    /// Set maximum coefficient count.
    pub fn max_degree(mut self, max_degree: usize) -> Self {
        self.max_degree = max_degree;
        self
    }

    /// Set absolute tolerance.
    pub fn absolute_tolerance(mut self, tol: f64) -> Self {
        self.abs_tol = tol;
        self
    }

    /// Set relative tolerance.
    pub fn relative_tolerance(mut self, tol: f64) -> Self {
        self.rel_tol = tol;
        self
    }
}

impl<F, X> AdaptiveFit<F, X>
where
    X: ChebyTime,
{
    /// Build the adaptive fit.
    pub fn build<T>(self) -> Result<AdaptiveFitResult<T>, ChebyError>
    where
        T: ChebyScalar,
        F: Fn(X) -> T,
    {
        if self.max_degree < 2 {
            return Err(ChebyError::InvalidDegree);
        }
        if !(self.abs_tol.is_finite()
            && self.abs_tol >= 0.0
            && self.rel_tol.is_finite()
            && self.rel_tol >= 0.0)
        {
            return Err(ChebyError::NonFiniteInput);
        }

        let degrees = degree_schedule(self.max_degree);
        for n in degrees {
            let coeffs = fit_dyn_from_fn_on::<T, X, F>(self.domain, n, &self.f);
            let err = error_estimation::estimated_truncation_error(&coeffs);
            let converged = error_estimation::is_converged(&coeffs, self.abs_tol, self.rel_tol);
            if converged || n == self.max_degree {
                return Ok(AdaptiveFitResult {
                    series: ChebySeriesDyn::new(coeffs)?,
                    report: FitReport {
                        degree: n,
                        estimated_error: err,
                        converged,
                    },
                });
            }
        }
        Err(ChebyError::InvalidDegree)
    }
}

fn degree_schedule(max_degree: usize) -> Vec<usize> {
    if max_degree < 8 {
        return vec![max_degree];
    }

    let mut degrees = Vec::new();
    let mut n = 8_usize;
    while n < max_degree {
        degrees.push(n);
        match n.checked_mul(2) {
            Some(next) => n = next,
            None => break,
        }
    }
    if degrees.last().copied() != Some(max_degree) {
        degrees.push(max_degree);
    }
    degrees
}

fn fit_dyn_from_fn_on<T, X, F>(domain: Domain<X>, n: usize, f: &F) -> Vec<T>
where
    T: ChebyScalar,
    X: ChebyTime,
    F: Fn(X) -> T,
{
    let nf = n as f64;
    // Sample at roots of T_n; cache the normalized node `x_k = cos(theta_k)`
    // so we can build T_j(x_k) by Chebyshev recurrence in `j` instead of
    // calling `cos` n*n times.
    let mut xs = vec![0.0_f64; n];
    let mut values: Vec<T> = Vec::with_capacity(n);
    for (k, x_slot) in xs.iter_mut().enumerate() {
        let x = (core::f64::consts::PI * (2.0 * k as f64 + 1.0) / (2.0 * nf)).cos();
        *x_slot = x;
        values.push(f(domain.denormalize(x)));
    }

    let mut coeffs = vec![T::zero(); n];
    for (k, &v) in values.iter().enumerate() {
        let x = xs[k];
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
    coeffs[0] = coeffs[0] / nf;
    let scale = 2.0 / nf;
    for c in coeffs.iter_mut().skip(1) {
        *c = *c * scale;
    }
    coeffs
}
