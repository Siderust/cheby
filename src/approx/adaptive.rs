//! Adaptive Chebyshev fitting.

use crate::approx::{error_estimation, fit};
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

macro_rules! build_degree {
    ($n:literal, $self:ident) => {{
        let fitted = fit::fit_from_fn_on::<T, X, $n>($self.domain, &$self.f);
        let coeffs = fitted.into_coeffs();
        let err = error_estimation::estimated_truncation_error(&coeffs);
        let converged = error_estimation::is_converged(&coeffs, $self.abs_tol, $self.rel_tol);
        let series = ChebySeriesDyn::new(alloc::vec::Vec::from(coeffs))?;
        return Ok(AdaptiveFitResult {
            series,
            report: FitReport {
                degree: $n,
                estimated_error: err,
                converged,
            },
        });
    }};
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
        for n in [8_usize, 16, 32, 64, 128] {
            if n > self.max_degree {
                break;
            }
            match n {
                8 => {
                    let fitted = fit::fit_from_fn_on::<T, X, 8>(self.domain, &self.f);
                    let coeffs = fitted.into_coeffs();
                    if error_estimation::is_converged(&coeffs, self.abs_tol, self.rel_tol)
                        || n == self.max_degree
                    {
                        let err = error_estimation::estimated_truncation_error(&coeffs);
                        return Ok(AdaptiveFitResult {
                            series: ChebySeriesDyn::new(alloc::vec::Vec::from(coeffs))?,
                            report: FitReport {
                                degree: n,
                                estimated_error: err,
                                converged: error_estimation::is_converged(
                                    &coeffs,
                                    self.abs_tol,
                                    self.rel_tol,
                                ),
                            },
                        });
                    }
                }
                16 => {
                    let fitted = fit::fit_from_fn_on::<T, X, 16>(self.domain, &self.f);
                    let coeffs = fitted.into_coeffs();
                    if error_estimation::is_converged(&coeffs, self.abs_tol, self.rel_tol)
                        || n == self.max_degree
                    {
                        let err = error_estimation::estimated_truncation_error(&coeffs);
                        return Ok(AdaptiveFitResult {
                            series: ChebySeriesDyn::new(alloc::vec::Vec::from(coeffs))?,
                            report: FitReport {
                                degree: n,
                                estimated_error: err,
                                converged: error_estimation::is_converged(
                                    &coeffs,
                                    self.abs_tol,
                                    self.rel_tol,
                                ),
                            },
                        });
                    }
                }
                32 => {
                    let fitted = fit::fit_from_fn_on::<T, X, 32>(self.domain, &self.f);
                    let coeffs = fitted.into_coeffs();
                    if error_estimation::is_converged(&coeffs, self.abs_tol, self.rel_tol)
                        || n == self.max_degree
                    {
                        let err = error_estimation::estimated_truncation_error(&coeffs);
                        return Ok(AdaptiveFitResult {
                            series: ChebySeriesDyn::new(alloc::vec::Vec::from(coeffs))?,
                            report: FitReport {
                                degree: n,
                                estimated_error: err,
                                converged: error_estimation::is_converged(
                                    &coeffs,
                                    self.abs_tol,
                                    self.rel_tol,
                                ),
                            },
                        });
                    }
                }
                64 => build_degree!(64, self),
                128 => build_degree!(128, self),
                _ => {}
            }
        }
        Err(ChebyError::InvalidDegree)
    }
}
