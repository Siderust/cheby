//! Fixed-size and dynamic Chebyshev series.

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

use super::{eval, ChebyError, ChebyScalar, ChebyTime, DifferentiateWith, Domain};

/// A fixed-size Chebyshev coefficient series on normalized `[-1, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChebySeries<T, const N: usize> {
    coeffs: [T; N],
}

#[cfg(feature = "serde")]
impl<T, const N: usize> serde::Serialize for ChebySeries<T, N>
where
    T: serde::Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.coeffs.as_slice().serialize(serializer)
    }
}

#[cfg(feature = "serde")]
impl<'de, T, const N: usize> serde::Deserialize<'de> for ChebySeries<T, N>
where
    T: serde::Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let coeffs = <alloc::vec::Vec<T>>::deserialize(deserializer)?;
        let coeffs: [T; N] = coeffs.try_into().map_err(|v: alloc::vec::Vec<T>| {
            serde::de::Error::invalid_length(v.len(), &"matching ChebySeries coefficient count")
        })?;
        Ok(Self { coeffs })
    }
}

impl<T: ChebyScalar, const N: usize> ChebySeries<T, N> {
    /// Create a series from coefficients in the first-kind Chebyshev basis.
    #[inline]
    pub fn new(coeffs: [T; N]) -> Self {
        Self { coeffs }
    }

    /// Borrow the coefficient array.
    #[inline]
    pub fn coeffs(&self) -> &[T; N] {
        &self.coeffs
    }

    /// Consume the series and return its coefficients.
    #[inline]
    pub fn into_coeffs(self) -> [T; N] {
        self.coeffs
    }

    /// Evaluate at normalized coordinate `x`.
    #[inline]
    pub fn evaluate(&self, x: f64) -> T {
        eval::evaluate(&self.coeffs, x)
    }

    /// Evaluate both value and normalized derivative.
    #[inline]
    pub fn evaluate_both(&self, x: f64) -> (T, T) {
        eval::evaluate_both(&self.coeffs, x)
    }

    /// Return a normalized-coordinate derivative series.
    ///
    /// The array length is preserved; unused high-order entries are zero.
    pub fn derivative(&self) -> Self {
        let mut out = [T::zero(); N];
        if N <= 1 {
            return Self::new(out);
        }
        out[N - 2] = self.coeffs[N - 1] * (2.0 * (N - 1) as f64);
        if N > 2 {
            for k in (1..=(N - 2)).rev() {
                out[k - 1] = out[k + 1] + self.coeffs[k] * (2.0 * k as f64);
            }
        }
        out[0] = out[0] / 2.0;
        Self::new(out)
    }

    /// Return a normalized-coordinate indefinite integral series.
    ///
    /// The array length is preserved; the single highest antiderivative term
    /// (`b_N`) is truncated. For a non-truncating dynamic version see
    /// [`ChebySeriesDyn::integral`].
    pub fn integral(&self, constant: T) -> Self {
        let mut out = [T::zero(); N];
        if N == 0 {
            return Self::new(out);
        }
        out[0] = constant;
        if N > 1 {
            out[1] = self.coeffs[0];
        }
        if N > 2 {
            out[2] = out[2] + self.coeffs[1] / 4.0;
        }
        // For k >= 2: ∫T_k = (T_{k+1}/(k+1) - T_{k-1}/(k-1)) / 2.
        // Apply both contributions; only `b_{k+1}` for k = N-1 is dropped.
        for k in 2..N {
            if k + 1 < N {
                out[k + 1] = out[k + 1] + self.coeffs[k] / (2.0 * (k + 1) as f64);
            }
            out[k - 1] = out[k - 1] - self.coeffs[k] / (2.0 * (k - 1) as f64);
        }
        Self::new(out)
    }
}

/// A fixed-size Chebyshev series tied to a physical domain.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChebySeriesOn<T, X, const N: usize> {
    domain: Domain<X>,
    series: ChebySeries<T, N>,
}

impl<T: ChebyScalar, X: ChebyTime, const N: usize> ChebySeriesOn<T, X, N> {
    /// Create a domain-aware series.
    #[inline]
    pub fn new(domain: Domain<X>, series: ChebySeries<T, N>) -> Self {
        Self { domain, series }
    }

    /// Domain of validity.
    #[inline]
    pub fn domain(&self) -> Domain<X> {
        self.domain
    }

    /// Underlying normalized-coordinate series.
    #[inline]
    pub fn series(&self) -> &ChebySeries<T, N> {
        &self.series
    }

    /// Consume and return raw coefficients.
    #[inline]
    pub fn into_coeffs(self) -> [T; N] {
        self.series.into_coeffs()
    }

    /// Evaluate inside the domain.
    #[inline]
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
        if !self.domain.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        Ok(self.series.evaluate(self.domain.normalize(x)))
    }

    /// Evaluate the dimensionally-correct physical derivative.
    #[inline]
    pub fn evaluate_derivative(
        &self,
        x: X,
    ) -> Result<<T as DifferentiateWith<X>>::Derivative, ChebyError>
    where
        T: DifferentiateWith<X>,
    {
        if !self.domain.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        let d_tau = self.series.evaluate_both(self.domain.normalize(x)).1;
        Ok((d_tau * 2.0).scale_derivative(self.domain.end() - self.domain.start()))
    }
}

/// A dynamic Chebyshev series.
#[cfg(feature = "alloc")]
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChebySeriesDyn<T> {
    coeffs: Vec<T>,
}

#[cfg(feature = "alloc")]
impl<T: ChebyScalar> ChebySeriesDyn<T> {
    /// Create a dynamic series from coefficients.
    #[inline]
    pub fn new(coeffs: Vec<T>) -> Result<Self, ChebyError> {
        if coeffs.is_empty() {
            return Err(ChebyError::EmptyCoefficientSet);
        }
        Ok(Self { coeffs })
    }

    /// Borrow coefficients.
    #[inline]
    pub fn coeffs(&self) -> &[T] {
        &self.coeffs
    }

    /// Number of coefficients (degree + 1).
    #[inline]
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Whether the series has no coefficients (always `false` once
    /// constructed via [`Self::new`]; included for API completeness).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Consume into coefficients.
    #[inline]
    pub fn into_coeffs(self) -> Vec<T> {
        self.coeffs
    }

    /// Evaluate at normalized coordinate `x`.
    #[inline]
    pub fn evaluate(&self, x: f64) -> T {
        eval::evaluate(&self.coeffs, x)
    }

    /// Evaluate value and normalized derivative `df/dtau` together.
    #[inline]
    pub fn evaluate_both(&self, x: f64) -> (T, T) {
        eval::evaluate_both(&self.coeffs, x)
    }

    /// Evaluate the normalized-coordinate derivative `df/dtau`.
    #[inline]
    pub fn evaluate_derivative(&self, x: f64) -> T {
        eval::evaluate_derivative(&self.coeffs, x)
    }

    /// Return a normalized-coordinate derivative series.
    ///
    /// The output has `max(len() - 1, 1)` coefficients.
    pub fn derivative(&self) -> Self {
        let n = self.coeffs.len();
        if n <= 1 {
            return Self {
                coeffs: alloc::vec![T::zero()],
            };
        }
        let mut out = alloc::vec![T::zero(); n - 1];
        out[n - 2] = self.coeffs[n - 1] * (2.0 * (n - 1) as f64);
        if n > 2 {
            for k in (1..=(n - 2)).rev() {
                let prev = if k + 1 < out.len() {
                    out[k + 1]
                } else {
                    T::zero()
                };
                out[k - 1] = prev + self.coeffs[k] * (2.0 * k as f64);
            }
        }
        out[0] = out[0] / 2.0;
        Self { coeffs: out }
    }

    /// Return a non-truncating normalized-coordinate indefinite integral
    /// series. Output length is `len() + 1`.
    pub fn integral(&self, constant: T) -> Self {
        let n = self.coeffs.len();
        let mut out = alloc::vec![T::zero(); n + 1];
        out[0] = constant;
        if n >= 1 {
            out[1] = out[1] + self.coeffs[0];
        }
        if n >= 2 {
            out[2] = out[2] + self.coeffs[1] / 4.0;
        }
        for k in 2..n {
            out[k + 1] = out[k + 1] + self.coeffs[k] / (2.0 * (k + 1) as f64);
            out[k - 1] = out[k - 1] - self.coeffs[k] / (2.0 * (k - 1) as f64);
        }
        Self { coeffs: out }
    }
}

#[cfg(feature = "alloc")]
impl ChebySeriesDyn<f64> {
    /// Return a copy with the constant (`T₀`) coefficient shifted by `delta`.
    ///
    /// Equivalent to evaluating `p(τ) + delta` without refitting.
    #[must_use]
    pub fn shifted_constant(&self, delta: f64) -> Self {
        let mut coeffs = self.coeffs.clone();
        if !coeffs.is_empty() && delta.is_finite() {
            coeffs[0] += delta;
        }
        Self { coeffs }
    }

    /// Shift the constant (`T₀`) coefficient in place by `delta`.
    pub fn shift_constant_in_place(&mut self, delta: f64) {
        if !self.coeffs.is_empty() && delta.is_finite() {
            self.coeffs[0] += delta;
        }
    }
}

/// A dynamic Chebyshev series tied to a physical domain.
#[cfg(feature = "alloc")]
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChebySeriesDynOn<T, X> {
    domain: Domain<X>,
    series: ChebySeriesDyn<T>,
}

#[cfg(feature = "alloc")]
impl<T: ChebyScalar, X: ChebyTime> ChebySeriesDynOn<T, X> {
    /// Create a domain-aware dynamic series.
    #[inline]
    pub fn new(domain: Domain<X>, series: ChebySeriesDyn<T>) -> Self {
        Self { domain, series }
    }

    /// Domain of validity.
    #[inline]
    pub fn domain(&self) -> Domain<X> {
        self.domain
    }

    /// Underlying normalized-coordinate series.
    #[inline]
    pub fn series(&self) -> &ChebySeriesDyn<T> {
        &self.series
    }

    /// Consume and return the underlying series.
    #[inline]
    pub fn into_series(self) -> ChebySeriesDyn<T> {
        self.series
    }

    /// Evaluate inside the domain.
    #[inline]
    pub fn evaluate(&self, x: X) -> Result<T, ChebyError> {
        if !self.domain.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        Ok(self.series.evaluate(self.domain.normalize(x)))
    }

    /// Evaluate the dimensionally-correct physical derivative at `x`.
    #[inline]
    pub fn evaluate_derivative(
        &self,
        x: X,
    ) -> Result<<T as DifferentiateWith<X>>::Derivative, ChebyError>
    where
        T: DifferentiateWith<X>,
    {
        if !self.domain.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        let d_tau = self.series.evaluate_both(self.domain.normalize(x)).1;
        Ok((d_tau * 2.0).scale_derivative(self.domain.end() - self.domain.start()))
    }
}
