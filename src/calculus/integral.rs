//! Integral helpers.

use crate::core::{ChebyError, ChebyScalar, ChebySeries, ChebySeriesOn, ChebyTime, Domain};

pub use crate::core::IntegrateWith;

/// A precomputed antiderivative shifted so that its value at the domain
/// start is zero.
///
/// Build one with [`ChebySeriesOn::definite_integral_from_start`] when you
/// will repeatedly query the definite integral `∫_{start}^{x} f(t) dt` for
/// many `x`. Each evaluation is then a single Clenshaw pass over the
/// `(N + 1)`-length antiderivative; no per-call coefficient construction
/// or normalization is performed.
#[derive(Debug, Clone, PartialEq)]
pub struct DefiniteIntegralFromStart<T, X, const N: usize>
where
    T: ChebyScalar,
    X: ChebyTime,
{
    domain: Domain<X>,
    /// Antiderivative coefficients, length `N + 1`, with `b_0` chosen so
    /// `A(-1) = 0` (i.e. evaluation at the domain start vanishes).
    coeffs: [T; N],
    coeffs_tail: T,
}

impl<T, X, const N: usize> DefiniteIntegralFromStart<T, X, N>
where
    T: ChebyScalar + IntegrateWith<X>,
    X: ChebyTime,
{
    /// Domain of validity.
    #[inline]
    pub fn domain(&self) -> Domain<X> {
        self.domain
    }

    /// Evaluate `∫_{start}^{x} f(t) dt`.
    pub fn evaluate(&self, x: X) -> Result<<T as IntegrateWith<X>>::Integral, ChebyError> {
        if !self.domain.contains(x) {
            return Err(ChebyError::EvaluationOutOfDomain);
        }
        let tau = self.domain.normalize(x);
        // Clenshaw with the `N + 1`-coefficient antiderivative `coeffs ++ [coeffs_tail]`.
        let normalized = clenshaw_with_tail::<T, N>(&self.coeffs, self.coeffs_tail, tau);
        Ok(normalized.scale_integral(self.domain.half_width()))
    }
}

#[inline]
fn clenshaw_with_tail<T: ChebyScalar, const N: usize>(coeffs: &[T; N], tail: T, x: f64) -> T {
    // Evaluate the `(N+1)`-coefficient series `[coeffs[0..N], tail]` via Clenshaw,
    // without allocating.
    let two_x = 2.0 * x;
    let mut b1 = T::zero();
    let mut b2 = T::zero();
    // Process the tail first (highest index = N).
    let mut b0 = b1 * two_x - b2 + tail;
    b2 = b1;
    b1 = b0;
    for k in (1..N).rev() {
        b0 = b1 * two_x - b2 + coeffs[k];
        b2 = b1;
        b1 = b0;
    }
    if N == 0 {
        // Series is just `[tail]`; b1 currently holds `tail`, b2 = 0.
        return b1;
    }
    coeffs[0] + b1 * x - b2
}

impl<T, X, const N: usize> ChebySeriesOn<T, X, N>
where
    T: ChebyScalar + IntegrateWith<X>,
    X: ChebyTime,
{
    /// Evaluate the definite integral from the domain start to `x`.
    ///
    /// Convenience for one-shot integration. For repeated queries prefer
    /// [`Self::definite_integral_from_start`], which precomputes the shifted
    /// antiderivative once.
    pub fn evaluate_integral_from_start(
        &self,
        x: X,
    ) -> Result<<T as IntegrateWith<X>>::Integral, ChebyError> {
        self.definite_integral_from_start().evaluate(x)
    }

    /// Build a reusable precomputed antiderivative for repeated definite
    /// integral queries from the domain start.
    pub fn definite_integral_from_start(&self) -> DefiniteIntegralFromStart<T, X, N> {
        let domain = self.domain();
        let coeffs_in = self.series().coeffs();
        // Build the full non-truncated antiderivative directly into a
        // length-(N+1) representation: store the first `N` entries inline
        // and the highest entry separately so we can stay const-generic
        // without `generic_const_exprs`.
        let mut out = [T::zero(); N];
        let mut tail = T::zero();
        if N >= 1 {
            // ∫T_0 = T_1: contributes c_0 to b_1.
            if N >= 2 {
                out[1] = out[1] + coeffs_in[0];
            } else {
                // Only one input coefficient; b_1 lives in the tail slot.
                tail = tail + coeffs_in[0];
            }
        }
        if N >= 2 {
            // ∫T_1 = T_2/4: contributes c_1/4 to b_2.
            if N >= 3 {
                out[2] = out[2] + coeffs_in[1] / 4.0;
            } else {
                tail = tail + coeffs_in[1] / 4.0;
            }
        }
        for k in 2..N {
            let plus = coeffs_in[k] / (2.0 * (k + 1) as f64);
            let minus = coeffs_in[k] / (2.0 * (k - 1) as f64);
            // b_{k+1} += plus  (k+1 may equal N => goes in tail)
            if k + 1 < N {
                out[k + 1] = out[k + 1] + plus;
            } else {
                tail = tail + plus;
            }
            // b_{k-1} -= minus  (k-1 < N always)
            out[k - 1] = out[k - 1] - minus;
        }
        // Choose b_0 so A(-1) = 0. With T_k(-1) = (-1)^k:
        //   A(-1) = sum_{k=0}^{N} b_k (-1)^k.
        let mut shift = T::zero();
        let mut sign = 1.0_f64;
        for &b in out.iter() {
            shift = shift + b * sign;
            sign = -sign;
        }
        // Tail contribution at index N: sign is (-1)^N.
        let tail_sign = if N.is_multiple_of(2) { 1.0 } else { -1.0 };
        shift = shift + tail * tail_sign;
        if N >= 1 {
            out[0] = out[0] - shift;
        } else {
            tail = tail - shift;
        }
        DefiniteIntegralFromStart {
            domain,
            coeffs: out,
            coeffs_tail: tail,
        }
    }
}

// Suppress an unused-import warning when the impl above is the only consumer
// of `ChebySeries`.
#[allow(dead_code)]
fn _force_chebyseries_use<T: ChebyScalar, const N: usize>(_: &ChebySeries<T, N>) {}
