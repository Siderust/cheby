// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Piecewise Chebyshev segment management.
//!
//! A [`ChebySegment`] stores Chebyshev coefficients and the domain of a
//! single interval, handling the `t → τ` normalisation internally.
//!
//! A [`ChebySegmentTable`] manages a sequence of uniform-duration segments
//! with automatic lookup, suitable for caching ephemeris-style data over
//! a time range.

use crate::eval;
use crate::fit;
use crate::scalar::{ChebyScalar, ChebyTime};

// ─────────────────────────────────────────────────────────────────────────
// ChebySegment — single segment
// ─────────────────────────────────────────────────────────────────────────

/// A single Chebyshev interpolation segment.
///
/// `T` is the value type (e.g., `f64`, `Kilometer`).
/// `Tt` is the time domain type (e.g., `f64`, `Quantity<Second>`).
/// `N` is the number of Chebyshev coefficients.
///
/// Stores `N` coefficients and the domain `[mid - half, mid + half]`.
/// The `eval*` methods handle the `t → τ` mapping automatically.
#[derive(Debug, Clone)]
pub struct ChebySegment<T: ChebyScalar, Tt: ChebyTime, const N: usize> {
    /// Chebyshev coefficients `c[0..N]`.
    pub coeffs: [T; N],
    /// Midpoint of the segment domain.
    pub mid: Tt,
    /// Half-width of the segment domain.
    pub half: Tt,
    /// Reciprocal of `half` as a raw `f64`, used for derivative chain-rule scaling.
    ///
    /// Stored so that `df/dt = (df/dτ) * half_inv` without needing to extract
    /// the value from `half` at evaluation time.
    pub half_inv: f64,
}

impl<T: ChebyScalar, Tt: ChebyTime, const N: usize> ChebySegment<T, Tt, N> {
    /// Create a segment from pre-computed coefficients and domain.
    ///
    /// `half_inv` is computed automatically from `half`.
    #[inline]
    pub fn new(coeffs: [T; N], mid: Tt, half: Tt) -> Self {
        let half_inv = half.recip_f64();
        Self {
            coeffs,
            mid,
            half,
            half_inv,
        }
    }

    /// Normalise `t` to `τ ∈ [-1, 1]` within this segment.
    ///
    /// For typed time (`Tt = Quantity<U>`), the division cancels the unit and
    /// returns a dimensionless `f64` — no `.value()` extraction needed.
    #[inline]
    pub fn normalise(&self, t: Tt) -> f64 {
        (t - self.mid) / self.half
    }

    /// Evaluate the Chebyshev polynomial at physical time `t`.
    #[inline]
    pub fn eval(&self, t: Tt) -> T {
        eval::evaluate(&self.coeffs, self.normalise(t))
    }

    /// Evaluate the derivative `df/dt` at physical time `t`.
    ///
    /// Accounts for the chain rule: `df/dt = (df/dτ) · (dτ/dt) = (df/dτ) · half_inv`.
    #[inline]
    pub fn eval_derivative(&self, t: Tt) -> T {
        eval::evaluate_derivative(&self.coeffs, self.normalise(t)) * self.half_inv
    }

    /// Evaluate both value and derivative `(f(t), df/dt)` in one pass.
    #[inline]
    pub fn eval_both(&self, t: Tt) -> (T, T) {
        let tau = self.normalise(t);
        let (v, d) = eval::evaluate_both(&self.coeffs, tau);
        (v, d * self.half_inv)
    }
}

/// A [`ChebySegment`] parameterised over `f64` time (backward-compatible alias).
pub type ChebySegmentF<T, const N: usize> = ChebySegment<T, f64, N>;

// ─────────────────────────────────────────────────────────────────────────
// ChebySegmentTable — uniform piecewise segments
// ─────────────────────────────────────────────────────────────────────────

/// A table of uniform-duration Chebyshev segments covering a time range.
///
/// `T` is the value type (e.g., `f64`, `Kilometer`).
/// `Tt` is the time domain type (e.g., `f64`, `Quantity<Second>`).
/// `N` is the number of Chebyshev coefficients per segment.
///
/// Each segment has the same duration; lookup is O(1) by index.
#[derive(Debug, Clone)]
pub struct ChebySegmentTable<T: ChebyScalar, Tt: ChebyTime, const N: usize> {
    /// Start of the first segment.
    start: Tt,
    /// Duration of each segment.
    segment_len: Tt,
    /// Segments, in chronological order.
    segments: Vec<ChebySegment<T, Tt, N>>,
}

impl<T: ChebyScalar, Tt: ChebyTime, const N: usize> ChebySegmentTable<T, Tt, N> {
    /// Build a segment table by sampling `f` at Chebyshev nodes within
    /// each segment.
    ///
    /// Divides `[start, end]` into segments of `segment_len` and fits
    /// Chebyshev coefficients in each.
    ///
    /// # Arguments
    ///
    /// - `f` — function to approximate; called at each Chebyshev node.
    /// - `start` — start of the domain.
    /// - `end` — end of the domain.
    /// - `segment_len` — duration of each segment.
    pub fn from_fn(f: impl Fn(Tt) -> T, start: Tt, end: Tt, segment_len: Tt) -> Self {
        let span_f64 = (end - start) / segment_len;
        let num_segments = (span_f64.ceil() as usize).max(1);
        let half = segment_len * 0.5;

        let mut segments = Vec::with_capacity(num_segments);
        for i in 0..num_segments {
            let seg_start = start + segment_len * (i as f64);
            let seg_end = seg_start + segment_len;
            let mid = seg_start + half;
            let coeffs = fit::fit_from_fn_t(&f, seg_start, seg_end);
            segments.push(ChebySegment::new(coeffs, mid, half));
        }

        Self {
            start,
            segment_len,
            segments,
        }
    }

    /// Build from pre-computed segments.
    pub fn from_segments(
        segments: Vec<ChebySegment<T, Tt, N>>,
        start: Tt,
        segment_len: Tt,
    ) -> Self {
        Self {
            start,
            segment_len,
            segments,
        }
    }

    /// Number of segments in the table.
    #[inline]
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Whether the table is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    /// Start of the covered domain.
    #[inline]
    pub fn start(&self) -> Tt {
        self.start
    }

    /// End of the covered domain.
    #[inline]
    pub fn end(&self) -> Tt {
        self.start + self.segment_len * (self.segments.len() as f64)
    }

    /// Duration of each segment.
    #[inline]
    pub fn segment_len(&self) -> Tt {
        self.segment_len
    }

    /// Look up the segment containing `t`, returning `None` if `t` is
    /// outside the table range.
    #[inline]
    pub fn get_segment(&self, t: Tt) -> Option<&ChebySegment<T, Tt, N>> {
        let offset = t - self.start;
        if offset < Tt::zero() {
            return None;
        }
        let idx = (offset / self.segment_len) as usize;
        self.segments.get(idx)
    }

    /// Evaluate at `t`, returning `None` if outside the table range.
    #[inline]
    pub fn eval(&self, t: Tt) -> Option<T> {
        self.get_segment(t).map(|s| s.eval(t))
    }

    /// Evaluate derivative at `t`, returning `None` if outside range.
    #[inline]
    pub fn eval_derivative(&self, t: Tt) -> Option<T> {
        self.get_segment(t).map(|s| s.eval_derivative(t))
    }

    /// Evaluate value and derivative at `t`, returning `None` if outside range.
    #[inline]
    pub fn eval_both(&self, t: Tt) -> Option<(T, T)> {
        self.get_segment(t).map(|s| s.eval_both(t))
    }

    /// Direct access to the underlying segments slice.
    #[inline]
    pub fn segments(&self) -> &[ChebySegment<T, Tt, N>] {
        &self.segments
    }
}

/// A [`ChebySegmentTable`] parameterised over `f64` time (backward-compatible alias).
pub type ChebySegmentTableF<T, const N: usize> = ChebySegmentTable<T, f64, N>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_eval_sin() {
        // Fit sin(t) on [0, π] as a single segment with 15 coefficients
        let coeffs: [f64; 15] = fit::fit_from_fn(f64::sin, 0.0, std::f64::consts::PI);
        let seg: ChebySegment<f64, f64, 15> = ChebySegment::new(
            coeffs,
            std::f64::consts::PI / 2.0,
            std::f64::consts::PI / 2.0,
        );

        // sin(π/2) = 1
        let val = seg.eval(std::f64::consts::FRAC_PI_2);
        assert!((val - 1.0).abs() < 1e-12, "sin(π/2) ≈ {val}");

        // sin(π/4) ≈ 0.7071
        let val2 = seg.eval(std::f64::consts::FRAC_PI_4);
        let exact = std::f64::consts::FRAC_PI_4.sin();
        assert!((val2 - exact).abs() < 1e-10, "sin(π/4) ≈ {val2}");
    }

    #[test]
    fn test_segment_derivative() {
        // d/dt sin(t) = cos(t)
        let coeffs: [f64; 15] = fit::fit_from_fn(f64::sin, 0.0, std::f64::consts::PI);
        let seg: ChebySegment<f64, f64, 15> = ChebySegment::new(
            coeffs,
            std::f64::consts::PI / 2.0,
            std::f64::consts::PI / 2.0,
        );

        let t = 1.0;
        let deriv = seg.eval_derivative(t);
        let exact = t.cos();
        assert!(
            (deriv - exact).abs() < 1e-10,
            "cos({t}) ≈ {deriv}, exact = {exact}"
        );
    }

    #[test]
    fn test_segment_eval_both() {
        let coeffs: [f64; 15] = fit::fit_from_fn(f64::sin, 0.0, std::f64::consts::PI);
        let seg: ChebySegment<f64, f64, 15> = ChebySegment::new(
            coeffs,
            std::f64::consts::PI / 2.0,
            std::f64::consts::PI / 2.0,
        );

        let t = 1.5;
        let (v, d) = seg.eval_both(t);
        assert!((v - seg.eval(t)).abs() < 1e-14);
        assert!((d - seg.eval_derivative(t)).abs() < 1e-14);
    }

    #[test]
    fn test_segment_normalise() {
        let seg: ChebySegment<f64, f64, 3> = ChebySegment::new([1.0_f64; 3], 5.0, 2.0);
        assert!((seg.normalise(5.0) - 0.0).abs() < 1e-15);
        assert!((seg.normalise(3.0) - (-1.0)).abs() < 1e-15);
        assert!((seg.normalise(7.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_table_from_fn() {
        // Approximate sin(t) on [0, 2π] with segments of length π/2
        let table: ChebySegmentTable<f64, f64, 9> = ChebySegmentTable::from_fn(
            f64::sin,
            0.0,
            2.0 * std::f64::consts::PI,
            std::f64::consts::FRAC_PI_2,
        );

        assert_eq!(table.len(), 4);

        // Check a few evaluation points
        for &t in &[0.1, 1.0, 2.0, 3.0, 5.0, 6.0] {
            let approx = table.eval(t).unwrap();
            let exact = t.sin();
            assert!(
                (approx - exact).abs() < 1e-8,
                "sin({t}): approx={approx}, exact={exact}"
            );
        }
    }

    #[test]
    fn test_table_derivative() {
        let table: ChebySegmentTable<f64, f64, 9> = ChebySegmentTable::from_fn(
            f64::sin,
            0.0,
            2.0 * std::f64::consts::PI,
            std::f64::consts::FRAC_PI_2,
        );

        let t = 2.0;
        let d = table.eval_derivative(t).unwrap();
        let exact = t.cos();
        assert!(
            (d - exact).abs() < 1e-8,
            "cos({t}): approx={d}, exact={exact}"
        );
    }

    #[test]
    fn test_table_metadata() {
        let table: ChebySegmentTable<f64, f64, 9> =
            ChebySegmentTable::from_fn(f64::sin, 1.0, 3.0, 0.5);
        assert_eq!(table.start(), 1.0);
        assert_eq!(table.segment_len(), 0.5);
        assert_eq!(table.end(), 3.0);
        assert_eq!(table.segments().len(), table.len());
    }

    #[test]
    fn test_table_out_of_range() {
        let table: ChebySegmentTable<f64, f64, 9> =
            ChebySegmentTable::from_fn(f64::sin, 0.0, 1.0, 0.5);
        assert!(table.eval(-0.1).is_none());
        // Just past the end
        assert!(table.eval(1.1).is_none());
    }

    #[test]
    fn test_segment_typed_time() {
        // Use qtty::Second as the time type, f64 as the value type.
        // Fit sin(t.value()) on [0, π] seconds.
        use qtty::Second;

        let start = Second::new(0.0);
        let end = Second::new(std::f64::consts::PI);
        let mid = Second::new(std::f64::consts::PI / 2.0);
        let half = Second::new(std::f64::consts::PI / 2.0);

        let coeffs: [f64; 15] = fit::fit_from_fn_t(&|t: Second| t.value().sin(), start, end);
        let seg: ChebySegment<f64, Second, 15> = ChebySegment::new(coeffs, mid, half);

        // sin(π/2) = 1
        let val = seg.eval(Second::new(std::f64::consts::FRAC_PI_2));
        assert!((val - 1.0).abs() < 1e-12, "typed time sin(π/2) ≈ {val}");

        // Derivative: d/dt sin(t) = cos(t)
        let t = Second::new(1.0);
        let deriv = seg.eval_derivative(t);
        let exact = 1.0_f64.cos();
        assert!(
            (deriv - exact).abs() < 1e-10,
            "typed time cos(1) ≈ {deriv}, exact = {exact}"
        );
    }

    #[test]
    fn test_table_typed_time() {
        // Use qtty::Second as the time type for a segment table.
        use qtty::Second;

        let table: ChebySegmentTable<f64, Second, 9> = ChebySegmentTable::from_fn(
            |t: Second| t.value().sin(),
            Second::new(0.0),
            Second::new(2.0 * std::f64::consts::PI),
            Second::new(std::f64::consts::FRAC_PI_2),
        );

        assert_eq!(table.len(), 4);

        for t_raw in &[0.1_f64, 1.0, 2.0, 3.0, 5.0, 6.0] {
            let approx = table.eval(Second::new(*t_raw)).unwrap();
            let exact = t_raw.sin();
            assert!(
                (approx - exact).abs() < 1e-8,
                "typed table sin({t_raw}): approx={approx}, exact={exact}"
            );
        }
    }
}
