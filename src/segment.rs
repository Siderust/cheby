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
use crate::scalar::ChebyScalar;

// ─────────────────────────────────────────────────────────────────────────
// ChebySegment — single segment
// ─────────────────────────────────────────────────────────────────────────

/// A single Chebyshev interpolation segment.
///
/// Stores `N` coefficients and the domain `[mid - half, mid + half]`.
/// The `eval*` methods handle the `t → τ` mapping automatically.
#[derive(Debug, Clone)]
pub struct ChebySegment<T: ChebyScalar, const N: usize> {
    /// Chebyshev coefficients `c[0..N]`.
    pub coeffs: [T; N],
    /// Midpoint of the segment domain.
    pub mid: f64,
    /// Half-width of the segment domain.
    pub half: f64,
}

impl<T: ChebyScalar, const N: usize> ChebySegment<T, N> {
    /// Create a segment from pre-computed coefficients and domain.
    #[inline]
    pub fn new(coeffs: [T; N], mid: f64, half: f64) -> Self {
        Self { coeffs, mid, half }
    }

    /// Normalise `t` to `τ ∈ [-1, 1]` within this segment.
    #[inline]
    pub fn normalise(&self, t: f64) -> f64 {
        (t - self.mid) / self.half
    }

    /// Evaluate the Chebyshev polynomial at physical time `t`.
    #[inline]
    pub fn eval(&self, t: f64) -> T {
        eval::evaluate(&self.coeffs, self.normalise(t))
    }

    /// Evaluate the derivative `df/dt` at physical time `t`.
    ///
    /// Accounts for the chain rule: `df/dt = (df/dτ) · (dτ/dt) = (df/dτ) / half`.
    #[inline]
    pub fn eval_derivative(&self, t: f64) -> T {
        let tau = self.normalise(t);
        eval::evaluate_derivative(&self.coeffs, tau) / self.half
    }

    /// Evaluate both value and derivative `(f(t), df/dt)` in one pass.
    #[inline]
    pub fn eval_both(&self, t: f64) -> (T, T) {
        let tau = self.normalise(t);
        let (v, d) = eval::evaluate_both(&self.coeffs, tau);
        (v, d / self.half)
    }
}

// ─────────────────────────────────────────────────────────────────────────
// ChebySegmentTable — uniform piecewise segments
// ─────────────────────────────────────────────────────────────────────────

/// A table of uniform-duration Chebyshev segments covering a time range.
///
/// Each segment has the same duration; lookup is O(1) by index.
#[derive(Debug, Clone)]
pub struct ChebySegmentTable<T: ChebyScalar, const N: usize> {
    /// Start of the first segment.
    start: f64,
    /// Duration of each segment.
    segment_len: f64,
    /// Segments, in chronological order.
    segments: Vec<ChebySegment<T, N>>,
}

impl<T: ChebyScalar, const N: usize> ChebySegmentTable<T, N> {
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
    pub fn from_fn(f: impl Fn(f64) -> T, start: f64, end: f64, segment_len: f64) -> Self {
        let span = end - start;
        let num_segments = ((span / segment_len).ceil() as usize).max(1);
        let half = segment_len * 0.5;

        let mut segments = Vec::with_capacity(num_segments);
        for i in 0..num_segments {
            let seg_start = start + i as f64 * segment_len;
            let seg_end = seg_start + segment_len;
            let mid = seg_start + half;
            let coeffs = fit::fit_from_fn(&f, seg_start, seg_end);
            segments.push(ChebySegment { coeffs, mid, half });
        }

        Self {
            start,
            segment_len,
            segments,
        }
    }

    /// Build from pre-computed segments.
    pub fn from_segments(segments: Vec<ChebySegment<T, N>>, start: f64, segment_len: f64) -> Self {
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
    pub fn start(&self) -> f64 {
        self.start
    }

    /// End of the covered domain.
    #[inline]
    pub fn end(&self) -> f64 {
        self.start + self.segments.len() as f64 * self.segment_len
    }

    /// Duration of each segment.
    #[inline]
    pub fn segment_len(&self) -> f64 {
        self.segment_len
    }

    /// Look up the segment containing `t`, returning `None` if `t` is
    /// outside the table range.
    #[inline]
    pub fn get_segment(&self, t: f64) -> Option<&ChebySegment<T, N>> {
        let offset = t - self.start;
        if offset < 0.0 {
            return None;
        }
        let idx = (offset / self.segment_len) as usize;
        self.segments.get(idx)
    }

    /// Evaluate at `t`, returning `None` if outside the table range.
    #[inline]
    pub fn eval(&self, t: f64) -> Option<T> {
        self.get_segment(t).map(|s| s.eval(t))
    }

    /// Evaluate derivative at `t`, returning `None` if outside range.
    #[inline]
    pub fn eval_derivative(&self, t: f64) -> Option<T> {
        self.get_segment(t).map(|s| s.eval_derivative(t))
    }

    /// Evaluate value and derivative at `t`, returning `None` if outside range.
    #[inline]
    pub fn eval_both(&self, t: f64) -> Option<(T, T)> {
        self.get_segment(t).map(|s| s.eval_both(t))
    }

    /// Direct access to the underlying segments slice.
    #[inline]
    pub fn segments(&self) -> &[ChebySegment<T, N>] {
        &self.segments
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_eval_sin() {
        // Fit sin(t) on [0, π] as a single segment with 15 coefficients
        let coeffs: [f64; 15] = fit::fit_from_fn(f64::sin, 0.0, std::f64::consts::PI);
        let seg = ChebySegment::new(
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
        let seg = ChebySegment::new(
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
        let seg = ChebySegment::new(
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
        let seg = ChebySegment::new([1.0_f64; 3], 5.0, 2.0);
        assert!((seg.normalise(5.0) - 0.0).abs() < 1e-15);
        assert!((seg.normalise(3.0) - (-1.0)).abs() < 1e-15);
        assert!((seg.normalise(7.0) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_table_from_fn() {
        // Approximate sin(t) on [0, 2π] with segments of length π/2
        let table: ChebySegmentTable<f64, 9> = ChebySegmentTable::from_fn(
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
        let table: ChebySegmentTable<f64, 9> = ChebySegmentTable::from_fn(
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
        let table: ChebySegmentTable<f64, 9> = ChebySegmentTable::from_fn(f64::sin, 1.0, 3.0, 0.5);
        assert_eq!(table.start(), 1.0);
        assert_eq!(table.segment_len(), 0.5);
        assert_eq!(table.end(), 3.0);
        assert_eq!(table.segments().len(), table.len());
    }

    #[test]
    fn test_table_out_of_range() {
        let table: ChebySegmentTable<f64, 9> = ChebySegmentTable::from_fn(f64::sin, 0.0, 1.0, 0.5);
        assert!(table.eval(-0.1).is_none());
        // Just past the end
        assert!(table.eval(1.1).is_none());
    }
}
