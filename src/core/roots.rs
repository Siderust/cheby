//! Root finding for dynamic Chebyshev series on the normalized interval `[-1, 1]`.
//!
//! # Algorithm
//!
//! Real roots of a Chebyshev expansion `p(τ) = Σ aₖ Tₖ(τ)` on `τ ∈ [-1, 1]` are
//! found in two stages:
//!
//! 1. **Primary pass:** trim trailing near-zero high coefficients, reject non-finite
//!    input, convert the Chebyshev expansion to a monic power-basis polynomial, and
//!    locate roots with a **Durand–Kerner** iteration. Candidates are kept when they
//!    are real (`|Im| ≤ unit_tol`), lie in `[-1, 1]`, and satisfy
//!    `|p(τ)| ≤ zero_tol` after Clenshaw evaluation.
//! 2. **Fallback pass:** if too few roots are found, scan Chebyshev nodes on `[-1, 1]`,
//!    refine sign-change brackets with **Brent's method**, and search for tangency
//!    minima with Newton and golden-section refinement. Fallback thresholds scale
//!    with [`RootOptions::zero_tol`].
//!
//! Linear series (`len == 2`) use a closed-form root. Constant series return an
//! empty list (see below). Repeated or tangent roots are only resolved to the extent
//! allowed by [`RootOptions::zero_tol`] and [`RootOptions::dedupe_eps`].
//!
//! # Constant and degenerate series
//!
//! A series with a single coefficient (constant), an identically-zero series, or a
//! near-zero constant within `zero_tol` has **no isolated roots** on `[-1, 1]` and
//! returns an empty list. A non-zero constant likewise returns an empty list.
//!
//! # Coordinate system
//!
//! [`ChebySeriesDyn::roots`] and [`ChebySeriesDyn::roots_with`] return roots in the
//! **normalized** coordinate `τ ∈ [-1, 1]`. Map to a physical domain with
//! [`Domain::denormalize`](crate::Domain::denormalize) or use [`ChebySeriesDynOn::roots`].
//!
//! # Limitations
//!
//! - Intended for **moderate degrees** (roughly up to a few dozen coefficients).
//!   Cost grows superlinearly with degree because all roots are found at once.
//! - Root finding converts to power basis internally; this can be **less stable than
//!   Clenshaw evaluation** at high degree even when residuals pass `zero_tol`.
//! - Roots are **numerical approximations** controlled by [`RootOptions`].
//! - Multiple or tangent roots are harder; tighten `zero_tol` or reduce degree when
//!   residuals are unsatisfactory.
//! - Non-finite coefficients yield an empty root list.

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

use super::{ChebyError, ChebySeriesDyn, ChebySeriesDynOn, ChebyTime};

const DEFAULT_UNIT_TOL: f64 = 1e-13;
const DEFAULT_ZERO_TOL: f64 = 1e-12;
const DEFAULT_DEDUPE_EPS: f64 = 1e-10;
const MIN_POSITIVE_TOL: f64 = f64::EPSILON;
const DURAND_KERNER_MAX_ITER: usize = 80;

/// Options for [`ChebySeriesDyn::roots_with`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RootOptions {
    /// Bracket width and root-polishing scale on `[-1, 1]`.
    pub unit_tol: f64,
    /// Residual magnitude `|p(τ)|` treated as a root.
    pub zero_tol: f64,
    /// Minimum separation between returned roots after deduplication.
    pub dedupe_eps: f64,
}

impl Default for RootOptions {
    fn default() -> Self {
        Self {
            unit_tol: DEFAULT_UNIT_TOL,
            zero_tol: DEFAULT_ZERO_TOL,
            dedupe_eps: DEFAULT_DEDUPE_EPS,
        }
    }
}

impl RootOptions {
    /// Validate tolerances, returning an error for non-finite, zero, or negative values.
    pub fn validated(self) -> Result<Self, ChebyError> {
        Ok(Self {
            unit_tol: positive_finite_tol(self.unit_tol)?,
            zero_tol: positive_finite_tol(self.zero_tol)?,
            dedupe_eps: positive_finite_tol(self.dedupe_eps)?,
        })
    }

    /// Return sanitized tolerances, substituting defaults for invalid fields.
    #[inline]
    pub fn effective(self) -> Self {
        self.validated().unwrap_or_default()
    }
}

fn positive_finite_tol(value: f64) -> Result<f64, ChebyError> {
    if !value.is_finite() || value <= 0.0 {
        Err(ChebyError::NonFiniteInput)
    } else {
        Ok(value.max(MIN_POSITIVE_TOL))
    }
}

fn coefficients_are_finite(coeffs: &[f64]) -> bool {
    coeffs.iter().all(|c| c.is_finite())
}

fn on_unit_interval(x: f64, unit_tol: f64) -> bool {
    x.is_finite() && (-1.0 - unit_tol..=1.0 + unit_tol).contains(&x)
}

fn clamp_unit(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

fn trim_trailing_coeffs(coeffs: &[f64], zero_tol: f64) -> &[f64] {
    let mut end = coeffs.len();
    while end > 1 && coeffs[end - 1].abs() <= zero_tol {
        end -= 1;
    }
    &coeffs[..end]
}

#[cfg(feature = "alloc")]
impl ChebySeriesDyn<f64> {
    /// Sum of the absolute values of the last `tail` coefficients.
    ///
    /// When `tail` exceeds the series length, all coefficients are used.
    pub fn tail_norm(&self, tail: usize) -> f64 {
        let coeffs = self.coeffs();
        let tail = tail.max(1).min(coeffs.len());
        coeffs.iter().rev().take(tail).map(|c| c.abs()).sum()
    }

    /// Find real roots in the normalized interval `[-1, 1]`.
    ///
    /// Returns roots as normalized coordinates `τ`, not physical domain values.
    /// See [`ChebySeriesDynOn::roots`] for domain-mapped roots.
    pub fn roots(&self) -> Vec<f64> {
        self.roots_with(RootOptions::default())
    }

    /// Find real roots in `[-1, 1]` with explicit tolerances.
    ///
    /// Invalid [`RootOptions`] fields are replaced with defaults via
    /// [`RootOptions::effective`]. Non-finite coefficients yield an empty list.
    ///
    /// Constant, identically-zero, and near-zero constant series return an empty
    /// list because they have no isolated roots on `[-1, 1]`.
    pub fn roots_with(&self, opts: RootOptions) -> Vec<f64> {
        let opts = opts.effective();
        if !coefficients_are_finite(self.coeffs()) {
            return Vec::new();
        }

        let coeffs = trim_trailing_coeffs(self.coeffs(), opts.zero_tol);
        if coeffs.len() <= 1 {
            return Vec::new();
        }
        if coeffs.len() == 2 {
            return self.linear_zero(coeffs, opts);
        }

        let power = chebyshev_to_power(coeffs);
        let mut roots = durand_kerner_real_roots(&power, opts);
        let expected = coeffs.len() - 1;
        if roots.len() < expected {
            supplement_roots_on_chebyshev_nodes(self, &mut roots, opts);
        }
        sort_dedup_f64(&mut roots, opts.dedupe_eps);
        roots.retain(|&r| self.evaluate(r).abs() <= opts.zero_tol);
        if roots.is_empty() {
            supplement_roots_on_chebyshev_nodes(self, &mut roots, opts);
            sort_dedup_f64(&mut roots, opts.dedupe_eps);
            roots.retain(|&r| self.evaluate(r).abs() <= opts.zero_tol);
        }
        roots
    }

    fn linear_zero(&self, coeffs: &[f64], opts: RootOptions) -> Vec<f64> {
        debug_assert_eq!(coeffs.len(), 2);
        let a = coeffs[1];
        if a.abs() <= opts.zero_tol {
            return Vec::new();
        }
        let x = -coeffs[0] / a;
        if on_unit_interval(x, opts.unit_tol) {
            let root = clamp_unit(x);
            if self.evaluate(root).abs() <= opts.zero_tol {
                return vec![root];
            }
        }
        Vec::new()
    }
}

#[cfg(feature = "alloc")]
impl<X: ChebyTime> ChebySeriesDynOn<f64, X> {
    /// Find real roots mapped into the series domain.
    pub fn roots(&self) -> Vec<X> {
        self.roots_with(RootOptions::default())
    }

    /// Find real roots mapped into the series domain with explicit tolerances.
    pub fn roots_with(&self, opts: RootOptions) -> Vec<X> {
        self.series()
            .roots_with(opts)
            .into_iter()
            .map(|tau| self.domain().denormalize(tau))
            .collect()
    }
}

fn chebyshev_to_power(coeffs: &[f64]) -> Vec<f64> {
    let degree = coeffs.len().saturating_sub(1);
    let mut power = vec![0.0; degree + 1];
    if degree == 0 {
        power[0] = coeffs[0];
        return power;
    }

    let mut t0 = vec![0.0; degree + 1];
    let mut t1 = vec![0.0; degree + 1];
    t0[0] = 1.0;
    t1[1] = 1.0;
    accumulate_scaled(&mut power, &t0, coeffs[0]);
    if degree >= 1 {
        accumulate_scaled(&mut power, &t1, coeffs[1]);
    }

    for (_, &ck) in coeffs.iter().enumerate().take(degree + 1).skip(2) {
        let mut tk = vec![0.0; degree + 1];
        for i in 0..=degree {
            if i >= 1 {
                tk[i] += 2.0 * t1[i - 1];
            }
            tk[i] -= t0[i];
        }
        accumulate_scaled(&mut power, &tk, ck);
        t0 = t1;
        t1 = tk;
    }

    power
}

fn accumulate_scaled(dst: &mut [f64], src: &[f64], scale: f64) {
    if scale == 0.0 {
        return;
    }
    for (d, &s) in dst.iter_mut().zip(src.iter()) {
        *d += scale * s;
    }
}

fn durand_kerner_real_roots(power: &[f64], opts: RootOptions) -> Vec<f64> {
    let degree = power.len().saturating_sub(1);
    if degree == 0 {
        return Vec::new();
    }

    let scale = power[degree];
    if !scale.is_finite() || scale.abs() <= opts.zero_tol {
        return Vec::new();
    }

    let mut monic = power.to_vec();
    for c in &mut monic {
        *c /= scale;
    }

    let mut roots = initial_durand_kerner_guess(degree);
    for _ in 0..DURAND_KERNER_MAX_ITER {
        let mut next = roots.clone();
        for i in 0..degree {
            let value = eval_power(&monic, roots[i].re, roots[i].im);
            let mut denom = Complex { re: 1.0, im: 0.0 };
            for (j, zj) in roots.iter().enumerate().take(degree) {
                if i != j {
                    denom = cmul(denom, csub(roots[i], *zj));
                }
            }
            if denom.abs() <= f64::EPSILON {
                continue;
            }
            next[i] = csub(roots[i], cdiv(value, denom));
        }
        roots = next;
    }

    let mut out = Vec::new();
    for z in roots {
        if z.im.abs() > opts.unit_tol {
            continue;
        }
        let x = clamp_unit(z.re);
        if on_unit_interval(x, opts.unit_tol) {
            out.push(x);
        }
    }
    out
}

fn initial_durand_kerner_guess(degree: usize) -> Vec<Complex> {
    let mut roots = Vec::with_capacity(degree);
    for k in 0..degree {
        let angle = core::f64::consts::TAU * (k as f64) / (degree as f64);
        roots.push(Complex {
            re: 0.8 * angle.cos(),
            im: 0.8 * angle.sin(),
        });
    }
    roots
}

#[derive(Clone, Copy, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn abs(self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
}

fn cadd(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re + b.re,
        im: a.im + b.im,
    }
}

fn csub(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re - b.re,
        im: a.im - b.im,
    }
}

fn cmul(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re * b.re - a.im * b.im,
        im: a.re * b.im + a.im * b.re,
    }
}

fn cdiv(a: Complex, b: Complex) -> Complex {
    let denom = b.re * b.re + b.im * b.im;
    Complex {
        re: (a.re * b.re + a.im * b.im) / denom,
        im: (a.im * b.re - a.re * b.im) / denom,
    }
}

fn eval_power(coeffs: &[f64], re: f64, im: f64) -> Complex {
    let z = Complex { re, im };
    let mut acc = Complex {
        re: coeffs[coeffs.len() - 1],
        im: 0.0,
    };
    for &c in coeffs[..coeffs.len() - 1].iter().rev() {
        acc = cadd(cmul(acc, z), Complex { re: c, im: 0.0 });
    }
    acc
}

fn supplement_roots_on_chebyshev_nodes(
    series: &ChebySeriesDyn<f64>,
    roots: &mut Vec<f64>,
    opts: RootOptions,
) {
    let degree = series.coeffs().len().saturating_sub(1).max(1);
    let nodes = degree.saturating_mul(8).clamp(64, 1024);
    let mut prev_x = -1.0;
    let mut prev_f = series.evaluate(prev_x);

    for k in 0..nodes {
        let x = if k + 1 == nodes {
            1.0
        } else {
            let theta =
                core::f64::consts::PI * (2.0 * (nodes - 1 - k) as f64 + 1.0) / (2.0 * nodes as f64);
            theta.cos()
        };
        let fx = series.evaluate(x);
        if fx.abs() <= opts.zero_tol {
            if !has_nearby_root(roots, x, opts.dedupe_eps) {
                roots.push(clamp_unit(x));
            }
        } else if fx.abs() <= opts.zero_tol * 10.0 {
            if let Some(root) = polish_newton_minimum(series, x, opts) {
                if !has_nearby_root(roots, root, opts.dedupe_eps) {
                    roots.push(root);
                }
            }
        } else if prev_f.abs() > opts.zero_tol && prev_f.signum() * fx.signum() < 0.0 {
            if let Some(root) = brent_on_unit(
                prev_x,
                x,
                prev_f,
                fx,
                |t| series.evaluate(t),
                opts.unit_tol,
                opts.zero_tol,
            ) {
                if !has_nearby_root(roots, root, opts.dedupe_eps) {
                    roots.push(clamp_unit(root));
                }
            }
        }
        prev_x = x;
        prev_f = fx;
    }

    let mut best_x = 0.0;
    let mut best_abs = f64::INFINITY;
    for k in 0..nodes {
        let x = if k + 1 == nodes {
            1.0
        } else {
            let theta =
                core::f64::consts::PI * (2.0 * (nodes - 1 - k) as f64 + 1.0) / (2.0 * nodes as f64);
            theta.cos()
        };
        let fx = series.evaluate(x).abs();
        if fx < best_abs {
            best_abs = fx;
            best_x = x;
        }
    }
    if best_abs <= opts.zero_tol * 100.0 {
        let delta = (2.0 / nodes as f64).max(opts.unit_tol * 4.0);
        let lo = clamp_unit(best_x - delta);
        let hi = clamp_unit(best_x + delta);
        if let Some(root) = refine_minimum_magnitude(series, lo, hi, opts) {
            if !has_nearby_root(roots, root, opts.dedupe_eps) {
                roots.push(root);
            }
        } else if best_abs <= opts.zero_tol && !has_nearby_root(roots, best_x, opts.dedupe_eps) {
            roots.push(clamp_unit(best_x));
        }
    }
}

fn refine_minimum_magnitude(
    series: &ChebySeriesDyn<f64>,
    mut lo: f64,
    mut hi: f64,
    opts: RootOptions,
) -> Option<f64> {
    for _ in 0..48 {
        if hi - lo <= opts.unit_tol {
            break;
        }
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if series.evaluate(m1).abs() <= series.evaluate(m2).abs() {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let x = clamp_unit(0.5 * (lo + hi));
    if series.evaluate(x).abs() <= opts.zero_tol {
        Some(x)
    } else {
        None
    }
}

fn has_nearby_root(roots: &[f64], x: f64, eps: f64) -> bool {
    roots.iter().any(|r| (*r - x).abs() <= eps)
}

fn polish_newton_minimum(
    series: &ChebySeriesDyn<f64>,
    start: f64,
    opts: RootOptions,
) -> Option<f64> {
    let mut x = clamp_unit(start);
    for _ in 0..32 {
        let (value, deriv) = series.evaluate_both(x);
        if value.abs() <= opts.zero_tol {
            return Some(clamp_unit(x));
        }
        if deriv.abs() <= f64::EPSILON {
            break;
        }
        x = clamp_unit(x - value / deriv);
    }
    let residual = series.evaluate(x).abs();
    if residual <= opts.zero_tol {
        Some(clamp_unit(x))
    } else {
        None
    }
}

fn sort_dedup_f64(values: &mut Vec<f64>, eps: f64) {
    values.retain(|v| v.is_finite());
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    values.dedup_by(|a, b| (*a - *b).abs() <= eps);
}

fn brent_on_unit<F>(
    lo: f64,
    hi: f64,
    f_lo: f64,
    f_hi: f64,
    mut f: F,
    unit_tol: f64,
    zero_tol: f64,
) -> Option<f64>
where
    F: FnMut(f64) -> f64,
{
    if !lo.is_finite() || !hi.is_finite() || hi < lo {
        return None;
    }
    if f_lo.abs() <= zero_tol {
        return verified_root(lo, &mut f, unit_tol, zero_tol);
    }
    if f_hi.abs() <= zero_tol {
        return verified_root(hi, &mut f, unit_tol, zero_tol);
    }
    if f_lo.signum() * f_hi.signum() > 0.0 {
        return None;
    }

    let mut a = lo;
    let mut b = hi;
    let mut fa = f_lo;
    let mut fb = f_hi;
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;
    const BRENT_MAX_ITER: usize = 100;

    for _ in 0..BRENT_MAX_ITER {
        if fb.signum() * fc.signum() > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        let tol1 = 2.0 * f64::EPSILON * b.abs() + unit_tol * 0.5;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 {
            return verified_root(b, &mut f, unit_tol, zero_tol);
        }
        if fb.abs() <= zero_tol {
            return verified_root(b, &mut f, unit_tol, zero_tol);
        }

        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q) = if (a - c).abs() <= f64::EPSILON {
                (2.0 * xm * s, 1.0 - s)
            } else {
                let q = fa / fc;
                let r = fb / fc;
                (
                    s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0)),
                    (q - 1.0) * (r - 1.0) * (s - 1.0),
                )
            };
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();

            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if 2.0 * p < min1.min(min2) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }

        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else {
            b += tol1.copysign(xm);
        }
        b = clamp_unit(b);
        fb = f(b);
        if !fb.is_finite() {
            return None;
        }
    }

    verified_root(b, &mut f, unit_tol, zero_tol)
}

fn verified_root<F>(x: f64, f: &mut F, unit_tol: f64, zero_tol: f64) -> Option<f64>
where
    F: FnMut(f64) -> f64,
{
    if !on_unit_interval(x, unit_tol) {
        return None;
    }
    let value = f(x);
    if value.is_finite() && value.abs() <= zero_tol {
        Some(clamp_unit(x))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::approx::fit::fit_dyn_from_fn;

    #[test]
    fn quadratic_roots_on_unit_interval() {
        let poly = fit_dyn_from_fn(12, |x| x * x - 0.25).unwrap();
        let roots = poly.roots();
        assert_eq!(roots.len(), 2, "{roots:?}");
        assert_roots_within_tol(&poly, &roots, RootOptions::default());
        assert!((roots[0] + 0.5).abs() < 1e-10, "{roots:?}");
        assert!((roots[1] - 0.5).abs() < 1e-10, "{roots:?}");
    }

    #[test]
    fn tail_norm_uses_last_coefficients() {
        let poly = fit_dyn_from_fn(8, |x| 2.0 * x * x - 1.0).unwrap();
        assert!(poly.tail_norm(4).is_finite());
        assert!(poly.tail_norm(4) >= 0.0);
    }

    #[test]
    fn shift_constant_moves_evaluation() {
        let poly = fit_dyn_from_fn(6, |x| x).unwrap();
        let shifted = poly.shifted_constant(2.0);
        assert!((shifted.evaluate(0.0) - poly.evaluate(0.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn triple_root_fitted() {
        let p = fit_dyn_from_fn(16, |x| (x + 0.2).powi(3)).unwrap();
        let opts = RootOptions {
            zero_tol: 1e-5,
            unit_tol: 1e-11,
            dedupe_eps: 1e-10,
        };
        let roots = p.roots_with(opts);
        assert!(!roots.is_empty(), "{roots:?}");
        assert!(roots.iter().any(|r| (*r + 0.2).abs() < 5e-3), "{roots:?}");
    }

    #[test]
    fn constant_zero_series_has_no_roots() {
        let p = ChebySeriesDyn::new(vec![0.0]).unwrap();
        assert!(p.roots_with(RootOptions::default()).is_empty());
    }

    #[test]
    fn chebyshev_polynomial_has_exact_degree_roots() {
        for n in [4usize, 8, 12] {
            let mut coeffs = vec![0.0_f64; n + 1];
            coeffs[n] = 1.0;
            let p = ChebySeriesDyn::new(coeffs).unwrap();
            let opts = RootOptions {
                zero_tol: 1e-9,
                unit_tol: 1e-11,
                dedupe_eps: 1e-10,
            };
            let roots = p.roots_with(opts);
            assert_eq!(roots.len(), n, "T_{n} should have {n} roots, got {roots:?}");
            assert_roots_within_tol(&p, &roots, opts);
        }
    }

    pub(crate) fn assert_roots_within_tol(
        poly: &ChebySeriesDyn<f64>,
        roots: &[f64],
        opts: RootOptions,
    ) {
        let opts = opts.effective();
        for &root in roots {
            assert!(
                poly.evaluate(root).abs() <= opts.zero_tol,
                "residual at {root} = {}",
                poly.evaluate(root)
            );
            assert!(on_unit_interval(root, opts.unit_tol));
        }
    }
}
