//! Root finding for dynamic Chebyshev series on `[-1, 1]`.

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

use super::ChebySeriesDyn;

const UNIT_ROOT_TOL: f64 = 1e-13;
const POLY_ZERO_TOL: f64 = 1e-12;
const DEFAULT_DEDUPE_EPS: f64 = 1e-10;

/// Options for [`ChebySeriesDyn::roots`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RootOptions {
    /// Tolerance for accepting a root on `[-1, 1]`.
    pub unit_tol: f64,
    /// Residual magnitude treated as zero when evaluating the polynomial.
    pub zero_tol: f64,
    /// Minimum separation between returned roots after deduplication.
    pub dedupe_eps: f64,
}

impl Default for RootOptions {
    fn default() -> Self {
        Self {
            unit_tol: UNIT_ROOT_TOL,
            zero_tol: POLY_ZERO_TOL,
            dedupe_eps: DEFAULT_DEDUPE_EPS,
        }
    }
}

#[cfg(feature = "alloc")]
impl ChebySeriesDyn<f64> {
    /// Sum of the absolute values of the last `tail` coefficients.
    ///
    /// When `tail` exceeds the series length, all coefficients are used.
    pub fn tail_norm(&self, tail: usize) -> f64 {
        let coeffs = self.coeffs();
        let tail = tail.max(1).min(coeffs.len());
        coeffs
            .iter()
            .rev()
            .take(tail)
            .map(|c| c.abs())
            .sum()
    }

    /// Find real roots in `[-1, 1]` using derivative segmentation and Brent
    /// refinement on each sign-change interval.
    pub fn roots(&self) -> Vec<f64> {
        self.roots_with(RootOptions::default())
    }

    /// Find real roots in `[-1, 1]` with explicit tolerances.
    pub fn roots_with(&self, opts: RootOptions) -> Vec<f64> {
        let mut roots = self.roots_recursive(0, opts);
        sort_dedup_f64(&mut roots, opts.dedupe_eps);
        roots
    }

    fn roots_recursive(&self, depth: usize, opts: RootOptions) -> Vec<f64> {
        let coeffs = self.coeffs();
        if coeffs.len() <= 1 {
            return Vec::new();
        }
        if coeffs.len() == 2 {
            let a = coeffs[1];
            if a.abs() <= opts.zero_tol {
                return Vec::new();
            }
            let x = -coeffs[0] / a;
            if (-1.0 - opts.unit_tol..=1.0 + opts.unit_tol).contains(&x) {
                return vec![x.clamp(-1.0, 1.0)];
            }
            return Vec::new();
        }

        let mut points = vec![-1.0, 1.0];
        if depth < coeffs.len() {
            points.extend(self.derivative().roots_recursive(depth + 1, opts));
        }
        sort_dedup_f64(&mut points, opts.dedupe_eps);

        let mut roots = Vec::new();
        for &x in &points {
            if self.evaluate(x).abs() <= opts.zero_tol {
                roots.push(x);
            }
        }

        for pair in points.windows(2) {
            let a = pair[0];
            let b = pair[1];
            if b - a <= opts.unit_tol {
                continue;
            }
            let fa = self.evaluate(a);
            let fb = self.evaluate(b);
            if fa.abs() <= opts.zero_tol || fb.abs() <= opts.zero_tol {
                continue;
            }
            if fa.signum() * fb.signum() < 0.0 {
                if let Some(root) =
                    brent_on_unit(a, b, fa, fb, |x| self.evaluate(x), opts.unit_tol)
                {
                    roots.push(root.clamp(-1.0, 1.0));
                }
            }
        }

        roots
    }
}

fn sort_dedup_f64(values: &mut Vec<f64>, eps: f64) {
    values.retain(|v| v.is_finite());
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    values.dedup_by(|a, b| (*a - *b).abs() <= eps);
}

fn brent_on_unit<F>(lo: f64, hi: f64, f_lo: f64, f_hi: f64, mut f: F, tol: f64) -> Option<f64>
where
    F: FnMut(f64) -> f64,
{
    if !lo.is_finite() || !hi.is_finite() || hi < lo {
        return None;
    }
    if f_lo.abs() <= POLY_ZERO_TOL {
        return Some(lo);
    }
    if f_hi.abs() <= POLY_ZERO_TOL {
        return Some(hi);
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

    for _ in 0..100 {
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

        let tol1 = 2.0 * f64::EPSILON * b.abs() + tol * 0.5;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb.abs() <= POLY_ZERO_TOL {
            return Some(b);
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
        fb = f(b);
        if !fb.is_finite() {
            return None;
        }
    }

    Some(b)
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
        assert!((roots[0] + 0.5).abs() < 1e-10, "{roots:?}");
        assert!((roots[1] - 0.5).abs() < 1e-10, "{roots:?}");
    }

    #[test]
    fn tail_norm_uses_last_coefficients() {
        let poly = fit_dyn_from_fn(8, |x| 2.0 * x * x - 1.0).unwrap();
        assert!(poly.tail_norm(4).is_finite());
        assert!(poly.tail_norm(4) >= 0.0);
    }
}
