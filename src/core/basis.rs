//! Chebyshev basis functions.

/// Evaluate `T_n(x)`, the Chebyshev polynomial of the first kind.
#[inline]
pub fn t(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut t0 = 1.0;
            let mut t1 = x;
            for _ in 2..=n {
                let next = 2.0 * x * t1 - t0;
                t0 = t1;
                t1 = next;
            }
            t1
        }
    }
}

/// Evaluate `U_n(x)`, the Chebyshev polynomial of the second kind.
#[inline]
pub fn u(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 2.0 * x,
        _ => {
            let mut u0 = 1.0;
            let mut u1 = 2.0 * x;
            for _ in 2..=n {
                let next = 2.0 * x * u1 - u0;
                u0 = u1;
                u1 = next;
            }
            u1
        }
    }
}
