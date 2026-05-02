//! Chebyshev differentiation matrices.

use crate::spectral::{collocation_points, Matrix};

/// Return the Chebyshev-Lobatto first differentiation matrix.
pub fn chebyshev_differentiation_matrix(n: usize) -> Matrix {
    if n == 0 {
        return Matrix::zeros(0, 0);
    }
    if n == 1 {
        return Matrix::zeros(1, 1);
    }

    let x = match n {
        2 => alloc::vec![1.0, -1.0],
        3 => alloc::vec![1.0, 0.0, -1.0],
        _ => (0..n)
            .map(|k| (core::f64::consts::PI * k as f64 / (n - 1) as f64).cos())
            .collect(),
    };
    let mut d = Matrix::zeros(n, n);
    let mut c = alloc::vec![1.0; n];
    c[0] = 2.0;
    c[n - 1] = 2.0;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
                d.set(i, j, sign * c[i] / c[j] / (x[i] - x[j]));
            }
        }
    }
    for i in 0..n {
        let row_sum: f64 = (0..n).filter(|&j| j != i).map(|j| d.get(i, j)).sum();
        d.set(i, i, -row_sum);
    }
    d
}

/// Fixed-size collocation helper for callers that need arrays.
#[inline]
pub fn fixed_collocation_points<const N: usize>() -> [f64; N] {
    collocation_points::<N>()
}
