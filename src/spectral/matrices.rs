//! Lightweight matrix storage.

use alloc::vec::Vec;

/// Dense row-major matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// Create a zero matrix.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: alloc::vec![0.0; rows * cols],
        }
    }

    /// Row count.
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Column count.
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Row-major backing data.
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    /// Get an entry.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.cols + col]
    }

    /// Set an entry.
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        self.data[row * self.cols + col] = value;
    }
}
