//! Segment lookup helpers.

use crate::core::ChebyTime;

/// Compute a uniform segment index for `x`.
#[inline]
pub fn uniform_index<X: ChebyTime>(start: X, segment_len: X, len: usize, x: X) -> Option<usize> {
    if len == 0 || x < start {
        return None;
    }
    let idx = ((x - start) / segment_len).floor() as usize;
    (idx < len).then_some(idx)
}
