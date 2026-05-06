//! Compact binary coefficient-table helpers for `f64` series.
//!
//! # Format (little-endian)
//!
//! ```text
//! offset  size  field
//! ------  ----  ----------------------------------
//!  0       4    u32 format_version  (== FORMAT_VERSION)
//!  4       8    u64 coefficient_count
//! 12       8*N  N * f64 coefficients
//! 12+8N    8    u64 FNV-1a checksum over [0, 12+8N)
//! ```
//!
//! [`decode_f64_series`] reports each malformation with its own typed
//! [`ChebyError`] variant, so callers can distinguish a truncated blob from
//! a version mismatch from a checksum failure.

use alloc::vec::Vec;

use crate::core::{ChebyError, ChebySeriesDyn};
use crate::io::metadata::FORMAT_VERSION;

const HEADER_LEN: usize = 12;
const CHECKSUM_LEN: usize = 8;
const MIN_LEN: usize = HEADER_LEN + CHECKSUM_LEN;

/// Encode an `f64` dynamic series using the [module format](self).
pub fn encode_f64_series(series: &ChebySeriesDyn<f64>) -> Vec<u8> {
    let mut out = Vec::with_capacity(MIN_LEN + 8 * series.coeffs().len());
    out.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
    out.extend_from_slice(&(series.coeffs().len() as u64).to_le_bytes());
    for coeff in series.coeffs() {
        out.extend_from_slice(&coeff.to_le_bytes());
    }
    let checksum = checksum(&out);
    out.extend_from_slice(&checksum.to_le_bytes());
    out
}

/// Decode an `f64` dynamic series produced by [`encode_f64_series`].
pub fn decode_f64_series(bytes: &[u8]) -> Result<ChebySeriesDyn<f64>, ChebyError> {
    if bytes.len() < MIN_LEN {
        return Err(ChebyError::BinaryTooShort);
    }
    let payload_len = bytes.len() - CHECKSUM_LEN;
    let expected = u64::from_le_bytes(bytes[payload_len..].try_into().unwrap());
    if checksum(&bytes[..payload_len]) != expected {
        return Err(ChebyError::BinaryChecksumMismatch);
    }
    let version = u32::from_le_bytes(bytes[0..4].try_into().unwrap());
    if version != FORMAT_VERSION {
        return Err(ChebyError::UnsupportedFormatVersion {
            found: version,
            expected: FORMAT_VERSION,
        });
    }
    let len = u64::from_le_bytes(bytes[4..HEADER_LEN].try_into().unwrap()) as usize;
    if payload_len != HEADER_LEN + 8 * len {
        return Err(ChebyError::BinaryLengthMismatch);
    }
    let mut coeffs = Vec::with_capacity(len);
    for chunk in bytes[HEADER_LEN..payload_len].chunks_exact(8) {
        coeffs.push(f64::from_le_bytes(chunk.try_into().unwrap()));
    }
    ChebySeriesDyn::new(coeffs)
}

fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf29ce484222325, |hash, b| {
        (hash ^ u64::from(*b)).wrapping_mul(0x100000001b3)
    })
}
