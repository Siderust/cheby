//! Compact binary coefficient-table helpers for `f64` series.

use alloc::vec::Vec;

use crate::core::{ChebyError, ChebySeriesDyn};
use crate::io::metadata::FORMAT_VERSION;

/// Encode an `f64` dynamic series.
pub fn encode_f64_series(series: &ChebySeriesDyn<f64>) -> Vec<u8> {
    let mut out = Vec::with_capacity(16 + 8 * series.coeffs().len());
    out.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
    out.extend_from_slice(&(series.coeffs().len() as u64).to_le_bytes());
    for coeff in series.coeffs() {
        out.extend_from_slice(&coeff.to_le_bytes());
    }
    let checksum = checksum(&out);
    out.extend_from_slice(&checksum.to_le_bytes());
    out
}

/// Decode an `f64` dynamic series.
pub fn decode_f64_series(bytes: &[u8]) -> Result<ChebySeriesDyn<f64>, ChebyError> {
    if bytes.len() < 20 {
        return Err(ChebyError::EmptyCoefficientSet);
    }
    let payload_len = bytes.len() - 8;
    let expected = u64::from_le_bytes(bytes[payload_len..].try_into().unwrap());
    if checksum(&bytes[..payload_len]) != expected {
        return Err(ChebyError::InvalidDomain);
    }
    let version = u32::from_le_bytes(bytes[0..4].try_into().unwrap());
    if version != FORMAT_VERSION {
        return Err(ChebyError::InvalidDegree);
    }
    let len = u64::from_le_bytes(bytes[4..12].try_into().unwrap()) as usize;
    if payload_len != 12 + 8 * len {
        return Err(ChebyError::InvalidDegree);
    }
    let mut coeffs = Vec::with_capacity(len);
    for chunk in bytes[12..payload_len].chunks_exact(8) {
        coeffs.push(f64::from_le_bytes(chunk.try_into().unwrap()));
    }
    ChebySeriesDyn::new(coeffs)
}

fn checksum(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf29ce484222325, |hash, b| {
        (hash ^ u64::from(*b)).wrapping_mul(0x100000001b3)
    })
}
