//! Encode/decode a Chebyshev series via the binary I/O subsystem.
//!
//! Run with `cargo run --example binary_roundtrip --features binary`.

use cheby::core::ChebySeriesDyn;
use cheby::io::binary::{decode_f64_series, encode_f64_series};

fn main() {
    let series = ChebySeriesDyn::<f64>::new(vec![1.0, -0.25, 0.0625, -0.03125, 0.015625])
        .expect("non-empty series");

    let bytes = encode_f64_series(&series);
    println!("encoded {} bytes", bytes.len());

    let decoded = decode_f64_series(&bytes).expect("round-trip");
    assert_eq!(series.coeffs(), decoded.coeffs());
    println!("round-trip OK on {} coefficients", series.coeffs().len());

    // Demonstrate a typed error variant: tampering flips the checksum.
    let mut tampered = bytes.clone();
    *tampered.last_mut().unwrap() ^= 0xFF;
    match decode_f64_series(&tampered) {
        Err(e) => println!("tampered payload rejected: {e:?}"),
        Ok(_) => unreachable!("checksum should have failed"),
    }
}
