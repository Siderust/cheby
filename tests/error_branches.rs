//! Coverage for every public [`ChebyError`] branch.

use cheby::core::{ChebyError, ChebySeriesDyn, Domain};

#[test]
fn domain_error_branches() {
    assert_eq!(
        Domain::<f64>::try_new(1.0, 1.0).unwrap_err(),
        ChebyError::EmptyDomain
    );
    assert_eq!(
        Domain::<f64>::try_new(1.0, 0.0).unwrap_err(),
        ChebyError::InvalidDomain
    );
    assert_eq!(
        Domain::<f64>::try_new(f64::NAN, 1.0).unwrap_err(),
        ChebyError::NonFiniteInput
    );
    assert_eq!(
        Domain::<f64>::try_new(0.0, f64::INFINITY).unwrap_err(),
        ChebyError::NonFiniteInput
    );
}

#[test]
fn dyn_series_empty_set() {
    assert_eq!(
        ChebySeriesDyn::<f64>::new(vec![]).unwrap_err(),
        ChebyError::EmptyCoefficientSet
    );
}

#[test]
fn series_on_evaluate_out_of_domain() {
    let domain = Domain::new(0.0, 1.0);
    let series = cheby::ChebySeriesOn::new(domain, cheby::ChebySeries::new([1.0, 0.0]));
    assert_eq!(
        series.evaluate(2.0).unwrap_err(),
        ChebyError::EvaluationOutOfDomain
    );
}

#[test]
#[cfg(feature = "binary")]
fn binary_error_branches() {
    use cheby::core::ChebySeriesDyn;
    use cheby::io::binary::{decode_f64_series, encode_f64_series};

    // Round-trip succeeds.
    let original = ChebySeriesDyn::new(vec![1.0, -2.5, 0.125, 7.0]).unwrap();
    let bytes = encode_f64_series(&original);
    let decoded = decode_f64_series(&bytes).unwrap();
    assert_eq!(decoded.coeffs(), original.coeffs());

    // Truncated.
    assert_eq!(
        decode_f64_series(&bytes[..5]).unwrap_err(),
        ChebyError::BinaryTooShort
    );

    // Corrupt checksum.
    let mut corrupted = bytes.clone();
    let last = corrupted.len() - 1;
    corrupted[last] ^= 0xFF;
    assert_eq!(
        decode_f64_series(&corrupted).unwrap_err(),
        ChebyError::BinaryChecksumMismatch
    );

    // Wrong version (recompute checksum so we hit version branch, not checksum).
    let mut v = bytes.clone();
    // Use a supported layout with an unsupported version byte.
    v[0] = 99;
    let payload_len = v.len() - 8;
    let new_ck = v[..payload_len]
        .iter()
        .fold(0xcbf29ce484222325_u64, |h, b| {
            (h ^ u64::from(*b)).wrapping_mul(0x100000001b3)
        });
    v[payload_len..].copy_from_slice(&new_ck.to_le_bytes());
    let err = decode_f64_series(&v).unwrap_err();
    assert!(matches!(
        err,
        ChebyError::UnsupportedFormatVersion {
            found: _,
            expected: _
        }
    ));

    // Length mismatch: tamper coefficient count, recompute checksum.
    let mut l = bytes;
    l[4..12].copy_from_slice(&999_u64.to_le_bytes());
    let payload_len = l.len() - 8;
    let new_ck = l[..payload_len]
        .iter()
        .fold(0xcbf29ce484222325_u64, |h, b| {
            (h ^ u64::from(*b)).wrapping_mul(0x100000001b3)
        });
    l[payload_len..].copy_from_slice(&new_ck.to_le_bytes());
    assert_eq!(
        decode_f64_series(&l).unwrap_err(),
        ChebyError::BinaryLengthMismatch
    );
}
