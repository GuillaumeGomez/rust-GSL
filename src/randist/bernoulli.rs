//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns either 0 or 1, the result of a Bernoulli trial with probability p. The probability distribution for a Bernoulli trial is,
///
/// p(0) = 1 - p
/// p(1) = p
pub fn bernoulli(r: &mut Rng, p: f64) -> u32 {
    unsafe { ffi::randist::gsl_ran_bernoulli(ffi::FFI::unwrap_unique(r), p) }
}

/// This function computes the probability p(k) of obtaining k from a Bernoulli distribution with probability parameter p, using the formula given above.
pub fn bernoulli_pdf(x: u32, p: f64) -> f64 {
    unsafe { ffi::randist::gsl_ran_bernoulli_pdf(x, p) }
}
