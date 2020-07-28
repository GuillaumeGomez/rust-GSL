//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random integer from the negative binomial distribution, the number of failures occurring before n successes in independent trials with
/// probability p of success. The probability distribution for negative binomial variates is,
///
/// p(k) = {\Gamma(n + k) \over \Gamma(k+1) \Gamma(n) } p^n (1-p)^k
///
/// Note that n is not required to be an integer.
pub fn negative_binomial(r: &mut Rng, p: f64, n: f64) -> u32 {
    unsafe { ffi::randist::gsl_ran_negative_binomial(ffi::FFI::unwrap_unique(r), p, n) }
}

/// This function computes the probability p(k) of obtaining k from a negative binomial distribution with parameters p and n, using the formula given above.
pub fn negative_binomial_pdf(k: u32, p: f64, n: f64) -> f64 {
    unsafe { ffi::randist::gsl_ran_negative_binomial_pdf(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the negative binomial distribution with parameters p and n.
pub fn negative_binomial_P(k: u32, p: f64, n: f64) -> f64 {
    unsafe { ffi::gsl_cdf_negative_binomial_P(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the negative binomial distribution with parameters p and n.
pub fn negative_binomial_Q(k: u32, p: f64, n: f64) -> f64 {
    unsafe { ffi::gsl_cdf_negative_binomial_Q(k, p, n) }
}
