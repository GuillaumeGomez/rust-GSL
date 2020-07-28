//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random integer from the binomial distribution, the number of successes in n independent trials with probability p. The probability distribution for binomial variates is,
///
/// p(k) = {n! \over k! (n-k)! } p^k (1-p)^{n-k}
///
/// for 0 <= k <= n.
pub fn binomial(r: &mut Rng, p: f64, n: u32) -> u32 {
    unsafe { ffi::randist::gsl_ran_binomial(ffi::FFI::unwrap_unique(r), p, n) }
}

/// This function computes the probability p(k) of obtaining k from a binomial distribution with parameters p and n, using the formula given above.
pub fn binomial_pdf(k: u32, p: f64, n: u32) -> f64 {
    unsafe { ffi::randist::gsl_ran_binomial_pdf(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the binomial distribution with parameters p and n.
pub fn binomial_P(k: u32, p: f64, n: u32) -> f64 {
    unsafe { ffi::gsl_cdf_binomial_P(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the binomial distribution with parameters p and n.
pub fn binomial_Q(k: u32, p: f64, n: u32) -> f64 {
    unsafe { ffi::gsl_cdf_binomial_Q(k, p, n) }
}
