//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random integer from the Poisson distribution with mean mu. The probability distribution for Poisson variates is,
/// 
/// p(k) = {\mu^k \over k!} \exp(-\mu)
/// 
/// for k >= 0.
pub fn poisson(r: &Rng, mu: f64) -> u32 {
    unsafe { ffi::gsl_ran_poisson(ffi::FFI::unwrap(r), mu) }
}

/// This function computes the probability p(k) of obtaining k from a Poisson distribution with mean mu, using the formula given above.
pub fn poisson_pdf(k: u32, mu: f64) -> f64 {
    unsafe { ffi::gsl_ran_poisson_pdf(k, mu) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the Poisson distribution with parameter mu.
pub fn poisson_P(k: u32, mu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_poisson_P(k, mu) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the Poisson distribution with parameter mu.
pub fn poisson_Q(k: u32, mu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_poisson_Q(k, mu) }
}