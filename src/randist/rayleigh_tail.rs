//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the tail of the Rayleigh distribution with scale parameter sigma and a lower limit of a. The distribution is,
/// 
/// p(x) dx = {x \over \sigma^2} \exp ((a^2 - x^2) /(2 \sigma^2)) dx
/// 
/// for x > a.
pub fn rayleigh_tail(r: &Rng, a: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_rayleigh_tail(ffi::FFI::unwrap(r) as *const ffi::gsl_rng, a, sigma) }
}

/// This function computes the probability density p(x) at x for a Rayleigh tail distribution with scale parameter sigma and lower limit a, using the formula given above.
pub fn rayleigh_tail_pdf(x: f64, a: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_rayleigh_tail_pdf(x, a, sigma) }
}