//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the flat (uniform) distribution from a to b. The distribution is,
/// 
/// p(x) dx = {1 \over (b-a)} dx
/// 
/// if a <= x < b and 0 otherwise.
pub fn flat(r: &Rng, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_ran_flat(ffi::FFI::unwrap(r), a, b) }
}

/// This function computes the probability density p(x) at x for a uniform distribution from a to b, using the formula given above.
pub fn flat_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_ran_flat_pdf(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for a uniform distribution from a to b.
pub fn flat_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_flat_P(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for a uniform distribution from a to b.
pub fn flat_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_flat_Q(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for a uniform distribution from a to b.
pub fn flat_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_flat_Pinv(P, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for a uniform distribution from a to b.
pub fn flat_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_flat_Qinv(Q, a, b) }
}