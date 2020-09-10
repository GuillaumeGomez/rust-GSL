//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random integer from the Pascal distribution. The Pascal distribution is simply a negative binomial distribution with an integer value of n.
///
/// p(k) = {(n + k - 1)! \over k! (n - 1)! } p^n (1-p)^k
///
/// for k >= 0
pub fn pascal(r: &mut Rng, p: f64, n: u32) -> u32 {
    unsafe { sys::gsl_ran_pascal(ffi::FFI::unwrap_unique(r), p, n) }
}

/// This function computes the probability p(k) of obtaining k from a Pascal distribution with parameters p and n, using the formula given above.
pub fn pascal_pdf(k: u32, p: f64, n: u32) -> f64 {
    unsafe { sys::gsl_ran_pascal_pdf(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the Pascal distribution with parameters p and n.
pub fn pascal_P(k: u32, p: f64, n: u32) -> f64 {
    unsafe { sys::gsl_cdf_pascal_P(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the Pascal distribution with parameters p and n.
pub fn pascal_Q(k: u32, p: f64, n: u32) -> f64 {
    unsafe { sys::gsl_cdf_pascal_Q(k, p, n) }
}
