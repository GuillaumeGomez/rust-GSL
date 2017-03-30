//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random integer from the geometric distribution, the number of independent trials with probability p until the first success.
/// The probability distribution for geometric variates is,
/// 
/// p(k) =  p (1-p)^(k-1)
/// 
/// for k >= 1. Note that the distribution begins with k=1 with this definition. There is another convention in which the exponent k-1 is replaced by k.
pub fn geometric(r: &mut Rng, p: f64) -> u32 {
    unsafe { ffi::gsl_ran_geometric(ffi::FFI::unwrap_unique(r), p) }
}

/// This function computes the probability p(k) of obtaining k from a geometric distribution with probability parameter p, using the formula given above.
pub fn geometric_pdf(k: u32, p: f64) -> f64 {
    unsafe { ffi::gsl_ran_geometric_pdf(k, p) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the geometric distribution with parameter p.
pub fn geometric_P(k: u32, p: f64) -> f64 {
    unsafe { ffi::gsl_cdf_geometric_P(k, p) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the geometric distribution with parameter p.
pub fn geometric_Q(k: u32, p: f64) -> f64 {
    unsafe { ffi::gsl_cdf_geometric_Q(k, p) }
}