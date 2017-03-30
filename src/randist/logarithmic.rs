//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random integer from the logarithmic distribution. The probability distribution for logarithmic random variates is,
/// 
/// p(k) = {-1 \over \log(1-p)} {(p^k \over k)}
/// 
/// for k >= 1.
pub fn logarithmic(r: &mut Rng, p: f64) -> u32 {
    unsafe { ffi::gsl_ran_logarithmic(ffi::FFI::unwrap_unique(r), p) }
}

/// This function computes the probability p(k) of obtaining k from a logarithmic distribution with probability parameter p, using the formula given above.
pub fn logarithmic_pdf(k: u32, p: f64) -> f64 {
    unsafe { ffi::gsl_ran_logarithmic_pdf(k, p) }
}