//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Weibull distribution. The distribution function is,
/// 
/// p(x) dx = {b \over a^b} x^{b-1}  \exp(-(x/a)^b) dx
/// 
/// for x >= 0.
pub fn weidbull(r: &Rng, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_ran_weidbull(ffi::FFI::unwrap(r) as *const ffi::gsl_rng, a, b) }
}

/// This function computes the probability density p(x) at x for a Weibull distribution with scale a and exponent b, using the formula given above.
pub fn weidbull_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_ran_weidbull_pdf(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weidbull_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_weidbull_P(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weidbull_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_weidbull_Q(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weidbull_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_weidbull_Pinv(P, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weidbull_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_weidbull_Qinv(Q, a, b) }
}