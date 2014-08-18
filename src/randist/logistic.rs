//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the logistic distribution. The distribution function is,
/// 
/// p(x) dx = { \exp(-x/a) \over a (1 + \exp(-x/a))^2 } dx
/// 
/// for -\infty < x < +\infty.
pub fn logistic(r: &Rng, a: f64) -> f64 {
    unsafe { ffi::gsl_ran_logistic(ffi::FFI::unwrap(r) as *const ffi::gsl_rng, a) }
}

/// This function computes the probability density p(x) at x for a logistic distribution with scale parameter a, using the formula given above.
pub fn logistic_pdf(x: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_ran_logistic_pdf(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
pub fn logistic_P(x: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_logistic_P(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
pub fn logistic_Q(x: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_logistic_Q(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
pub fn logistic_Pinv(P: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_logistic_Pinv(P, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
pub fn logistic_Qinv(Q: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_logistic_Qinv(Q, a) }
}