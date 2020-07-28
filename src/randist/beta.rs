//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the beta distribution. The distribution function is,
///
/// p(x) dx = {Gamma(a+b) over Gamma(a) Gamma(b)} x^{a-1} (1-x)^{b-1} dx
///
/// for 0 <= x <= 1.
pub fn beta(r: &mut Rng, a: f64, b: f64) -> f64 {
    unsafe { ffi::randist::gsl_ran_beta(ffi::FFI::unwrap_unique(r), a, b) }
}

/// This function computes the probability density p(x) at x for a beta distribution with parameters a and b, using the formula given above.
pub fn beta_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::randist::gsl_ran_beta_pdf(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_beta_P(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_beta_Q(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_beta_Pinv(P, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_beta_Qinv(Q, a, b) }
}
