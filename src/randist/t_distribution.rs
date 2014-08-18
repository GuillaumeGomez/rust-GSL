//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The t-distribution arises in statistics. If Y_1 has a normal distribution and Y_2 has a chi-squared distribution with \nu degrees of freedom then the ratio,

X = { Y_1 \over \sqrt{Y_2 / \nu} }

has a t-distribution t(x;\nu) with \nu degrees of freedom.
!*/

use ffi;
use types::Rng;

/// This function returns a random variate from the t-distribution. The distribution function is,
/// 
/// p(x) dx = {Gamma((\nu + 1)/2) \over \sqrt{\pi \nu} Gamma(\nu/2)}
/// 
///    (1 + x^2/\nu)^{-(\nu + 1)/2} dx
/// 
/// for -\infty < x < +\infty.
pub fn tdist(r: &Rng, nu: f64) -> f64 {
    unsafe { ffi::gsl_ran_tdist(ffi::FFI::unwrap(r) as *const ffi::gsl_rng, nu) }
}

/// This function computes the probability density p(x) at x for a t-distribution with nu degrees of freedom, using the formula given above.
pub fn tdist_pdf(x: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_ran_tdist_pdf(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the t-distribution with nu degrees of freedom.
pub fn tdist_P(x: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_tdist_P(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the t-distribution with nu degrees of freedom.
pub fn tdist_Q(x: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_tdist_Q(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the t-distribution with nu degrees of freedom.
pub fn tdist_Pinv(P: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_tdist_Pinv(P, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the t-distribution with nu degrees of freedom.
pub fn tdist_Qinv(Q: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_tdist_Qinv(Q, nu) }
}