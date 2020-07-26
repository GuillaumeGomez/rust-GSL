//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The chi-squared distribution arises in statistics. If Y_i are n independent Gaussian random variates with unit variance then the sum-of-squares,

X_i = \sum_i Y_i^2

has a chi-squared distribution with n degrees of freedom.
!*/

use ffi;
use types::Rng;

/// This function returns a random variate from the chi-squared distribution with nu degrees of freedom. The distribution function is,
///
/// p(x) dx = {1 \over 2 Gamma(\nu/2) } (x/2)^{\nu/2 - 1} \exp(-x/2) dx
///
/// for x >= 0.
pub fn chisq(r: &mut Rng, nu: f64) -> f64 {
    unsafe { ffi::gsl_ran_chisq(ffi::FFI::unwrap_unique(r), nu) }
}

/// This function computes the probability density p(x) at x for a chi-squared distribution with nu degrees of freedom, using the formula given above.
pub fn chisq_pdf(x: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_ran_chisq_pdf(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
pub fn chisq_P(x: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_chisq_P(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
pub fn chisq_Q(x: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_chisq_Q(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
pub fn chisq_Pinv(P: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_chisq_Pinv(P, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
pub fn chisq_Qinv(Q: f64, nu: f64) -> f64 {
    unsafe { ffi::gsl_cdf_chisq_Qinv(Q, nu) }
}
