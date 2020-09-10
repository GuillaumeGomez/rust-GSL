//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the lognormal distribution. The distribution function is,
///
/// p(x) dx = {1 \over x \sqrt{2 \pi \sigma^2} } \exp(-(\ln(x) - \zeta)^2/2 \sigma^2) dx
///
/// for x > 0.
pub fn lognormal(r: &mut Rng, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_ran_lognormal(ffi::FFI::unwrap_unique(r), zeta, sigma) }
}

/// This function computes the probability density p(x) at x for a lognormal distribution with parameters zeta and sigma, using the formula given above.
pub fn lognormal_pdf(x: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_ran_lognormal_pdf(x, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
pub fn lognormal_P(x: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_P(x, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
pub fn lognormal_Q(x: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_Q(x, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
pub fn lognormal_Pinv(P: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_Pinv(P, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
pub fn lognormal_Qinv(Q: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_Qinv(Q, zeta, sigma) }
}
