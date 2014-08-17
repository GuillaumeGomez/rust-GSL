//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Rayleigh distribution with scale parameter sigma. The distribution is,
/// 
/// p(x) dx = {x \over \sigma^2} \exp(- x^2/(2 \sigma^2)) dx
/// 
/// for x > 0.
pub fn rayleigh(r: &Rng, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_rayleigh(ffi::FFI::unwrap(r) as *const ffi::gsl_rng, sigma) }
}

/// This function computes the probability density p(x) at x for a Rayleigh distribution with scale parameter sigma, using the formula given above.
pub fn rayleigh_pdf(x: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_rayleigh_pdf(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
pub fn rayleigh_P(x: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_rayleigh_P(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
pub fn rayleigh_Q(x: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_rayleigh_Q(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
pub fn rayleigh_Pinv(P: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_rayleigh_Pinv(P, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
pub fn rayleigh_Qinv(Q: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_rayleigh_Qinv(Q, sigma) }
}