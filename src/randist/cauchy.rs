//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Cauchy distribution with scale parameter a. The probability distribution for Cauchy random variates is,
///
/// p(x) dx = {1 \over a\pi (1 + (x/a)^2) } dx
///
/// for x in the range -\infty to +\infty. The Cauchy distribution is also known as the Lorentz distribution.
pub fn cauchy(r: &mut Rng, a: f64) -> f64 {
    unsafe { ffi::randist::gsl_ran_cauchy(ffi::FFI::unwrap_unique(r), a) }
}

/// This function computes the probability density p(x) at x for a Cauchy distribution with scale parameter a, using the formula given above.
pub fn cauchy_pdf(x: f64, a: f64) -> f64 {
    unsafe { ffi::randist::gsl_ran_cauchy_pdf(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
pub fn cauchy_P(x: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_cauchy_P(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
pub fn cauchy_Q(x: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_cauchy_Q(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
pub fn cauchy_Pinv(P: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_cauchy_Pinv(P, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
pub fn cauchy_Qinv(Q: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_cdf_cauchy_Qinv(Q, a) }
}
