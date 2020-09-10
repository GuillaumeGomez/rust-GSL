//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the gamma distribution. The distribution function is,
///
/// p(x) dx = {1 over Gamma(a) b^a} x^{a-1} e^{-x/b} dx
///
/// for x > 0.
///
/// The gamma distribution with an integer parameter a is known as the Erlang distribution.
///
/// The variates are computed using the Marsaglia-Tsang fast gamma method. This function for this method was previously called gsl_ran_gamma_mt and can still be accessed using this name.
pub fn gamma(r: &mut Rng, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_gamma(ffi::FFI::unwrap_unique(r), a, b) }
}

/// This function returns a gamma variate using the algorithms from Knuth (vol 2).
pub fn gamma_knuth(r: &mut Rng, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_gamma_knuth(ffi::FFI::unwrap_unique(r), a, b) }
}

/// This function computes the probability density p(x) at x for a gamma distribution with parameters a and b, using the formula given above.
pub fn gamma_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_gamma_pdf(x, a, b) }
}

/// This function computes This function computes the probability density p(x) at x for a gamma distribution with parameters a and b, using the formula given above.
pub fn gamma_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_gamma_P(x, a, b) }
}

/// This function computes This function computes the probability density p(x) at x for a gamma distribution with parameters a and b, using the formula given above.
pub fn gamma_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_gamma_Q(x, a, b) }
}

/// This function computes This function computes the probability density p(x) at x for a gamma distribution with parameters a and b, using the formula given above.
pub fn gamma_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_gamma_Pinv(P, a, b) }
}

/// This function computes This function computes the probability density p(x) at x for a gamma distribution with parameters a and b, using the formula given above.
pub fn gamma_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_gamma_Qinv(Q, a, b) }
}
