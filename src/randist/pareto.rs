//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Pareto distribution of order a. The distribution function is,
/// 
/// p(x) dx = (a/b) / (x/b)^{a+1} dx
/// 
/// for x >= b.
pub fn pareto(r: &Rng, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_ran_pareto(ffi::FFI::unwrap(r), a, b) }
}

/// This function computes the probability density p(x) at x for a Pareto distribution with exponent a and scale b, using the formula given above.
pub fn pareto_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_ran_pareto_pdf(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
pub fn pareto_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_pareto_P(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
pub fn pareto_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_pareto_Q(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
pub fn pareto_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_pareto_Pinv(P, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
pub fn pareto_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { ffi::gsl_cdf_pareto_Qinv(Q, a, b) }
}