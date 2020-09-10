//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Laplace distribution with width a. The distribution is,
///
/// p(x) dx = {1 \over 2 a}  \exp(-|x/a|) dx
///
/// for -\infty < x < \infty.
pub fn laplace(r: &mut Rng, a: f64) -> f64 {
    unsafe { sys::gsl_ran_laplace(ffi::FFI::unwrap_unique(r), a) }
}

/// This function computes the probability density p(x) at x for a Laplace distribution with width a, using the formula given above.
pub fn laplace_pdf(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_ran_laplace_pdf(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Laplace distribution with width a.
pub fn laplace_P(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_laplace_P(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Laplace distribution with width a.
pub fn laplace_Q(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_laplace_Q(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Laplace distribution with width a.
pub fn laplace_Pinv(P: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_laplace_Pinv(P, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Laplace distribution with width a.
pub fn laplace_Qinv(Q: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_laplace_Qinv(Q, a) }
}
