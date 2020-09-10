//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the exponential power distribution with scale parameter a and exponent b. The distribution is,
///
/// p(x) dx = {1 \over 2 a Gamma(1+1/b)} \exp(-|x/a|^b) dx
///
/// for x >= 0. For b = 1 this reduces to the Laplace distribution. For b = 2 it has the same form as a Gaussian distribution, but with a = \sqrt{2} \sigma.
pub fn exppow(r: &mut Rng, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_exppow(ffi::FFI::unwrap_unique(r), a, b) }
}

/// This function computes the probability density p(x) at x for an exponential power distribution with scale parameter a and exponent b, using the formula given above.
pub fn exppow_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_exppow_pdf(x, a, b) }
}

/// This function computes tthe cumulative distribution functions P(x), Q(x) for the exponential power distribution with parameters a and b.
pub fn exppow_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_exppow_P(x, a, b) }
}

/// This function computes tthe cumulative distribution functions P(x), Q(x) for the exponential power distribution with parameters a and b.
pub fn exppow_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_exppow_Q(x, a, b) }
}
