//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Landau distribution. The probability distribution for Landau random variates is defined analytically by the complex integral,
///
/// p(x) = (1/(2 \pi i)) \int_{c-i\infty}^{c+i\infty} ds exp(s log(s) + x s)
///
/// For numerical purposes it is more convenient to use the following equivalent form of the integral,
///
/// p(x) = (1/\pi) \int_0^\infty dt \exp(-t \log(t) - x t) \sin(\pi t).
pub fn landau(r: &mut Rng) -> f64 {
    unsafe { ffi::gsl_ran_landau(ffi::FFI::unwrap_unique(r)) }
}

/// This function computes the probability density p(x) at x for the Landau distribution using an approximation to the formula given above.
pub fn landau_pdf(x: f64) -> f64 {
    unsafe { ffi::gsl_ran_landau_pdf(x) }
}
