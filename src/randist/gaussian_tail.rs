//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function provides random variates from the upper tail of a Gaussian distribution with standard deviation sigma.
/// The values returned are larger than the lower limit a, which must be positive. The method is based on Marsaglia’s famous rectangle-wedge-tail algorithm (Ann. Math. Stat. 32, 894–899 (1961)), with this aspect explained in Knuth, v2, 3rd ed, p139,586 (exercise 11).
/// 
/// The probability distribution for Gaussian tail random variates is,
/// 
/// p(x) dx = {1 \over N(a;\sigma) \sqrt{2 \pi \sigma^2}} \exp (- x^2/(2 \sigma^2)) dx
/// 
/// for x > a where N(a;\sigma) is the normalization constant,
/// 
/// N(a;\sigma) = (1/2) erfc(a / sqrt(2 sigma^2)).
pub fn gaussian_tail(r: &Rng, a: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_gaussian_tail(ffi::FFI::unwrap(r), a, sigma) }
}

/// This function computes the probability density p(x) at x for a Gaussian tail distribution with standard deviation sigma and lower limit a, using the formula given above.
pub fn gaussian_tail_pdf(x: f64, a: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_gaussian_tail_pdf(x, a, sigma) }
}

/// This function computes results for the tail of a unit Gaussian distribution. They are equivalent to the functions above with a standard deviation of one, sigma = 1.
pub fn ugaussian_tail(r: &Rng, a: f64) -> f64 {
    unsafe { ffi::gsl_ran_ugaussian_tail(ffi::FFI::unwrap(r), a) }
}

/// This function computes results for the tail of a unit Gaussian distribution. They are equivalent to the functions above with a standard deviation of one, sigma = 1.
pub fn ugaussian_tail_pdf(x: f64, a: f64) -> f64 {
    unsafe { ffi::gsl_ran_ugaussian_tail_pdf(x, a) }
}