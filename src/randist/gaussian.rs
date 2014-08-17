//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a Gaussian random variate, with mean zero and standard deviation sigma.
/// The probability distribution for Gaussian random variates is,
/// 
/// p(x) dx = {1 \over \sqrt{2 \pi \sigma^2}} \exp (-x^2 / 2\sigma^2) dx
/// for x in the range -\infty to +\infty. Use the transformation z = \mu + x on the numbers returned by gsl_ran_gaussian to obtain a Gaussian distribution with mean \mu.
/// This function uses the Box-Muller algorithm which requires two calls to the random number generator r.
pub fn gaussian(r: &Rng, sigma: f64) -> f64 {
	unsafe { ffi::gsl_ran_gaussian(ffi::FFI::unwrap(r) as *const ffi::gsl_rng, sigma) }
}