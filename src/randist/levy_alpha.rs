//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Levy symmetric stable distribution with scale c and exponent alpha. The symmetric stable probability distribution is defined by a Fourier transform,
///
/// p(x) = {1 \over 2 \pi} \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^alpha)
///
/// There is no explicit solution for the form of p(x) and the library does not define a corresponding pdf function. For \alpha = 1 the distribution reduces to the Cauchy distribution. For \alpha = 2 it is a Gaussian distribution with \sigma = \sqrt{2} c. For \alpha < 1 the tails of the distribution become extremely wide.
///
/// The algorithm only works for 0 < alpha <= 2.
pub fn levy(r: &mut Rng, c: f64, alpha: f64) -> f64 {
    unsafe { ffi::randist::gsl_ran_levy(ffi::FFI::unwrap_unique(r), c, alpha) }
}
