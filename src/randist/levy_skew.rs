//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random variate from the Levy skew stable distribution with scale c, exponent alpha and skewness parameter beta.
/// The skewness parameter must lie in the range [-1,1]. The Levy skew stable probability distribution is defined by a Fourier transform,
///
/// p(x) = {1 \over 2 \pi} \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^alpha (1-i beta sign(t) tan(pi alpha/2)))
///
/// When \alpha = 1 the term \tan(\pi \alpha/2) is replaced by -(2/\pi)\log|t|. There is no explicit solution for the form of p(x) and the library does not define a corresponding pdf function.
/// For \alpha = 2 the distribution reduces to a Gaussian distribution with \sigma = \sqrt{2} c and the skewness parameter has no effect. For \alpha < 1 the tails of the distribution become extremely wide.
/// The symmetric distribution corresponds to \beta = 0.
///
/// The algorithm only works for 0 < alpha <= 2.
///
/// The Levy alpha-stable distributions have the property that if N alpha-stable variates are drawn from the distribution p(c, \alpha, \beta) then the sum Y = X_1 + X_2 + \dots + X_N will also be distributed as an alpha-stable variate, p(N^(1/\alpha) c, \alpha, \beta).
pub fn levy_skew(r: &mut Rng, c: f64, alpha: f64, beta: f64) -> f64 {
    unsafe { sys::gsl_ran_levy_skew(ffi::FFI::unwrap_unique(r), c, alpha, beta) }
}
