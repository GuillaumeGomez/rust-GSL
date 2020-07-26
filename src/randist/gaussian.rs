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
pub fn gaussian(r: &mut Rng, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_gaussian(ffi::FFI::unwrap_unique(r), sigma) }
}

/// This function computes the probability density p(x) at x for a Gaussian distribution with standard deviation sigma, using the formula given above.
pub fn gaussian_pdf(x: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_gaussian_pdf(x, sigma) }
}

pub fn gaussian_ziggurat(r: &mut Rng, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_gaussian_ziggurat(ffi::FFI::unwrap_unique(r), sigma) }
}

/// This function computes a Gaussian random variate using the alternative Marsaglia-Tsang ziggurat and Kinderman-Monahan-Leva ratio methods.
/// The Ziggurat algorithm is the fastest available algorithm in most cases.
pub fn gaussian_ratio_method(r: &mut Rng, sigma: f64) -> f64 {
    unsafe { ffi::gsl_ran_gaussian_ratio_method(ffi::FFI::unwrap_unique(r), sigma) }
}

/// This function computes results for the unit Gaussian distribution.
/// They are equivalent to the functions above with a standard deviation of one, sigma = 1.
pub fn ugaussian(r: &mut Rng) -> f64 {
    unsafe { ffi::gsl_ran_ugaussian(ffi::FFI::unwrap_unique(r)) }
}

/// This function computes results for the unit Gaussian distribution.
/// They are equivalent to the functions above with a standard deviation of one, sigma = 1.
pub fn ugaussian_pdf(x: f64) -> f64 {
    unsafe { ffi::gsl_ran_ugaussian_pdf(x) }
}

/// This function computes results for the unit Gaussian distribution.
/// They are equivalent to the functions above with a standard deviation of one, sigma = 1.
pub fn ugaussian_ratio_method(r: &mut Rng) -> f64 {
    unsafe { ffi::gsl_ran_ugaussian_ratio_method(ffi::FFI::unwrap_unique(r)) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_P(x: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_gaussian_P(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_Q(x: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_gaussian_Q(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_Pinv(P: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_gaussian_Pinv(P, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_Qinv(Q: f64, sigma: f64) -> f64 {
    unsafe { ffi::gsl_cdf_gaussian_Qinv(Q, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_P(x: f64) -> f64 {
    unsafe { ffi::gsl_cdf_ugaussian_P(x) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_Q(x: f64) -> f64 {
    unsafe { ffi::gsl_cdf_ugaussian_Q(x) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_Pinv(P: f64) -> f64 {
    unsafe { ffi::gsl_cdf_ugaussian_Pinv(P) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_Qinv(Q: f64) -> f64 {
    unsafe { ffi::gsl_cdf_ugaussian_Qinv(Q) }
}
