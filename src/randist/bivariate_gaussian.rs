//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function generates a pair of correlated Gaussian variates, with mean zero, correlation coefficient rho and standard deviations sigma_x and sigma_y in the x and y directions.
/// The probability distribution for bivariate Gaussian random variates is,
/// 
/// p(x,y) dx dy = {1 \over 2 \pi \sigma_x \sigma_y \sqrt{1-\rho^2}} \exp (-(x^2/\sigma_x^2 + y^2/\sigma_y^2 - 2 \rho x y/(\sigma_x\sigma_y))/2(1-\rho^2)) dx dy
/// 
/// for x,y in the range -\infty to +\infty. The correlation coefficient rho should lie between 1 and -1.
pub fn gaussian_tail(r: &Rng, sigma_x: f64, sigma_y: f64, rho: f64, x: &mut f64, y: &mut f64) {
    unsafe { ffi::gsl_ran_bivariate_gaussian(ffi::FFI::unwrap(r), sigma_x, sigma_y, rho, x, y) }
}

/// This function computes the probability density p(x,y) at (x,y) for a bivariate Gaussian distribution with standard deviations sigma_x, sigma_y and correlation coefficient rho, using the formula given above.
pub fn gaussian_tail_pdf(x: f64, y: f64, sigma_x: f64, sigma_y: f64, rho: f64) -> f64 {
    unsafe { ffi::gsl_ran_bivariate_gaussian_pdf(x, y, sigma_x, sigma_y, rho) }
}