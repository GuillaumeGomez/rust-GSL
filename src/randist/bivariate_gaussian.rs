//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x,y) at (x,y) for a bivariate Gaussian distribution with standard deviations sigma_x, sigma_y and correlation coefficient rho, using the formula given above.
pub fn gaussian_tail_pdf(x: f64, y: f64, sigma_x: f64, sigma_y: f64, rho: f64) -> f64 {
    unsafe { sys::gsl_ran_bivariate_gaussian_pdf(x, y, sigma_x, sigma_y, rho) }
}
