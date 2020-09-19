//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a Gaussian tail distribution with standard deviation sigma and lower limit a, using the formula given above.
pub fn gaussian_tail_pdf(x: f64, a: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_ran_gaussian_tail_pdf(x, a, sigma) }
}

/// This function computes results for the tail of a unit Gaussian distribution. They are equivalent to the functions above with a standard deviation of one, sigma = 1.
pub fn ugaussian_tail_pdf(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_ran_ugaussian_tail_pdf(x, a) }
}
