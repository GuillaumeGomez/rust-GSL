//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a Gaussian distribution with standard deviation sigma, using the formula given above.
pub fn gaussian_pdf(x: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_ran_gaussian_pdf(x, sigma) }
}

/// This function computes results for the unit Gaussian distribution.
/// They are equivalent to the functions above with a standard deviation of one, sigma = 1.
pub fn ugaussian_pdf(x: f64) -> f64 {
    unsafe { sys::gsl_ran_ugaussian_pdf(x) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_P(x: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_gaussian_P(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_Q(x: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_gaussian_Q(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_Pinv(P: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_gaussian_Pinv(P, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Gaussian distribution with standard deviation sigma.
pub fn gaussian_Qinv(Q: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_gaussian_Qinv(Q, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_P(x: f64) -> f64 {
    unsafe { sys::gsl_cdf_ugaussian_P(x) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_Q(x: f64) -> f64 {
    unsafe { sys::gsl_cdf_ugaussian_Q(x) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_Pinv(P: f64) -> f64 {
    unsafe { sys::gsl_cdf_ugaussian_Pinv(P) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the unit Gaussian distribution.
pub fn ugaussian_Qinv(Q: f64) -> f64 {
    unsafe { sys::gsl_cdf_ugaussian_Qinv(Q) }
}
