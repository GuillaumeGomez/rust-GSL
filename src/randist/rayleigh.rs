//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a Rayleigh distribution with scale parameter sigma, using the formula given above.
#[doc(alias = "gsl_ran_rayleigh_pdf")]
pub fn rayleigh_pdf(x: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_ran_rayleigh_pdf(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
#[doc(alias = "gsl_cdf_rayleigh_P")]
pub fn rayleigh_P(x: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_rayleigh_P(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
#[doc(alias = "gsl_cdf_rayleigh_Q")]
pub fn rayleigh_Q(x: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_rayleigh_Q(x, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
#[doc(alias = "gsl_cdf_rayleigh_Pinv")]
pub fn rayleigh_Pinv(P: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_rayleigh_Pinv(P, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Rayleigh distribution with scale parameter sigma.
#[doc(alias = "gsl_cdf_rayleigh_Qinv")]
pub fn rayleigh_Qinv(Q: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_rayleigh_Qinv(Q, sigma) }
}
