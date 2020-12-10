//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a lognormal distribution with parameters zeta and sigma, using the formula given above.
#[doc(alias = "gsl_ran_lognormal_pdf")]
pub fn lognormal_pdf(x: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_ran_lognormal_pdf(x, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
#[doc(alias = "gsl_cdf_lognormal_P")]
pub fn lognormal_P(x: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_P(x, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
#[doc(alias = "gsl_cdf_lognormal_Q")]
pub fn lognormal_Q(x: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_Q(x, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
#[doc(alias = "gsl_cdf_lognormal_Pinv")]
pub fn lognormal_Pinv(P: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_Pinv(P, zeta, sigma) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the lognormal distribution with parameters zeta and sigma.
#[doc(alias = "gsl_cdf_lognormal_Qinv")]
pub fn lognormal_Qinv(Q: f64, zeta: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_cdf_lognormal_Qinv(Q, zeta, sigma) }
}
