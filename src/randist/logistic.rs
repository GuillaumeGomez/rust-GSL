//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a logistic distribution with scale parameter a, using the formula given above.
#[doc(alias = "gsl_ran_logistic_pdf")]
pub fn logistic_pdf(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_ran_logistic_pdf(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
#[doc(alias = "gsl_cdf_logistic_P")]
pub fn logistic_P(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_logistic_P(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
#[doc(alias = "gsl_cdf_logistic_Q")]
pub fn logistic_Q(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_logistic_Q(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
#[doc(alias = "gsl_cdf_logistic_Pinv")]
pub fn logistic_Pinv(P: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_logistic_Pinv(P, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the logistic distribution with scale parameter a.
#[doc(alias = "gsl_cdf_logistic_Qinv")]
pub fn logistic_Qinv(Q: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_logistic_Qinv(Q, a) }
}
