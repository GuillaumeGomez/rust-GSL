//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a Cauchy distribution with scale parameter a, using the formula given above.
#[doc(alias = "gsl_ran_cauchy_pdf")]
pub fn cauchy_pdf(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_ran_cauchy_pdf(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
#[doc(alias = "gsl_cdf_cauchy_P")]
pub fn cauchy_P(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_cauchy_P(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
#[doc(alias = "gsl_cdf_cauchy_Q")]
pub fn cauchy_Q(x: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_cauchy_Q(x, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
#[doc(alias = "gsl_cdf_cauchy_Pinv")]
pub fn cauchy_Pinv(P: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_cauchy_Pinv(P, a) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Cauchy distribution with scale parameter a.
#[doc(alias = "gsl_cdf_cauchy_Qinv")]
pub fn cauchy_Qinv(Q: f64, a: f64) -> f64 {
    unsafe { sys::gsl_cdf_cauchy_Qinv(Q, a) }
}
