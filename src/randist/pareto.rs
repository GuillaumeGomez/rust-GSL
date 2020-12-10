//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a Pareto distribution with exponent a and scale b, using the formula given above.
#[doc(alias = "gsl_ran_pareto_pdf")]
pub fn pareto_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_pareto_pdf(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
#[doc(alias = "gsl_cdf_pareto_P")]
pub fn pareto_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_pareto_P(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
#[doc(alias = "gsl_cdf_pareto_Q")]
pub fn pareto_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_pareto_Q(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
#[doc(alias = "gsl_cdf_pareto_Pinv")]
pub fn pareto_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_pareto_Pinv(P, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Pareto distribution with exponent a and scale b.
#[doc(alias = "gsl_cdf_pareto_Qinv")]
pub fn pareto_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_pareto_Qinv(Q, a, b) }
}
