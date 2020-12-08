//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for an exponential power distribution with scale parameter a and exponent b, using the formula given above.
#[doc(alias = "gsl_ran_exppow_pdf")]
pub fn exppow_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_exppow_pdf(x, a, b) }
}

/// This function computes tthe cumulative distribution functions P(x), Q(x) for the exponential power distribution with parameters a and b.
#[doc(alias = "gsl_cdf_exppow_P")]
pub fn exppow_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_exppow_P(x, a, b) }
}

/// This function computes tthe cumulative distribution functions P(x), Q(x) for the exponential power distribution with parameters a and b.
#[doc(alias = "gsl_cdf_exppow_Q")]
pub fn exppow_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_exppow_Q(x, a, b) }
}
