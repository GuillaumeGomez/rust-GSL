//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a Weibull distribution with scale a and exponent b, using the formula given above.
pub fn weibull_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_weibull_pdf(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weibull_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_weibull_P(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weibull_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_weibull_Q(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weibull_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_weibull_Pinv(P, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Weibull distribution with scale a and exponent b.
pub fn weibull_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_weibull_Qinv(Q, a, b) }
}
