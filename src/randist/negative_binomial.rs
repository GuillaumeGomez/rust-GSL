//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability p(k) of obtaining k from a negative binomial distribution with parameters p and n, using the formula given above.
pub fn negative_binomial_pdf(k: u32, p: f64, n: f64) -> f64 {
    unsafe { sys::gsl_ran_negative_binomial_pdf(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the negative binomial distribution with parameters p and n.
pub fn negative_binomial_P(k: u32, p: f64, n: f64) -> f64 {
    unsafe { sys::gsl_cdf_negative_binomial_P(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the negative binomial distribution with parameters p and n.
pub fn negative_binomial_Q(k: u32, p: f64, n: f64) -> f64 {
    unsafe { sys::gsl_cdf_negative_binomial_Q(k, p, n) }
}
