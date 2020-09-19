//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability p(k) of obtaining k from a binomial distribution with parameters p and n, using the formula given above.
pub fn binomial_pdf(k: u32, p: f64, n: u32) -> f64 {
    unsafe { sys::gsl_ran_binomial_pdf(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the binomial distribution with parameters p and n.
pub fn binomial_P(k: u32, p: f64, n: u32) -> f64 {
    unsafe { sys::gsl_cdf_binomial_P(k, p, n) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the binomial distribution with parameters p and n.
pub fn binomial_Q(k: u32, p: f64, n: u32) -> f64 {
    unsafe { sys::gsl_cdf_binomial_Q(k, p, n) }
}
