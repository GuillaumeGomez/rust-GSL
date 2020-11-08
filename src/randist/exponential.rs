//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for an exponential distribution with mean mu, using the formula given above.
pub fn exponential_pdf(x: f64, mu: f64) -> f64 {
    unsafe { sys::gsl_ran_exponential_pdf(x, mu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the exponential distribution with mean mu.
pub fn exponential_P(x: f64, mu: f64) -> f64 {
    unsafe { sys::gsl_cdf_exponential_P(x, mu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the exponential distribution with mean mu.
pub fn exponential_Q(x: f64, mu: f64) -> f64 {
    unsafe { sys::gsl_cdf_exponential_Q(x, mu) }
}

pub fn exponential_Pinv(P: f64, mu: f64) -> f64 {
    unsafe { sys::gsl_cdf_exponential_Pinv(P, mu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the exponential distribution with mean mu.
pub fn exponential_Qinv(Q: f64, mu: f64) -> f64 {
    unsafe { sys::gsl_cdf_exponential_Qinv(Q, mu) }
}
