//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a beta distribution with parameters a and b, using the formula given above.
pub fn beta_pdf(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_ran_beta_pdf(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_P(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_beta_P(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_Q(x: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_beta_Q(x, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_Pinv(P: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_beta_Pinv(P, a, b) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the beta distribution with parameters a and b.
pub fn beta_Qinv(Q: f64, a: f64, b: f64) -> f64 {
    unsafe { sys::gsl_cdf_beta_Qinv(Q, a, b) }
}
