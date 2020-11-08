//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability p(k) of obtaining k from a geometric distribution with probability parameter p, using the formula given above.
pub fn geometric_pdf(k: u32, p: f64) -> f64 {
    unsafe { sys::gsl_ran_geometric_pdf(k, p) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the geometric distribution with parameter p.
pub fn geometric_P(k: u32, p: f64) -> f64 {
    unsafe { sys::gsl_cdf_geometric_P(k, p) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the geometric distribution with parameter p.
pub fn geometric_Q(k: u32, p: f64) -> f64 {
    unsafe { sys::gsl_cdf_geometric_Q(k, p) }
}
