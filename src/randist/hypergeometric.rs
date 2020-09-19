//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability p(k) of obtaining k from a hypergeometric distribution with parameters n1, n2, t, using the formula given above.
pub fn hypergeometric_pdf(k: u32, n1: u32, n2: u32, t: u32) -> f64 {
    unsafe { sys::gsl_ran_hypergeometric_pdf(k, n1, n2, t) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the hypergeometric distribution with parameters n1, n2 and t.
pub fn hypergeometric_P(k: u32, n1: u32, n2: u32, t: u32) -> f64 {
    unsafe { sys::gsl_cdf_hypergeometric_P(k, n1, n2, t) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the hypergeometric distribution with parameters n1, n2 and t.
pub fn hypergeometric_Q(k: u32, n1: u32, n2: u32, t: u32) -> f64 {
    unsafe { sys::gsl_cdf_hypergeometric_Q(k, n1, n2, t) }
}
