//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability P(n_1, n_2, ..., n_K) of sampling `n[K]` from a
/// multinomial distribution with parameters `p[K]`, using the formula given above.
#[doc(alias = "gsl_ran_multinomial_pdf")]
pub fn multinomial_pdf(p: &[f64], n: &[u32]) -> f64 {
    unsafe { sys::gsl_ran_multinomial_pdf(p.len() as _, p.as_ptr(), n.as_ptr()) }
}

/// This function returns the logarithm of the probability for the multinomial
/// distribution P(n_1, n_2, ..., n_K) with parameters `p[K]`.
#[doc(alias = "gsl_ran_multinomial_lnpdf")]
pub fn multinomial_lnpdf(p: &[f64], n: &[u32]) -> f64 {
    unsafe { sys::gsl_ran_multinomial_lnpdf(p.len() as _, p.as_ptr(), n.as_ptr()) }
}
