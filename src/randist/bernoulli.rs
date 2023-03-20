//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability p(k) of obtaining k from a Bernoulli distribution with probability parameter p, using the formula given above.
#[doc(alias = "gsl_ran_bernoulli_pdf")]
pub unsafe fn bernoulli_pdf(x: u32, p: f64) -> f64 {
    sys::gsl_ran_bernoulli_pdf(x, p)
}
