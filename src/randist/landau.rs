//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for the Landau distribution using an approximation to the formula given above.
#[doc(alias = "gsl_ran_landau_pdf")]
pub fn landau_pdf(x: f64) -> f64 {
    unsafe { sys::gsl_ran_landau_pdf(x) }
}
