//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(x) at x for a Rayleigh tail distribution with scale parameter sigma and lower limit a, using the formula given above.
pub fn rayleigh_tail_pdf(x: f64, a: f64, sigma: f64) -> f64 {
    unsafe { sys::gsl_ran_rayleigh_tail_pdf(x, a, sigma) }
}
