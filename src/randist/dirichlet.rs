//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This function computes the probability density p(\theta_1, ... , \theta_K) at `theta[K]`
/// for a Dirichlet distribution with parameters `alpha[K]`, using the formula given above.
pub fn dirichlet_pdf(alpha: &[f64], theta: &[f64]) -> f64 {
    unsafe { sys::gsl_ran_dirichlet_pdf(alpha.len() as _, alpha.as_ptr(), theta.as_ptr()) }
}

/// This function computes the logarithm of the probability density p(\theta_1, ... , \theta_K)
/// for a Dirichlet distribution with parameters `alpha[K]`.
pub fn dirichlet_lnpdf(alpha: &[f64], theta: &[f64]) -> f64 {
    unsafe { sys::gsl_ran_dirichlet_lnpdf(alpha.len() as _, alpha.as_ptr(), theta.as_ptr()) }
}
