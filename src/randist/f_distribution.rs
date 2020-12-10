//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The F-distribution arises in statistics. If Y_1 and Y_2 are chi-squared deviates with \nu_1 and \nu_2 degrees of freedom then the ratio,

X = { (Y_1 / \nu_1) \over (Y_2 / \nu_2) }

has an F-distribution F(x;\nu_1,\nu_2).
!*/

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the F-distribution with nu1 and nu2 degrees of freedom.
#[doc(alias = "gsl_ran_fdist_pdf")]
pub fn fdist_pdf(x: f64, nu1: f64, nu2: f64) -> f64 {
    unsafe { sys::gsl_ran_fdist_pdf(x, nu1, nu2) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the F-distribution with nu1 and nu2 degrees of freedom.
#[doc(alias = "gsl_cdf_fdist_P")]
pub fn fdist_P(x: f64, nu1: f64, nu2: f64) -> f64 {
    unsafe { sys::gsl_cdf_fdist_P(x, nu1, nu2) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the F-distribution with nu1 and nu2 degrees of freedom.
#[doc(alias = "gsl_cdf_fdist_Q")]
pub fn fdist_Q(x: f64, nu1: f64, nu2: f64) -> f64 {
    unsafe { sys::gsl_cdf_fdist_Q(x, nu1, nu2) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the F-distribution with nu1 and nu2 degrees of freedom.
#[doc(alias = "gsl_cdf_fdist_Pinv")]
pub fn fdist_Pinv(P: f64, nu1: f64, nu2: f64) -> f64 {
    unsafe { sys::gsl_cdf_fdist_Pinv(P, nu1, nu2) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the F-distribution with nu1 and nu2 degrees of freedom.
#[doc(alias = "gsl_cdf_fdist_Qinv")]
pub fn fdist_Qinv(Q: f64, nu1: f64, nu2: f64) -> f64 {
    unsafe { sys::gsl_cdf_fdist_Qinv(Q, nu1, nu2) }
}
