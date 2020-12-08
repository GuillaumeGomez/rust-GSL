//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The chi-squared distribution arises in statistics. If Y_i are n independent Gaussian random variates with unit variance then the sum-of-squares,

X_i = \sum_i Y_i^2

has a chi-squared distribution with n degrees of freedom.
!*/

/// This function computes the probability density p(x) at x for a chi-squared distribution with nu degrees of freedom, using the formula given above.
#[doc(alias = "gsl_ran_chisq_pdf")]
pub fn chisq_pdf(x: f64, nu: f64) -> f64 {
    unsafe { sys::gsl_ran_chisq_pdf(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
#[doc(alias = "gsl_cdf_chisq_P")]
pub fn chisq_P(x: f64, nu: f64) -> f64 {
    unsafe { sys::gsl_cdf_chisq_P(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
#[doc(alias = "gsl_cdf_chisq_Q")]
pub fn chisq_Q(x: f64, nu: f64) -> f64 {
    unsafe { sys::gsl_cdf_chisq_Q(x, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
#[doc(alias = "gsl_cdf_chisq_Pinv")]
pub fn chisq_Pinv(P: f64, nu: f64) -> f64 {
    unsafe { sys::gsl_cdf_chisq_Pinv(P, nu) }
}

/// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the chi-squared distribution with nu degrees of freedom.
#[doc(alias = "gsl_cdf_chisq_Qinv")]
pub fn chisq_Qinv(Q: f64, nu: f64) -> f64 {
    unsafe { sys::gsl_cdf_chisq_Qinv(Q, nu) }
}
