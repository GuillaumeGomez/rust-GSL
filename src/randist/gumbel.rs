//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod type_1 {
    /// This function computes the probability density p(x) at x for a Type-1 Gumbel distribution with parameters a and b, using the formula given above.
    #[doc(alias = "gsl_ran_gumbel1_pdf")]
    pub fn gumbel1_pdf(x: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_ran_gumbel1_pdf(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel1_P")]
    pub fn gumbel1_P(x: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel1_P(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel1_Q")]
    pub fn gumbel1_Q(x: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel1_Q(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel1_Pinv")]
    pub fn gumbel1_Pinv(P: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel1_Pinv(P, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel1_Qinv")]
    pub fn gumbel1_Qinv(Q: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel1_Qinv(Q, a, b) }
    }
}

pub mod type_2 {
    /// This function computes the probability density p(x) at x for a Type-2 Gumbel distribution with parameters a and b, using the formula given above.
    #[doc(alias = "gsl_ran_gumbel2_pdf")]
    pub fn gumbel2_pdf(x: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_ran_gumbel2_pdf(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel2_P")]
    pub fn gumbel2_P(x: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel2_P(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel2_Q")]
    pub fn gumbel2_Q(x: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel2_Q(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel2_Pinv")]
    pub fn gumbel2_Pinv(P: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel2_Pinv(P, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    #[doc(alias = "gsl_cdf_gumbel2_Qinv")]
    pub fn gumbel2_Qinv(Q: f64, a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_cdf_gumbel2_Qinv(Q, a, b) }
    }
}
