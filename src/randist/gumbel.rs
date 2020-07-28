//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod type_1 {
    use ffi;
    use types::Rng;

    /// This function returns a random variate from the Type-1 Gumbel distribution. The Type-1 Gumbel distribution function is,
    ///
    /// p(x) dx = a b \exp(-(b \exp(-ax) + ax)) dx
    ///
    /// for -\infty < x < \infty.
    pub fn gumbel1(r: &mut Rng, a: f64, b: f64) -> f64 {
        unsafe { ffi::randist::gsl_ran_gumbel1(ffi::FFI::unwrap_unique(r), a, b) }
    }

    /// This function computes the probability density p(x) at x for a Type-1 Gumbel distribution with parameters a and b, using the formula given above.
    pub fn gumbel1_pdf(x: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::randist::gsl_ran_gumbel1_pdf(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    pub fn gumbel1_P(x: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel1_P(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    pub fn gumbel1_Q(x: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel1_Q(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    pub fn gumbel1_Pinv(P: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel1_Pinv(P, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-1 Gumbel distribution with parameters a and b.
    pub fn gumbel1_Qinv(Q: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel1_Qinv(Q, a, b) }
    }
}

pub mod type_2 {
    use ffi;
    use types::Rng;

    /// This function returns a random variate from the Type-2 Gumbel distribution. The Type-2 Gumbel distribution function is,
    ///
    /// p(x) dx = a b x^{-a-1} \exp(-b x^{-a}) dx
    ///
    /// for 0 < x < \infty.
    pub fn gumbel2(r: &mut Rng, a: f64, b: f64) -> f64 {
        unsafe { ffi::randist::gsl_ran_gumbel2(ffi::FFI::unwrap_unique(r), a, b) }
    }

    /// This function computes the probability density p(x) at x for a Type-2 Gumbel distribution with parameters a and b, using the formula given above.
    pub fn gumbel2_pdf(x: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::randist::gsl_ran_gumbel2_pdf(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    pub fn gumbel2_P(x: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel2_P(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    pub fn gumbel2_Q(x: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel2_Q(x, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    pub fn gumbel2_Pinv(P: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel2_Pinv(P, a, b) }
    }

    /// This function computes the cumulative distribution functions P(x), Q(x) and their inverses for the Type-2 Gumbel distribution with parameters a and b.
    pub fn gumbel2_Qinv(Q: f64, a: f64, b: f64) -> f64 {
        unsafe { ffi::gsl_cdf_gumbel2_Qinv(Q, a, b) }
    }
}
