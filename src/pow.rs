//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub trait Pow {
    /// This routine computes the power x^n for integer n. The power is computed efficiently—for example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
    /// A version of this function which also computes the numerical error in the result is available as gsl_sf_pow_int_e.
    fn pow_int(&self, n: i32) -> Self;
    /// This routine computes the power x^n for integer n. The power is computed efficiently—for example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
    /// A version of this function which also computes the numerical error in the result is available as gsl_sf_pow_int_e.
    fn pow_uint(&self, n: u32) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_2(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_3(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_4(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_5(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_6(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_7(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_8(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow_9(&self) -> Self;
}

impl Pow for f64 {
    #[doc(alias = "gsl_pow_int")]
    fn pow_int(&self, n: i32) -> f64 {
        unsafe { ::sys::gsl_pow_int(*self, n) }
    }

    #[doc(alias = "gsl_pow_uint")]
    fn pow_uint(&self, n: u32) -> f64 {
        unsafe { ::sys::gsl_pow_uint(*self, n) }
    }

    #[doc(alias = "gsl_pow_2")]
    fn pow_2(&self) -> f64 {
        unsafe { ::sys::gsl_pow_2(*self) }
    }

    #[doc(alias = "gsl_pow_3")]
    fn pow_3(&self) -> f64 {
        unsafe { ::sys::gsl_pow_3(*self) }
    }

    #[doc(alias = "gsl_pow_4")]
    fn pow_4(&self) -> f64 {
        unsafe { ::sys::gsl_pow_4(*self) }
    }

    #[doc(alias = "gsl_pow_5")]
    fn pow_5(&self) -> f64 {
        unsafe { ::sys::gsl_pow_5(*self) }
    }

    #[doc(alias = "gsl_pow_6")]
    fn pow_6(&self) -> f64 {
        unsafe { ::sys::gsl_pow_6(*self) }
    }

    #[doc(alias = "gsl_pow_7")]
    fn pow_7(&self) -> f64 {
        unsafe { ::sys::gsl_pow_7(*self) }
    }

    #[doc(alias = "gsl_pow_8")]
    fn pow_8(&self) -> f64 {
        unsafe { ::sys::gsl_pow_8(*self) }
    }

    #[doc(alias = "gsl_pow_9")]
    fn pow_9(&self) -> f64 {
        unsafe { ::sys::gsl_pow_9(*self) }
    }
}
