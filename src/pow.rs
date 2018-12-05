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
    fn pow2(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow3(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow4(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow5(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow6(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow7(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow8(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn pow9(&self) -> Self;
}

impl Pow for f64 {
    fn pow_int(&self, n: i32) -> f64 {
        unsafe { ::ffi::gsl_pow_int(*self, n) }
    }

    fn pow_uint(&self, n: u32) -> f64 {
        unsafe { ::ffi::gsl_pow_uint(*self, n) }
    }

    fn pow2(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_2(*self) }
    }

    fn pow3(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_3(*self) }
    }

    fn pow4(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_4(*self) }
    }

    fn pow5(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_5(*self) }
    }

    fn pow6(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_6(*self) }
    }

    fn pow7(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_7(*self) }
    }

    fn pow8(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_8(*self) }
    }

    fn pow9(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_9(*self) }
    }
}
