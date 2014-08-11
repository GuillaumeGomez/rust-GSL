//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub trait Pow {
    /// This routine computes the power x^n for integer n. The power is computed efficiently—for example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
    /// A version of this function which also computes the numerical error in the result is available as gsl_sf_pow_int_e.
    fn _int(&self, n: i32) -> Self;
    /// This routine computes the power x^n for integer n. The power is computed efficiently—for example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
    /// A version of this function which also computes the numerical error in the result is available as gsl_sf_pow_int_e.
    fn _uint(&self, n: u32) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _2(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _3(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _4(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _5(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _6(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _7(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _8(&self) -> Self;
    /// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
    fn _9(&self) -> Self;
}

impl Pow for f64 {
    fn _int(&self, n: i32) -> f64 {
        unsafe { ::ffi::gsl_pow_int(*self, n) }
    }

    fn _uint(&self, n: u32) -> f64 {
        unsafe { ::ffi::gsl_pow_uint(*self, n) }
    }

    fn _2(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_2(*self) }
    }

    fn _3(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_3(*self) }
    }

    fn _4(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_4(*self) }
    }

    fn _5(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_5(*self) }
    }

    fn _6(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_6(*self) }
    }

    fn _7(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_7(*self) }
    }

    fn _8(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_8(*self) }
    }

    fn _9(&self) -> f64 {
        unsafe { ::ffi::gsl_pow_9(*self) }
    }
}