/*
* A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
*/

use ffi;

pub trait Elementary {
    /// This function computes the value of __log(1+x)__ in a way that is accurate for small x. It provides an alternative to the BSD math function log1p(x).
    fn log1p(&self) -> Self;
    /// This function computes the value of __exp(x)-1__ in a way that is accurate for small x. It provides an alternative to the BSD math function expm1(x).
    fn expm1(&self) -> Self;
    /// This function computes the value of __sqrt{x^2 + y^2}__ in a way that avoids overflow. It provides an alternative to the BSD math function hypot(x,y).
    fn hypot(&self, y: f64) -> Self;
    /// This function computes the value of __sqrt{x^2 + y^2 + z^2}__ in a way that avoids overflow.
    fn hypot3(&self, y: f64, z: f64) -> Self;
    /// This function computes the value of __arccosh(x)__. It provides an alternative to the standard math function acosh(x).
    fn acosh(&self) -> Self;
    /// This function computes the value of __arcsinh(x)__. It provides an alternative to the standard math function asinh(x).
    fn asinh(&self) -> Self;
    /// This function computes the value of __arctanh(x)__. It provides an alternative to the standard math function atanh(x).
    fn atanh(&self) -> Self;
    /// This function computes the value of __x * 2^e__. It provides an alternative to the standard math function ldexp(x,e).
    fn ldexp(&self, e: i32) -> Self;
    /// This function splits the number x into its normalized fraction f and exponent e, such that x = f * 2^e and 0.5 <= f < 1. The function returns f and stores the exponent in e.
    /// If x is zero, both f and e are set to zero. This function provides an alternative to the standard math function frexp(x, e).
    fn frexp(&self, e: &mut i32) -> Self;
}

impl Elementary for f64 {
    fn log1p(&self) -> f64 {
        unsafe { ffi::gsl_log1p(*self) }
    }

    fn expm1(&self) -> f64 {
        unsafe { ffi::gsl_expm1(*self) }
    }

    fn hypot(&self, y: f64) -> f64 {
        unsafe { ffi::gsl_hypot(*self, y) }
    }

    fn hypot3(&self, y: f64, z: f64) -> f64 {
        unsafe { ffi::gsl_hypot3(*self, y, z) }
    }

    fn acosh(&self) -> f64 {
        unsafe { ffi::gsl_acosh(*self) }
    }

    fn asinh(&self) -> f64 {
        unsafe { ffi::gsl_asinh(*self) }
    }

    fn atanh(&self) -> f64 {
        unsafe { ffi::gsl_atanh(*self) }
    }

    fn ldexp(&self, e: i32) -> f64 {
        unsafe { ffi::gsl_ldexp(*self, e) }
    }

    fn frexp(&self, e: &mut i32) -> f64 {
        unsafe { ffi::gsl_frexp(*self, e) }
    }
}