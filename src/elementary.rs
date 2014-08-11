/*
* A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
*/

use ffi;

/// This function computes the value of __log(1+x)__ in a way that is accurate for small x. It provides an alternative to the BSD math function log1p(x).
pub fn log1p(x: f64) -> f64 {
    unsafe { ffi::gsl_log1p(x) }
}

/// This function computes the value of __exp(x)-1__ in a way that is accurate for small x. It provides an alternative to the BSD math function expm1(x).
pub fn expm1(x: f64) -> f64 {
    unsafe { ffi::gsl_expm1(x) }
}

/// This function computes the value of __sqrt{x^2 + y^2}__ in a way that avoids overflow. It provides an alternative to the BSD math function hypot(x,y).
pub fn hypot(x: f64, y: f64) -> f64 {
    unsafe { ffi::gsl_hypot(x, y) }
}

/// This function computes the value of __sqrt{x^2 + y^2 + z^2}__ in a way that avoids overflow.
pub fn hypot3(x: f64, y: f64, z: f64) -> f64 {
    unsafe { ffi::gsl_hypot3(x, y, z) }
}

/// This function computes the value of __arccosh(x)__. It provides an alternative to the standard math function acosh(x).
pub fn acosh(x: f64) -> f64 {
    unsafe { ffi::gsl_acosh(x) }
}

/// This function computes the value of __arcsinh(x)__. It provides an alternative to the standard math function asinh(x).
pub fn asinh(x: f64) -> f64 {
    unsafe { ffi::gsl_asinh(x) }
}

/// This function computes the value of __arctanh(x)__. It provides an alternative to the standard math function atanh(x).
pub fn atanh(x: f64) -> f64 {
    unsafe { ffi::gsl_atanh(x) }
}

/// This function computes the value of __x * 2^e__. It provides an alternative to the standard math function ldexp(x,e).
pub fn ldexp(x: f64, e: i32) -> f64 {
    unsafe { ffi::gsl_ldexp(x, e) }
}

/// This function splits the number x into its normalized fraction f and exponent e, such that x = f * 2^e and 0.5 <= f < 1. The function returns f and stores the exponent in e.
/// If x is zero, both f and e are set to zero. This function provides an alternative to the standard math function frexp(x, e).
pub fn frexp(x: f64, e: &mut i32) -> f64 {
    unsafe { ffi::gsl_frexp(x, e) }
}