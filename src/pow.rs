//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// This routine computes the power x^n for integer n. The power is computed efficiently—for example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
/// A version of this function which also computes the numerical error in the result is available as gsl_sf_pow_int_e.
pub fn _int(x: f64, n: i32) -> f64 {
    unsafe { ::ffi::gsl_pow_int(x, n) }
}

/// This routine computes the power x^n for integer n. The power is computed efficiently—for example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
/// A version of this function which also computes the numerical error in the result is available as gsl_sf_pow_int_e.
pub fn _uint(x: f64, n: u32) -> f64 {
    unsafe { ::ffi::gsl_pow_uint(x, n) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _2(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_2(x) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _3(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_3(x) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _4(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_4(x) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _5(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_5(x) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _6(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_6(x) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _7(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_7(x) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _8(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_8(x) }
}

/// This function can be used to compute small integer powers x^2, x^3, etc. efficiently.
pub fn _9(x: f64) -> f64 {
    unsafe { ::ffi::gsl_pow_9(x) }
}