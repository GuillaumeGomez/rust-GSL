//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Information on the properties of the Logarithm function can be found in Abramowitz & Stegun, Chapter 4.

use crate::Value;
use std::mem::MaybeUninit;

/// This routine computes the logarithm of x, \log(x), for x > 0.
pub fn log(x: f64) -> f64 {
    unsafe { sys::gsl_sf_log(x) }
}

/// This routine computes the logarithm of x, \log(x), for x > 0.
pub fn log_e(x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_log_e(x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the logarithm of the magnitude of x, \log(|x|), for x \ne 0.
pub fn log_abs(x: f64) -> f64 {
    unsafe { sys::gsl_sf_log_abs(x) }
}

/// This routine computes the logarithm of the magnitude of x, \log(|x|), for x \ne 0.
pub fn log_abs_e(x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_log_abs_e(x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the complex logarithm of z = z_r + i z_i.
/// The results are returned as lnr, theta such that \exp(lnr + i \theta) = z_r + i z_i, where \theta lies in the range [-\pi,\pi].
pub fn complex_log_e(zr: f64, zi: f64) -> (Value, ::types::Result, ::types::Result) {
    let mut lnr = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let mut theta = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_complex_log_e(zr, zi, lnr.as_mut_ptr(), theta.as_mut_ptr()) };

    (
        Value::from(ret),
        unsafe { lnr.assume_init() }.into(),
        unsafe { theta.assume_init() }.into(),
    )
}

/// This routine computes \log(1 + x) for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx(x: f64) -> f64 {
    unsafe { sys::gsl_sf_log_1plusx(x) }
}

/// This routine computes \log(1 + x) for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx_e(x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_log_1plusx_e(x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes \log(1 + x) - x for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx_mx(x: f64) -> f64 {
    unsafe { sys::gsl_sf_log_1plusx_mx(x) }
}

/// This routine computes \log(1 + x) - x for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx_mx_e(x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_log_1plusx_mx_e(x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}
