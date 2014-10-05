//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Information on the properties of the Logarithm function can be found in Abramowitz & Stegun, Chapter 4.

use std::mem::zeroed;
use enums;
use ffi;

/// This routine computes the logarithm of x, \log(x), for x > 0.
pub fn log(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_log(x) }
}

/// This routine computes the logarithm of x, \log(x), for x > 0.
pub fn log_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_log_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the logarithm of the magnitude of x, \log(|x|), for x \ne 0.
pub fn log_abs(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_log_abs(x) }
}

/// This routine computes the logarithm of the magnitude of x, \log(|x|), for x \ne 0.
pub fn log_abs_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_log_abs_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the complex logarithm of z = z_r + i z_i.
/// The results are returned as lnr, theta such that \exp(lnr + i \theta) = z_r + i z_i, where \theta lies in the range [-\pi,\pi].
pub fn complex_log_e(zr: f64, zi: f64) -> (enums::value::Value, ::types::Result, ::types::Result) {
    let mut lnr = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let mut theta = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_complex_log_e(zr, zi, &mut lnr, &mut theta) };

    (ret, ::types::Result{val: lnr.val, err: lnr.err}, ::types::Result{val: theta.val, err: theta.err})
}

/// This routine computes \log(1 + x) for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_log_1plusx(x) }
}

/// This routine computes \log(1 + x) for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_log_1plusx_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes \log(1 + x) - x for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx_mx(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_log_1plusx_mx(x) }
}

/// This routine computes \log(1 + x) - x for x > -1 using an algorithm that is accurate for small x.
pub fn log_1plusx_mx_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_log_1plusx_mx_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}