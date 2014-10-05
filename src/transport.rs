//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The transport functions J(n,x) are defined by the integral representations J(n,x) := \int_0^x dt t^n e^t /(e^t - 1)^2.

use ffi;
use std::mem::zeroed;
use enums;

/// This routine computes the transport function J(2,x).
pub fn transport_2(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_transport_2(x) }
}

/// This routine computes the transport function J(2,x).
pub fn transport_2_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_transport_2_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the transport function J(3,x).
pub fn transport_3(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_transport_3(x) }
}

/// This routine computes the transport function J(3,x).
pub fn transport_3_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_transport_3_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the transport function J(4,x).
pub fn transport_4(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_transport_4(x) }
}

/// This routine computes the transport function J(4,x).
pub fn transport_4_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_transport_4_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the transport function J(5,x).
pub fn transport_5(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_transport_5(x) }
}

/// This routine computes the transport function J(5,x).
pub fn transport_5_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_transport_5_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}