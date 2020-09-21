//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The transport functions J(n,x) are defined by the integral representations J(n,x) := \int_0^x dt t^n e^t /(e^t - 1)^2.

use enums;
use std::mem::MaybeUninit;

/// This routine computes the transport function J(2,x).
pub fn transport_2(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_2(x) }
}

/// This routine computes the transport function J(2,x).
pub fn transport_2_e(x: f64) -> Result<::types::Result, enums::Value> {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_transport_2_e(x, result.as_mut_ptr()) };

    result!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the transport function J(3,x).
pub fn transport_3(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_3(x) }
}

/// This routine computes the transport function J(3,x).
pub fn transport_3_e(x: f64) -> Result<::types::Result, enums::Value> {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_transport_3_e(x, result.as_mut_ptr()) };

    result!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the transport function J(4,x).
pub fn transport_4(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_4(x) }
}

/// This routine computes the transport function J(4,x).
pub fn transport_4_e(x: f64) -> Result<::types::Result, enums::Value> {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_transport_4_e(x, result.as_mut_ptr()) };

    result!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the transport function J(5,x).
pub fn transport_5(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_5(x) }
}

/// This routine computes the transport function J(5,x).
pub fn transport_5_e(x: f64) -> Result<::types::Result, enums::Value> {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_transport_5_e(x, result.as_mut_ptr()) };

    result!(ret, unsafe { result.assume_init() }.into())
}
