//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The transport functions J(n,x) are defined by the integral representations J(n,x) := \int_0^x dt t^n e^t /(e^t - 1)^2.

use crate::{types, Error};
use std::mem::MaybeUninit;

/// This routine computes the transport function J(2,x).
#[doc(alias = "gsl_sf_transport_2")]
pub fn transport_2(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_2(x) }
}

/// This routine computes the transport function J(2,x).
#[doc(alias = "gsl_sf_transport_2_e")]
pub fn transport_2_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_transport_2_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the transport function J(3,x).
#[doc(alias = "gsl_sf_transport_3")]
pub fn transport_3(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_3(x) }
}

/// This routine computes the transport function J(3,x).
#[doc(alias = "gsl_sf_transport_3_e")]
pub fn transport_3_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_transport_3_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the transport function J(4,x).
#[doc(alias = "gsl_sf_transport_4")]
pub fn transport_4(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_4(x) }
}

/// This routine computes the transport function J(4,x).
#[doc(alias = "gsl_sf_transport_4_e")]
pub fn transport_4_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_transport_4_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the transport function J(5,x).
#[doc(alias = "gsl_sf_transport_5")]
pub fn transport_5(x: f64) -> f64 {
    unsafe { sys::gsl_sf_transport_5(x) }
}

/// This routine computes the transport function J(5,x).
#[doc(alias = "gsl_sf_transport_5_e")]
pub fn transport_5_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_transport_5_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
