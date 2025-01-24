//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{types, Value};
use std::mem::MaybeUninit;

/// This routine computes the exponential integral E_1(x),
///
/// E_1(x) := \Re \int_1^\infty dt \exp(-xt)/t.
#[doc(alias = "gsl_sf_expint_E1")]
pub fn E1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_expint_E1(x) }
}

/// This routine computes the exponential integral E_1(x),
///
/// E_1(x) := \Re \int_1^\infty dt \exp(-xt)/t.
#[doc(alias = "gsl_sf_expint_E1_e")]
pub fn E1_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_expint_E1_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the second-order exponential integral E_2(x),
///
/// E_2(x) := \Re \int_1^\infty dt \exp(-xt)/t^2.
#[doc(alias = "gsl_sf_expint_E2")]
pub fn E2(x: f64) -> f64 {
    unsafe { sys::gsl_sf_expint_E2(x) }
}

/// This routine computes the second-order exponential integral E_2(x),
///
/// E_2(x) := \Re \int_1^\infty dt \exp(-xt)/t^2.
#[doc(alias = "gsl_sf_expint_E2_e")]
pub fn E2_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_expint_E2_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the exponential integral E_n(x) of order n,
///
/// E_n(x) := \Re \int_1^\infty dt \exp(-xt)/t^n.
#[doc(alias = "gsl_sf_expint_En")]
pub fn En(n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_expint_En(n, x) }
}

/// This routine computes the exponential integral E_n(x) of order n,
///
/// E_n(x) := \Re \int_1^\infty dt \exp(-xt)/t^n.
#[doc(alias = "gsl_sf_expint_En_e")]
pub fn En_e(n: i32, x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_expint_En_e(n, x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the exponential integral Ei(x),
///
/// Ei(x) := - PV(\int_{-x}^\infty dt \exp(-t)/t)
///
/// where PV denotes the principal value of the integral.
#[doc(alias = "gsl_sf_expint_Ei")]
pub fn Ei(x: f64) -> f64 {
    unsafe { sys::gsl_sf_expint_Ei(x) }
}

/// This routine computes the exponential integral Ei(x),
///
/// Ei(x) := - PV(\int_{-x}^\infty dt \exp(-t)/t)
///
/// where PV denotes the principal value of the integral.
#[doc(alias = "gsl_sf_expint_Ei_e")]
pub fn Ei_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_expint_Ei_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the integral Shi(x) = \int_0^x dt \sinh(t)/t.
#[doc(alias = "gsl_sf_Shi")]
pub fn Shi(x: f64) -> f64 {
    unsafe { sys::gsl_sf_Shi(x) }
}

/// This routine computes the integral Shi(x) = \int_0^x dt \sinh(t)/t.
#[doc(alias = "gsl_sf_Shi_e")]
pub fn Shi_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_Shi_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the integral Chi(x) := \Re[ \gamma_E + \log(x) + \int_0^x dt (\cosh(t)-1)/t] , where \gamma_E is the Euler constant (available as the macro M_EULER).
#[doc(alias = "gsl_sf_Chi")]
pub fn Chi(x: f64) -> f64 {
    unsafe { sys::gsl_sf_Chi(x) }
}

/// This routine computes the integral Chi(x) := \Re[ \gamma_E + \log(x) + \int_0^x dt (\cosh(t)-1)/t] , where \gamma_E is the Euler constant (available as the macro M_EULER).
#[doc(alias = "gsl_sf_Chi_e")]
pub fn Chi_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_Chi_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the third-order exponential integral Ei_3(x) = \int_0^xdt \exp(-t^3) for x >= 0.
#[doc(alias = "gsl_sf_expint_3")]
pub fn _3(x: f64) -> f64 {
    unsafe { sys::gsl_sf_expint_3(x) }
}

/// This routine computes the third-order exponential integral Ei_3(x) = \int_0^xdt \exp(-t^3) for x >= 0.
#[doc(alias = "gsl_sf_expint_3_e")]
pub fn _3_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_expint_3_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Sine integral Si(x) = \int_0^x dt \sin(t)/t.
#[doc(alias = "gsl_sf_Si")]
pub fn Si(x: f64) -> f64 {
    unsafe { sys::gsl_sf_Si(x) }
}

/// This routine computes the Sine integral Si(x) = \int_0^x dt \sin(t)/t.
#[doc(alias = "gsl_sf_Si_e")]
pub fn Si_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_Si_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Cosine integral Ci(x) = -\int_x^\infty dt \cos(t)/t for x > 0.
#[doc(alias = "gsl_sf_Ci")]
pub fn Ci(x: f64) -> f64 {
    unsafe { sys::gsl_sf_Ci(x) }
}

/// This routine computes the Cosine integral Ci(x) = -\int_x^\infty dt \cos(t)/t for x > 0.
#[doc(alias = "gsl_sf_Ci_e")]
pub fn Ci_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_Ci_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Arctangent integral, which is defined as AtanInt(x) = \int_0^x dt \arctan(t)/t.
#[doc(alias = "gsl_sf_atanint")]
pub fn atanint(x: f64) -> f64 {
    unsafe { sys::gsl_sf_atanint(x) }
}

/// This routine computes the Arctangent integral, which is defined as AtanInt(x) = \int_0^x dt \arctan(t)/t.
#[doc(alias = "gsl_sf_atanint_e")]
pub fn atanint_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_atanint_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}
