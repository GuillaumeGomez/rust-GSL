//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::Value;
use std::mem::MaybeUninit;

/// This routine computes the Airy function Ai(x) with an accuracy specified by mode.
pub fn Ai(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Ai(x, mode.into()) }
}

/// This routine computes the Airy function Ai(x) with an accuracy specified by mode.
pub fn Ai_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Ai_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the Airy function Bi(x) with an accuracy specified by mode.
pub fn Bi(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Bi(x, mode.into()) }
}

/// This routine computes the Airy function Bi(x) with an accuracy specified by mode.
pub fn Bi_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Bi_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes a scaled version of the Airy function S_A(x) Ai(x). For x>0 the scaling factor S_A(x) is \exp(+(2/3) x^(3/2)), and is 1 for x<0.
pub fn Ai_scaled(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Ai_scaled(x, mode.into()) }
}

/// This routine computes a scaled version of the Airy function S_A(x) Ai(x). For x>0 the scaling factor S_A(x) is \exp(+(2/3) x^(3/2)), and is 1 for x<0.
pub fn Ai_scaled_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Ai_scaled_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes a scaled version of the Airy function S_B(x) Bi(x). For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
pub fn Bi_scaled(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Bi_scaled(x, mode.into()) }
}

/// This routine computes a scaled version of the Airy function S_B(x) Bi(x). For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
pub fn Bi_scaled_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Bi_scaled_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the Airy function derivative Ai'(x) with an accuracy specified by mode.
pub fn Ai_deriv(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Ai_deriv(x, mode.into()) }
}

/// This routine computes the Airy function derivative Ai'(x) with an accuracy specified by mode.
pub fn Ai_deriv_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Ai_deriv_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the Airy function derivative Bi'(x) with an accuracy specified by mode.
pub fn Bi_deriv(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Bi_deriv(x, mode.into()) }
}

/// This routine computes the Airy function derivative Bi'(x) with an accuracy specified by mode.
pub fn Bi_deriv_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Bi_deriv_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the scaled Airy function derivative S_A(x) Ai'(x). For x>0 the scaling factor S_A(x) is \exp(+(2/3) x^(3/2)), and is 1 for x<0.
pub fn Ai_deriv_scaled(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Ai_deriv_scaled(x, mode.into()) }
}

/// This routine computes the scaled Airy function derivative S_A(x) Ai'(x). For x>0 the scaling factor S_A(x) is \exp(+(2/3) x^(3/2)), and is 1 for x<0.
pub fn Ai_deriv_scaled_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Ai_deriv_scaled_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the scaled Airy function derivative S_B(x) Bi'(x). For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
pub fn Bi_deriv_scaled(x: f64, mode: ::Mode) -> f64 {
    unsafe { sys::gsl_sf_airy_Bi_deriv_scaled(x, mode.into()) }
}

/// This routine computes the scaled Airy function derivative S_B(x) Bi'(x). For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
pub fn Bi_deriv_scaled_e(x: f64, mode: ::Mode) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_Bi_deriv_scaled_e(x, mode.into(), result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the location of the s-th zero of the Airy function Ai(x).
pub fn zero_Ai(s: u32) -> f64 {
    unsafe { sys::gsl_sf_airy_zero_Ai(s) }
}

/// This routine computes the location of the s-th zero of the Airy function Ai(x).
pub fn zero_Ai_e(s: u32) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_zero_Ai_e(s, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the location of the s-th zero of the Airy function Bi(x).
pub fn zero_Bi(s: u32) -> f64 {
    unsafe { sys::gsl_sf_airy_zero_Bi(s) }
}

/// This routine computes the location of the s-th zero of the Airy function Bi(x).
pub fn zero_Bi_e(s: u32) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_zero_Bi_e(s, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the location of the s-th zero of the Airy function derivative Ai'(x).
pub fn zero_Ai_deriv(s: u32) -> f64 {
    unsafe { sys::gsl_sf_airy_zero_Ai_deriv(s) }
}

/// This routine computes the location of the s-th zero of the Airy function derivative Ai'(x).
pub fn zero_Ai_deriv_e(s: u32) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_zero_Ai_deriv_e(s, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the location of the s-th zero of the Airy function derivative Bi'(x).
pub fn zero_Bi_deriv(s: u32) -> f64 {
    unsafe { sys::gsl_sf_airy_zero_Bi_deriv(s) }
}

/// This routine computes the location of the s-th zero of the Airy function derivative Bi'(x).
pub fn zero_Bi_deriv_e(s: u32) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_airy_zero_Bi_deriv_e(s, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}
