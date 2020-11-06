//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The generalized Laguerre polynomials are defined in terms of confluent hypergeometric functions as L^a_n(x) = ((a+1)_n / n!) 1F1(-n,a+1,x), and are sometimes referred to as the associated Laguerre polynomials.
They are related to the plain Laguerre polynomials L_n(x) by L^0_n(x) = L_n(x) and L^k_n(x) = (-1)^k (d^k/dx^k) L_(n+k)(x). For more information see Abramowitz & Stegun, Chapter 22.
!*/

use crate::Value;
use std::mem::MaybeUninit;

/// This function evaluates the generalized Laguerre polynomials L^a_1(x), L^a_2(x), L^a_3(x) using explicit representations.
pub fn laguerre_1(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_laguerre_1(a, x) }
}

/// This function evaluates the generalized Laguerre polynomials L^a_1(x), L^a_2(x), L^a_3(x) using explicit representations.
pub fn laguerre_2(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_laguerre_2(a, x) }
}

/// This function evaluates the generalized Laguerre polynomials L^a_1(x), L^a_2(x), L^a_3(x) using explicit representations.
pub fn laguerre_3(a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_laguerre_3(a, x) }
}

/// This function evaluates the generalized Laguerre polynomials L^a_1(x), L^a_2(x), L^a_3(x) using explicit representations.
pub fn laguerre_1_e(a: f64, x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_laguerre_1_e(a, x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This function evaluates the generalized Laguerre polynomials L^a_1(x), L^a_2(x), L^a_3(x) using explicit representations.
pub fn laguerre_2_e(a: f64, x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_laguerre_2_e(a, x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This function evaluates the generalized Laguerre polynomials L^a_1(x), L^a_2(x), L^a_3(x) using explicit representations.
pub fn laguerre_3_e(a: f64, x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_laguerre_3_e(a, x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// the generalized Laguerre polynomials L^a_n(x) for a > -1, n >= 0.
pub fn laguerre_n(n: i32, a: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_laguerre_n(n, a, x) }
}

/// the generalized Laguerre polynomials L^a_n(x) for a > -1, n >= 0.
pub fn laguerre_n_e(n: i32, a: f64, x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_laguerre_n_e(n, a, x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}
