//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The Gegenbauer polynomials are defined in Abramowitz & Stegun, Chapter 22, where they are known as Ultraspherical polynomials.

use ffi;
use std::mem::zeroed;
use enums;

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
pub fn gegenpoly_1(lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_gegenpoly_1(lambda, x) }
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
pub fn gegenpoly_2(lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_gegenpoly_2(lambda, x) }
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
pub fn gegenpoly_3(lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_gegenpoly_3(lambda, x) }
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
pub fn gegenpoly_1_e(lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_gegenpoly_1_e(lambda, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
pub fn gegenpoly_2_e(lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_gegenpoly_2_e(lambda, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
pub fn gegenpoly_3_e(lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_gegenpoly_3_e(lambda, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This function evaluates the Gegenbauer polynomial C^{(\lambda)}_n(x) for a specific value of n, lambda, x subject to \lambda > -1/2, n >= 0.
pub fn gegenpoly_n(n: i32, lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_gegenpoly_n(n, lambda, x) }
}

/// This function evaluates the Gegenbauer polynomial C^{(\lambda)}_n(x) for a specific value of n, lambda, x subject to \lambda > -1/2, n >= 0.
pub fn gegenpoly_n_e(n: i32, lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_gegenpoly_n_e(n, lambda, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This function computes an array of Gegenbauer polynomials C^{(\lambda)}_n(x) for n = 0, 1, 2, \dots, nmax, subject to \lambda > -1/2, nmax >= 0.
pub fn gegenpoly_array(lambda: f64, x: f64, result_array: &mut [f64]) -> enums::Value {
    unsafe { ffi::gsl_sf_gegenpoly_array(result_array.len() as i32, lambda, x, result_array.as_mut_ptr()) }
}
