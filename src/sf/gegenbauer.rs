//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The Gegenbauer polynomials are defined in Abramowitz & Stegun, Chapter 22, where they are known as Ultraspherical polynomials.

use crate::{types, Error};
use std::mem::MaybeUninit;

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
#[doc(alias = "gsl_sf_gegenpoly_1")]
pub fn gegenpoly_1(lambda: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_gegenpoly_1(lambda, x) }
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
#[doc(alias = "gsl_sf_gegenpoly_2")]
pub fn gegenpoly_2(lambda: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_gegenpoly_2(lambda, x) }
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
#[doc(alias = "gsl_sf_gegenpoly_3")]
pub fn gegenpoly_3(lambda: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_gegenpoly_3(lambda, x) }
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
#[doc(alias = "gsl_sf_gegenpoly_1_e")]
pub fn gegenpoly_1_e(lambda: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gegenpoly_1_e(lambda, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
#[doc(alias = "gsl_sf_gegenpoly_2_e")]
pub fn gegenpoly_2_e(lambda: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gegenpoly_2_e(lambda, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This function evaluates the Gegenbauer polynomials C^{(\lambda)}_n(x) using explicit representations for n =1, 2, 3.
#[doc(alias = "gsl_sf_gegenpoly_3_e")]
pub fn gegenpoly_3_e(lambda: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gegenpoly_3_e(lambda, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This function evaluates the Gegenbauer polynomial C^{(\lambda)}_n(x) for a specific value of n, lambda, x subject to \lambda > -1/2, n >= 0.
#[doc(alias = "gsl_sf_gegenpoly_n")]
pub fn gegenpoly_n(n: i32, lambda: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_gegenpoly_n(n, lambda, x) }
}

/// This function evaluates the Gegenbauer polynomial C^{(\lambda)}_n(x) for a specific value of n, lambda, x subject to \lambda > -1/2, n >= 0.
#[doc(alias = "gsl_sf_gegenpoly_n_e")]
pub fn gegenpoly_n_e(n: i32, lambda: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_gegenpoly_n_e(n, lambda, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This function computes an array of Gegenbauer polynomials C^{(\lambda)}_n(x) for n = 0, 1, 2, \dots, nmax, subject to \lambda > -1/2, nmax >= 0.
#[doc(alias = "gsl_sf_gegenpoly_array")]
pub fn gegenpoly_array(lambda: f64, x: f64, result_array: &mut [f64]) -> Result<(), Error> {
    let ret = unsafe {
        sys::gsl_sf_gegenpoly_array(
            result_array.len() as _,
            lambda,
            x,
            result_array.as_mut_ptr(),
        )
    };
    Error::handle(ret, ())
}
