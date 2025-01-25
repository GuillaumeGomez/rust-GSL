//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Hypergeometric functions are described in Abramowitz & Stegun, Chapters 13 and 15.

use crate::{types, Error};
use std::mem::MaybeUninit;

/// This routine computes the hypergeometric function 0F1(c,x).
#[doc(alias = "gsl_sf_hyperg_0F1")]
pub fn hyperg_0F1(c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_0F1(c, x) }
}

/// This routine computes the hypergeometric function 0F1(c,x).
#[doc(alias = "gsl_sf_hyperg_0F1_e")]
pub fn hyperg_0F1_e(c: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_0F1_e(c, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the confluent hypergeometric function 1F1(m,n,x) = M(m,n,x) for integer parameters m, n.
#[doc(alias = "gsl_sf_hyperg_1F1_int")]
pub fn hyperg_1F1_int(m: i32, n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_1F1_int(m, n, x) }
}

/// This routine computes the confluent hypergeometric function 1F1(m,n,x) = M(m,n,x) for integer parameters m, n.
#[doc(alias = "gsl_sf_hyperg_1F1_int_e")]
pub fn hyperg_1F1_int_e(m: i32, n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_1F1_int_e(m, n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the confluent hypergeometric function 1F1(a,b,x) = M(a,b,x) for general parameters a, b.
#[doc(alias = "gsl_sf_hyperg_1F1")]
pub fn hyperg_1F1(a: f64, b: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_1F1(a, b, x) }
}

/// This routine computes the confluent hypergeometric function 1F1(a,b,x) = M(a,b,x) for general parameters a, b.
#[doc(alias = "gsl_sf_hyperg_1F1_e")]
pub fn hyperg_1F1_e(a: f64, b: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_1F1_e(a, b, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the confluent hypergeometric function U(m,n,x) for integer parameters m, n.
#[doc(alias = "gsl_sf_hyperg_U_int")]
pub fn hyperg_U_int(m: i32, n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_U_int(m, n, x) }
}

/// This routine computes the confluent hypergeometric function U(m,n,x) for integer parameters m, n.
#[doc(alias = "gsl_sf_hyperg_U_int_e")]
pub fn hyperg_U_int_e(m: i32, n: i32, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_U_int_e(m, n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the confluent hypergeometric function U(m,n,x) for integer parameters m, n using the
/// [`ResultE10]`(types/result/struct.ResultE10.html) type to return a result with extended range.
#[doc(alias = "gsl_sf_hyperg_U_int_e10_e")]
pub fn hyperg_U_int_e10_e(m: i32, n: i32, x: f64) -> Result<types::ResultE10, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result_e10>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_U_int_e10_e(m, n, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the confluent hypergeometric function U(a,b,x).
#[doc(alias = "gsl_sf_hyperg_U")]
pub fn hyperg_U(a: f64, b: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_U(a, b, x) }
}

/// This routine computes the confluent hypergeometric function U(a,b,x).
#[doc(alias = "gsl_sf_hyperg_U_e")]
pub fn hyperg_U_e(a: f64, b: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_U_e(a, b, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the confluent hypergeometric function U(a,b,x) using the
/// [`ResultE10]`(types/result/struct.ResultE10.html) type to return a result with extended range.
#[doc(alias = "gsl_sf_hyperg_U_e10_e")]
pub fn hyperg_U_e10_e(a: f64, b: f64, x: f64) -> Result<types::ResultE10, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result_e10>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_U_e10_e(a, b, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Gauss hypergeometric function 2F1(a,b,c,x) = F(a,b,c,x) for |x| < 1.
///
/// If the arguments (a,b,c,x) are too close to a singularity then the function can return the error code
/// [`MaxIter`](enums/type.Value.html) when the series approximation converges too slowly.
/// This occurs in the region of x=1, c - a - b = m for integer m.
#[doc(alias = "gsl_sf_hyperg_2F1")]
pub fn hyperg_2F1(a: f64, b: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1(a, b, c, x) }
}

/// This routine computes the Gauss hypergeometric function 2F1(a,b,c,x) = F(a,b,c,x) for |x| < 1.
///
/// If the arguments (a,b,c,x) are too close to a singularity then the function can return the error code
/// [`MaxIter`](enums/type.Value.html) when the series approximation converges too slowly.
/// This occurs in the region of x=1, c - a - b = m for integer m.
#[doc(alias = "gsl_sf_hyperg_2F1_e")]
pub fn hyperg_2F1_e(a: f64, b: f64, c: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_2F1_e(a, b, c, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) with complex parameters for |x| < 1.
#[doc(alias = "gsl_sf_hyperg_2F1_conj")]
pub fn hyperg_2F1_conj(aR: f64, aI: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1_conj(aR, aI, c, x) }
}

/// This routine computes the Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) with complex parameters for |x| < 1.
#[doc(alias = "gsl_sf_hyperg_2F1_conj_e")]
pub fn hyperg_2F1_conj_e(aR: f64, aI: f64, c: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_2F1_conj_e(aR, aI, c, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a,b,c,x) / \Gamma(c) for |x| < 1.
#[doc(alias = "gsl_sf_hyperg_2F1_renorm")]
pub fn hyperg_2F1_renorm(a: f64, b: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1_renorm(a, b, c, x) }
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a,b,c,x) / \Gamma(c) for |x| < 1.
#[doc(alias = "gsl_sf_hyperg_2F1_renorm_e")]
pub fn hyperg_2F1_renorm_e(a: f64, b: f64, c: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_2F1_renorm_e(a, b, c, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c) for |x| < 1.
#[doc(alias = "gsl_sf_hyperg_2F1_conj_renorm")]
pub fn hyperg_2F1_conj_renorm(aR: f64, aI: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1_conj_renorm(aR, aI, c, x) }
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c) for |x| < 1.
#[doc(alias = "gsl_sf_hyperg_2F1_conj_renorm_e")]
pub fn hyperg_2F1_conj_renorm_e(aR: f64, aI: f64, c: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_2F1_conj_renorm_e(aR, aI, c, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the hypergeometric function 2F0(a,b,x). The series representation is a divergent hypergeometric series.
/// However, for x < 0 we have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
#[doc(alias = "gsl_sf_hyperg_2F0")]
pub fn hyperg_2F0(a: f64, b: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F0(a, b, x) }
}

/// This routine computes the hypergeometric function 2F0(a,b,x). The series representation is a divergent hypergeometric series.
/// However, for x < 0 we have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
#[doc(alias = "gsl_sf_hyperg_2F0_e")]
pub fn hyperg_2F0_e(a: f64, b: f64, x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hyperg_2F0_e(a, b, x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
