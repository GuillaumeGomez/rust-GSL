//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Hypergeometric functions are described in Abramowitz & Stegun, Chapters 13 and 15.

use enums;
use ffi;
use std::mem::zeroed;

/// This routine computes the hypergeometric function 0F1(c,x).
pub fn hyperg_0F1(c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_0F1(c, x) }
}

/// This routine computes the hypergeometric function 0F1(c,x).
pub fn hyperg_0F1_e(c: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_0F1_e(c, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the confluent hypergeometric function 1F1(m,n,x) = M(m,n,x) for integer parameters m, n.
pub fn hyperg_1F1_int(m: i32, n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_1F1_int(m, n, x) }
}

/// This routine computes the confluent hypergeometric function 1F1(m,n,x) = M(m,n,x) for integer parameters m, n.
pub fn hyperg_1F1_int_e(m: i32, n: i32, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_1F1_int_e(m, n, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the confluent hypergeometric function 1F1(a,b,x) = M(a,b,x) for general parameters a, b.
pub fn hyperg_1F1(a: f64, b: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_1F1(a, b, x) }
}

/// This routine computes the confluent hypergeometric function 1F1(a,b,x) = M(a,b,x) for general parameters a, b.
pub fn hyperg_1F1_e(a: f64, b: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_1F1_e(a, b, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the confluent hypergeometric function U(m,n,x) for integer parameters m, n.
pub fn hyperg_1F1_U_int(m: i32, n: i32, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_U_int(m, n, x) }
}

/// This routine computes the confluent hypergeometric function U(m,n,x) for integer parameters m, n.
pub fn hyperg_1F1_U_int_e(m: i32, n: i32, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_U_int_e(m, n, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the confluent hypergeometric function U(m,n,x) for integer parameters m, n using the
/// [`ResultE10]`(types/result/struct.ResultE10.html) type to return a result with extended range.
pub fn hyperg_1F1_U_int_e10_e(m: i32, n: i32, x: f64) -> (enums::Value, ::types::ResultE10) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result_e10>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_U_int_e10_e(m, n, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::ResultE10 {
            val: result.val,
            err: result.err,
            e10: result.e10,
        },
    )
}

/// This routine computes the confluent hypergeometric function U(a,b,x).
pub fn hyperg_U(a: f64, b: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_U(a, b, x) }
}

/// This routine computes the confluent hypergeometric function U(a,b,x).
pub fn hyperg_U_e(a: f64, b: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_U_e(a, b, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the confluent hypergeometric function U(a,b,x) using the
/// [`ResultE10]`(types/result/struct.ResultE10.html) type to return a result with extended range.
pub fn hyperg_U_e10_e(a: f64, b: f64, x: f64) -> (enums::Value, ::types::ResultE10) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result_e10>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_U_e10_e(a, b, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::ResultE10 {
            val: result.val,
            err: result.err,
            e10: result.e10,
        },
    )
}

/// This routine computes the Gauss hypergeometric function 2F1(a,b,c,x) = F(a,b,c,x) for |x| < 1.
///
/// If the arguments (a,b,c,x) are too close to a singularity then the function can return the error code
/// [`MaxIter`](enums/type.Value.html) when the series approximation converges too slowly.
/// This occurs in the region of x=1, c - a - b = m for integer m.
pub fn hyperg_2F1(a: f64, b: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1(a, b, c, x) }
}

/// This routine computes the Gauss hypergeometric function 2F1(a,b,c,x) = F(a,b,c,x) for |x| < 1.
///
/// If the arguments (a,b,c,x) are too close to a singularity then the function can return the error code
/// [`MaxIter`](enums/type.Value.html) when the series approximation converges too slowly.
/// This occurs in the region of x=1, c - a - b = m for integer m.
pub fn hyperg_2F1_e(a: f64, b: f64, c: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_2F1_e(a, b, c, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) with complex parameters for |x| < 1.
pub fn hyperg_2F1_conj(aR: f64, aI: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1_conj(aR, aI, c, x) }
}

/// This routine computes the Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) with complex parameters for |x| < 1.
pub fn hyperg_2F1_conj_e(aR: f64, aI: f64, c: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_2F1_conj_e(aR, aI, c, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a,b,c,x) / \Gamma(c) for |x| < 1.
pub fn hyperg_2F1_renorm(a: f64, b: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1_renorm(a, b, c, x) }
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a,b,c,x) / \Gamma(c) for |x| < 1.
pub fn hyperg_2F1_renorm_e(a: f64, b: f64, c: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_2F1_renorm_e(a, b, c, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c) for |x| < 1.
pub fn hyperg_2F1_conj_renorm(aR: f64, aI: f64, c: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F1_conj_renorm(aR, aI, c, x) }
}

/// This routine computes the renormalized Gauss hypergeometric function 2F1(a_R + i a_I, a_R - i a_I, c, x) / \Gamma(c) for |x| < 1.
pub fn hyperg_2F1_conj_renorm_e(
    aR: f64,
    aI: f64,
    c: f64,
    x: f64,
) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_2F1_conj_renorm_e(aR, aI, c, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the hypergeometric function 2F0(a,b,x). The series representation is a divergent hypergeometric series.
/// However, for x < 0 we have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
pub fn hyperg_2F0(a: f64, b: f64, x: f64) -> f64 {
    unsafe { sys::gsl_sf_hyperg_2F0(a, b, x) }
}

/// This routine computes the hypergeometric function 2F0(a,b,x). The series representation is a divergent hypergeometric series.
/// However, for x < 0 we have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
pub fn hyperg_2F0_e(a: f64, b: f64, x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_hyperg_2F0_e(a, b, x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}
