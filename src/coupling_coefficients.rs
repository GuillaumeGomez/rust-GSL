//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Wigner 3-j, 6-j and 9-j symbols give the coupling coefficients for combined angular momentum vectors.
Since the arguments of the standard coupling coefficient functions are integer or half-integer, the arguments of the following functions
are, by convention, integers equal to twice the actual spin value.
!*/

use crate::Value;
use std::mem::MaybeUninit;

/// This routine computes the Wigner 3-j coefficient,
///
/// (ja jb jc
///  ma mb mc)
///
/// where the arguments are given in half-integer units, ja = two_ja/2, ma = two_ma/2, etc.
pub fn _3j(two_ja: i32, two_jb: i32, two_jc: i32, two_ma: i32, two_mb: i32, two_mc: i32) -> f64 {
    unsafe { ::sys::gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc) }
}

/// This routine computes the Wigner 3-j coefficient,
///
/// (ja jb jc
///  ma mb mc)
///
/// where the arguments are given in half-integer units, ja = two_ja/2, ma = two_ma/2, etc.
pub fn _3j_e(
    two_ja: i32,
    two_jb: i32,
    two_jc: i32,
    two_ma: i32,
    two_mb: i32,
    two_mc: i32,
) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe {
        ::sys::gsl_sf_coupling_3j_e(
            two_ja,
            two_jb,
            two_jc,
            two_ma,
            two_mb,
            two_mc,
            result.as_mut_ptr(),
        )
    };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the Wigner 6-j coefficient,
///
/// {ja jb jc
/// jd je jf}
///
/// where the arguments are given in half-integer units, ja = two_ja/2, ma = two_ma/2, etc.
pub fn _6j(two_ja: i32, two_jb: i32, two_jc: i32, two_jd: i32, two_je: i32, two_jf: i32) -> f64 {
    unsafe { ::sys::gsl_sf_coupling_6j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf) }
}

/// This routine computes the Wigner 6-j coefficient,
///
/// {ja jb jc
/// jd je jf}
///
/// where the arguments are given in half-integer units, ja = two_ja/2, ma = two_ma/2, etc.
pub fn _6j_e(
    two_ja: i32,
    two_jb: i32,
    two_jc: i32,
    two_jd: i32,
    two_je: i32,
    two_jf: i32,
) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe {
        ::sys::gsl_sf_coupling_6j_e(
            two_ja,
            two_jb,
            two_jc,
            two_jd,
            two_je,
            two_jf,
            result.as_mut_ptr(),
        )
    };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the Wigner 9-j coefficient,
///
/// {ja jb jc
/// jd je jf
/// jg jh ji}
/// where the arguments are given in half-integer units, ja = two_ja/2, ma = two_ma/2, etc.
pub fn _9j(
    two_ja: i32,
    two_jb: i32,
    two_jc: i32,
    two_jd: i32,
    two_je: i32,
    two_jf: i32,
    two_jg: i32,
    two_jh: i32,
    two_ji: i32,
) -> f64 {
    unsafe {
        ::sys::gsl_sf_coupling_9j(
            two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji,
        )
    }
}

/// This routine computes the Wigner 9-j coefficient,
///
/// {ja jb jc
/// jd je jf
/// jg jh ji}
/// where the arguments are given in half-integer units, ja = two_ja/2, ma = two_ma/2, etc.
pub fn _9j_e(
    two_ja: i32,
    two_jb: i32,
    two_jc: i32,
    two_jd: i32,
    two_je: i32,
    two_jf: i32,
    two_jg: i32,
    two_jh: i32,
    two_ji: i32,
) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe {
        ::sys::gsl_sf_coupling_9j_e(
            two_ja,
            two_jb,
            two_jc,
            two_jd,
            two_je,
            two_jf,
            two_jg,
            two_jh,
            two_ji,
            result.as_mut_ptr(),
        )
    };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}
