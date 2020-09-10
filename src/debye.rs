//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Debye functions D_n(x) are defined by the following integral,

D_n(x) = n/x^n \int_0^x dt (t^n/(e^t - 1))

For further information see Abramowitz & Stegun, Section 27.1.
!*/

use enums;
use std::mem::zeroed;

/// This routine computes the first-order Debye function D_1(x) = (1/x) \int_0^x dt (t/(e^t - 1)).
pub fn _1(x: f64) -> f64 {
    unsafe { ::sys::gsl_sf_debye_1(x) }
}

/// This routine computes the first-order Debye function D_1(x) = (1/x) \int_0^x dt (t/(e^t - 1)).
pub fn _1_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_debye_1_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the second-order Debye function D_2(x) = (2/x^2) \int_0^x dt (t^2/(e^t - 1)).
pub fn _2(x: f64) -> f64 {
    unsafe { ::sys::gsl_sf_debye_2(x) }
}

/// This routine computes the second-order Debye function D_2(x) = (2/x^2) \int_0^x dt (t^2/(e^t - 1)).
pub fn _2_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_debye_2_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the third-order Debye function D_3(x) = (3/x^3) \int_0^x dt (t^3/(e^t - 1)).
pub fn _3(x: f64) -> f64 {
    unsafe { ::sys::gsl_sf_debye_3(x) }
}

/// This routine computes the third-order Debye function D_3(x) = (3/x^3) \int_0^x dt (t^3/(e^t - 1)).
pub fn _3_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_debye_3_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the fourth-order Debye function D_4(x) = (4/x^4) \int_0^x dt (t^4/(e^t - 1)).
pub fn _4(x: f64) -> f64 {
    unsafe { ::sys::gsl_sf_debye_4(x) }
}

/// This routine computes the fourth-order Debye function D_4(x) = (4/x^4) \int_0^x dt (t^4/(e^t - 1)).
pub fn _4_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_debye_4_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the fifth-order Debye function D_5(x) = (5/x^5) \int_0^x dt (t^5/(e^t - 1)).
pub fn _5(x: f64) -> f64 {
    unsafe { ::sys::gsl_sf_debye_5(x) }
}

/// This routine computes the fifth-order Debye function D_5(x) = (5/x^5) \int_0^x dt (t^5/(e^t - 1)).
pub fn _5_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_debye_5_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This routine computes the sixth-order Debye function D_6(x) = (6/x^6) \int_0^x dt (t^6/(e^t - 1)).
pub fn _6(x: f64) -> f64 {
    unsafe { ::sys::gsl_sf_debye_6(x) }
}

/// This routine computes the sixth-order Debye function D_6(x) = (6/x^6) \int_0^x dt (t^6/(e^t - 1)).
pub fn _6_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
    let ret = unsafe { ::sys::gsl_sf_debye_6_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}
