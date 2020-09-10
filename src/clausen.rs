//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Clausen function is defined by the following integral,

Cl_2(x) = - \int_0^x dt \log(2 \sin(t/2))

It is related to the dilogarithm by Cl_2(\theta) = \Im Li_2(\exp(i\theta)).
!*/

use enums;
use ffi;
use std::mem::zeroed;

/// This routine computes the Clausen integral Cl_2(x).
pub fn clausen(x: f64) -> f64 {
    unsafe { sys::gsl_sf_clausen(x) }
}

/// This routine computes the Clausen integral Cl_2(x).
pub fn clausen_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
    let ret = unsafe { sys::gsl_sf_clausen_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}
