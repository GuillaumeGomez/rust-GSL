//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Dawson integral is defined by \exp(-x^2) \int_0^x dt \exp(t^2).
A table of Dawson’s integral can be found in Abramowitz & Stegun, Table 7.5.
!*/

use enums;
use std::mem::zeroed;

/// This routine computes the value of Dawson’s integral for x.
pub fn dawson(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_dawson(x) }
}

/// This routine computes the value of Dawson’s integral for x.
pub fn dawson_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
    let ret = unsafe { ::ffi::gsl_sf_dawson_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}
