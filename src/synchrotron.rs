//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use std::mem::zeroed;
use enums;

/// This routine computes the first synchrotron function x \int_x^\infty dt K_{5/3}(t) for x >= 0.
pub fn synchrotron_1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_synchrotron_1(x) }
}

/// This routine computes the first synchrotron function x \int_x^\infty dt K_{5/3}(t) for x >= 0.
pub fn synchrotron_1_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_synchrotron_1_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the second synchrotron function x K_{2/3}(x) for x >= 0.
pub fn synchrotron_2(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_synchrotron_2(x) }
}

/// This routine computes the second synchrotron function x K_{2/3}(x) for x >= 0.
pub fn synchrotron_2_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_synchrotron_2_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}