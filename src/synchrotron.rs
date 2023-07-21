//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{types, Value};
use std::mem::MaybeUninit;

/// This routine computes the first synchrotron function x \int_x^\infty dt K_{5/3}(t) for x >= 0.
#[doc(alias = "gsl_sf_synchrotron_1")]
pub fn synchrotron_1(x: f64) -> f64 {
    unsafe { sys::gsl_sf_synchrotron_1(x) }
}

/// This routine computes the first synchrotron function x \int_x^\infty dt K_{5/3}(t) for x >= 0.
#[doc(alias = "gsl_sf_synchrotron_1_e")]
pub fn synchrotron_1_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_synchrotron_1_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the second synchrotron function x K_{2/3}(x) for x >= 0.
#[doc(alias = "gsl_sf_synchrotron_2")]
pub fn synchrotron_2(x: f64) -> f64 {
    unsafe { sys::gsl_sf_synchrotron_2(x) }
}

/// This routine computes the second synchrotron function x K_{2/3}(x) for x >= 0.
#[doc(alias = "gsl_sf_synchrotron_2_e")]
pub fn synchrotron_2_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_synchrotron_2_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}
