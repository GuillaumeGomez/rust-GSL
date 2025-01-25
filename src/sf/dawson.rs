//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Dawson integral is defined by \exp(-x^2) \int_0^x dt \exp(t^2).
A table of Dawson’s integral can be found in Abramowitz & Stegun, Table 7.5.
!*/

use crate::{types, Error};
use std::mem::MaybeUninit;

/// This routine computes the value of Dawson’s integral for x.
#[doc(alias = "gsl_sf_dawson")]
pub fn dawson(x: f64) -> f64 {
    unsafe { sys::gsl_sf_dawson(x) }
}

/// This routine computes the value of Dawson’s integral for x.
#[doc(alias = "gsl_sf_dawson_e")]
pub fn dawson_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_dawson_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
