//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Clausen function is defined by the following integral,

Cl_2(x) = - \int_0^x dt \log(2 \sin(t/2))

It is related to the dilogarithm by Cl_2(\theta) = \Im Li_2(\exp(i\theta)).
!*/

use crate::Value;
use std::mem::MaybeUninit;

/// This routine computes the Clausen integral Cl_2(x).
pub fn clausen(x: f64) -> f64 {
    unsafe { sys::gsl_sf_clausen(x) }
}

/// This routine computes the Clausen integral Cl_2(x).
pub fn clausen_e(x: f64) -> (Value, ::types::Result) {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_clausen_e(x, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}
