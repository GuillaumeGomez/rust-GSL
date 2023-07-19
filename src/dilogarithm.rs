//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{sys, types, Value};
use std::mem::MaybeUninit;

/// These routines compute the dilogarithm for a real argument. In Lewin’s notation this is Li_2(x), the real part of the dilogarithm of a real x.
/// It is defined by the integral representation Li_2(x) = - \Re \int_0^x ds \log(1-s) / s. Note that \Im(Li_2(x)) = 0 for x <= 1, and -\pi\log(x) for x > 1.
///
/// Note that Abramowitz & Stegun refer to the Spence integral S(x)=Li_2(1-x) as the dilogarithm rather than Li_2(x).
#[doc(alias = "gsl_sf_dilog")]
pub fn dilog(x: f64) -> f64 {
    unsafe { sys::gsl_sf_dilog(x) }
}

/// These routines compute the dilogarithm for a real argument. In Lewin’s notation this is Li_2(x), the real part of the dilogarithm of a real x.
/// It is defined by the integral representation Li_2(x) = - \Re \int_0^x ds \log(1-s) / s. Note that \Im(Li_2(x)) = 0 for x <= 1, and -\pi\log(x) for x > 1.
///
/// Note that Abramowitz & Stegun refer to the Spence integral S(x)=Li_2(1-x) as the dilogarithm rather than Li_2(x).
#[doc(alias = "gsl_sf_dilog_e")]
pub fn dilog_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_dilog_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This function computes the full complex-valued dilogarithm for the complex argument z = r \exp(i \theta).
/// The real and imaginary parts of the result are returned in result_re, result_im.
#[doc(alias = "gsl_sf_complex_dilog_e")]
pub fn complex_dilog_e(r: f64, theta: f64) -> Result<(types::Result, types::Result), Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let mut result_im = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe {
        sys::gsl_sf_complex_dilog_e(r, theta, result.as_mut_ptr(), result_im.as_mut_ptr())
    };

    result_handler!(
        ret,
        (
            unsafe { result.assume_init() }.into(),
            unsafe { result_im.assume_init() }.into(),
        )
    )
}
