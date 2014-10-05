//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::mem::zeroed;
use enums;

/// These routines compute the dilogarithm for a real argument. In Lewin’s notation this is Li_2(x), the real part of the dilogarithm of a real x.
/// It is defined by the integral representation Li_2(x) = - \Re \int_0^x ds \log(1-s) / s. Note that \Im(Li_2(x)) = 0 for x <= 1, and -\pi\log(x) for x > 1.
/// 
/// Note that Abramowitz & Stegun refer to the Spence integral S(x)=Li_2(1-x) as the dilogarithm rather than Li_2(x).
pub fn dilog(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_dilog(x) }
}

/// These routines compute the dilogarithm for a real argument. In Lewin’s notation this is Li_2(x), the real part of the dilogarithm of a real x.
/// It is defined by the integral representation Li_2(x) = - \Re \int_0^x ds \log(1-s) / s. Note that \Im(Li_2(x)) = 0 for x <= 1, and -\pi\log(x) for x > 1.
/// 
/// Note that Abramowitz & Stegun refer to the Spence integral S(x)=Li_2(1-x) as the dilogarithm rather than Li_2(x).
pub fn dilog_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
    let ret = unsafe { ::ffi::gsl_sf_dilog_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This function computes the full complex-valued dilogarithm for the complex argument z = r \exp(i \theta).
/// The real and imaginary parts of the result are returned in result_re, result_im.
pub fn complex_dilog_e(r: f64, theta: f64) -> (enums::value::Value, ::types::Result, ::types::Result) {
    let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
    let mut result_im = unsafe { zeroed::<::ffi::gsl_sf_result>() };
    let ret = unsafe { ::ffi::gsl_sf_complex_dilog_e(r, theta, &mut result, &mut result_im) };

    (ret, ::types::Result{val: result.val, err: result.err}, ::types::Result{val: result_im.val, err: result_im.err})
}