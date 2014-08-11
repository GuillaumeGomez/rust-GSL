//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use gsl;
use std::mem::zeroed;
use enums;

/// This function multiplies x and y storing the product and its associated error in result.
pub fn multiply_e(x: f64, y: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
    let ret = unsafe { ::ffi::gsl_sf_multiply_e(x, y, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}

/// This function multiplies x and y with associated absolute errors dx and dy.
/// The product xy +/- xy \sqrt((dx/x)^2 +(dy/y)^2) is stored in result.
pub fn multiply_err_e(x: f64, dx: f64, y: f64, dy: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
    let ret = unsafe { ::ffi::gsl_sf_multiply_err_e(x, dx, y, dy, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}