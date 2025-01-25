//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use crate::Error;

/// This function returns the index i of the array x_array such that `x_array[i] <= x < x_array[i+1]`.
/// The index is searched for in the range `[index_lo,index_hi]`.
#[doc(alias = "gsl_interp_bsearch")]
pub fn bsearch(x_array: &[f64], x: f64, index_lo: usize, index_hi: usize) -> usize {
    unsafe { sys::gsl_interp_bsearch(x_array.as_ptr(), x, index_lo, index_hi) }
}

/// This function returns the interpolated value of y for a given point x, using the interpolation
/// object interp, data arrays xa and ya and the accelerator acc. When x is outside the range of xa,
/// the error code crate::Dom is returned with a value of rgsl::NAN for y.
#[doc(alias = "gsl_interp_eval")]
pub fn eval(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    x: f64,
    acc: &mut crate::InterpAccel,
) -> f64 {
    unsafe {
        sys::gsl_interp_eval(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            x,
            &mut acc.0,
        )
    }
}

/// This function returns the interpolated value of y for a given point x, using the interpolation
/// object interp, data arrays xa and ya and the accelerator acc. When x is outside the range of xa,
/// the error code crate::Dom is returned with a value of rgsl::NAN for y.
///
/// Returns `y`.
#[doc(alias = "gsl_interp_eval_e")]
pub fn eval_e(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    x: f64,
    acc: &mut crate::InterpAccel,
) -> Result<f64, Error> {
    let mut y = 0.;
    let ret = unsafe {
        sys::gsl_interp_eval_e(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            x,
            &mut acc.0,
            &mut y,
        )
    };
    Error::handle(ret, y)
}

/// This function returns the derivative d of an interpolated function for a given point x, using
/// the interpolation object interp, data arrays xa and ya and the accelerator acc.
#[doc(alias = "gsl_interp_eval_deriv")]
pub fn eval_deriv(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    x: f64,
    acc: &mut crate::InterpAccel,
) -> f64 {
    unsafe {
        sys::gsl_interp_eval_deriv(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            x,
            &mut acc.0,
        )
    }
}

/// This function returns the derivative d of an interpolated function for a given point x, using
/// the interpolation object interp, data arrays xa and ya and the accelerator acc.
///
/// Returns `(Value, d)`.
#[doc(alias = "gsl_interp_eval_deriv_e")]
pub fn eval_deriv_e(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    x: f64,
    acc: &mut crate::InterpAccel,
) -> Result<f64, Error> {
    let mut d = 0.;
    let ret = unsafe {
        sys::gsl_interp_eval_deriv_e(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            x,
            &mut acc.0,
            &mut d,
        )
    };
    Error::handle(ret, d)
}

/// This function returns the second derivative d2 of an interpolated function for a given point x,
/// using the interpolation object interp, data arrays xa and ya and the accelerator acc.
#[doc(alias = "gsl_interp_eval_deriv2")]
pub fn eval_deriv2(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    x: f64,
    acc: &mut crate::InterpAccel,
) -> f64 {
    unsafe {
        sys::gsl_interp_eval_deriv2(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            x,
            &mut acc.0,
        )
    }
}

/// This function returns the second derivative d2 of an interpolated function for a given point x,
/// using the interpolation object interp, data arrays xa and ya and the accelerator acc.
///
/// Returns `(Value, d2)`.
#[doc(alias = "gsl_interp_eval_deriv2_e")]
pub fn eval_deriv2_e(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    x: f64,
    acc: &mut crate::InterpAccel,
) -> Result<f64, Error> {
    let mut d2 = 0.;
    let ret = unsafe {
        sys::gsl_interp_eval_deriv2_e(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            x,
            &mut acc.0,
            &mut d2,
        )
    };
    Error::handle(ret, d2)
}

/// This function returns the numerical integral result of an interpolated function over the range
/// [a, b], using the interpolation object interp, data arrays xa and ya and the accelerator acc.
#[doc(alias = "gsl_interp_eval_integ")]
pub fn eval_integ(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    a: f64,
    b: f64,
    acc: &mut crate::InterpAccel,
) -> f64 {
    unsafe {
        sys::gsl_interp_eval_integ(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            a,
            b,
            &mut acc.0,
        )
    }
}

/// This function returns the numerical integral result of an interpolated function over the range
/// [a, b], using the interpolation object interp, data arrays xa and ya and the accelerator acc.
///
/// Returns `(Value, result)`.
#[doc(alias = "gsl_interp_eval_integ_e")]
pub fn eval_integ_e(
    interp: &crate::Interp,
    xa: &[f64],
    ya: &[f64],
    a: f64,
    b: f64,
    acc: &mut crate::InterpAccel,
) -> Result<f64, Error> {
    let mut result = 0.;
    let ret = unsafe {
        sys::gsl_interp_eval_integ_e(
            interp.unwrap_shared(),
            xa.as_ptr(),
            ya.as_ptr(),
            a,
            b,
            &mut acc.0,
            &mut result,
        )
    };
    Error::handle(ret, result)
}
