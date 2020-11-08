//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{MatrixF64, Value, VectorF64};
use ffi::FFI;

pub fn applyW(
    x: &MatrixF64,
    w: &VectorF64,
    y: &VectorF64,
    wx: &mut MatrixF64,
    wy: &mut VectorF64,
) -> Value {
    unsafe {
        Value::from(sys::gsl_multifit_linear_applyW(
            x.unwrap_shared(),
            w.unwrap_shared(),
            y.unwrap_shared(),
            wx.unwrap_unique(),
            wy.unwrap_unique(),
        ))
    }
}

pub fn L_decomp(l: &mut MatrixF64, tau: &mut VectorF64) -> Value {
    unsafe {
        Value::from(sys::gsl_multifit_linear_L_decomp(
            l.unwrap_unique(),
            tau.unwrap_unique(),
        ))
    }
}

pub fn lreg(smin: f64, smax: f64, reg_param: &mut VectorF64) -> Value {
    unsafe {
        Value::from(sys::gsl_multifit_linear_lreg(
            smin,
            smax,
            reg_param.unwrap_unique(),
        ))
    }
}

/// Returns `(Value, idx)`.
pub fn lcorner(rho: &VectorF64, eta: &VectorF64) -> (Value, usize) {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner(rho.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    (Value::from(ret), idx)
}

/// Returns `(Value, idx)`.
pub fn lcorner2(reg_param: &VectorF64, eta: &VectorF64) -> (Value, usize) {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner2(reg_param.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    (Value::from(ret), idx)
}

pub fn Lk(p: usize, k: usize, l: &mut MatrixF64) -> Value {
    unsafe { Value::from(sys::gsl_multifit_linear_Lk(p, k, l.unwrap_unique())) }
}

/// Returns `(Value, y, y_err)`.
pub fn linear_est(x: &VectorF64, c: &VectorF64, cov: &MatrixF64) -> (Value, f64, f64) {
    let mut y = 0.;
    let mut y_err = 0.;
    let ret = unsafe {
        sys::gsl_multifit_linear_est(
            x.unwrap_shared(),
            c.unwrap_shared(),
            cov.unwrap_shared(),
            &mut y,
            &mut y_err,
        )
    };
    (Value::from(ret), y, y_err)
}

pub fn linear_residuals(x: &MatrixF64, y: &VectorF64, c: &VectorF64, r: &mut VectorF64) -> Value {
    unsafe {
        Value::from(sys::gsl_multifit_linear_residuals(
            x.unwrap_shared(),
            y.unwrap_shared(),
            c.unwrap_shared(),
            r.unwrap_unique(),
        ))
    }
}
