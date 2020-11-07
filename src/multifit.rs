//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{MatrixF64, Value, VectorF64};
use ffi::FFI;

/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
pub fn covar(J: &MatrixF64, epsrel: f64, covar: &mut MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_multifit_covar(J.unwrap_shared(), epsrel, covar.unwrap_unique())
    })
}

pub fn test_delta(dx: &VectorF64, x: &VectorF64, epsabs: f64, epsrel: f64) -> Value {
    Value::from(unsafe {
        sys::gsl_multifit_test_delta(dx.unwrap_shared(), x.unwrap_shared(), epsabs, epsrel)
    })
}

pub fn gradient(J: &MatrixF64, f: &VectorF64, g: &mut VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_multifit_gradient(J.unwrap_shared(), f.unwrap_shared(), g.unwrap_unique())
    })
}

pub fn linear_lreg(smin: f64, smax: f64, reg_param: &mut VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_multifit_linear_lreg(smin, smax, reg_param.unwrap_unique()) })
}

/// Returns `(Value, idx)`.
pub fn linear_lcorner(rho: &VectorF64, eta: &VectorF64) -> (Value, usize) {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner(rho.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    (Value::from(ret), idx)
}

/// Returns `(Value, idx)`.
pub fn linear_lcorner2(rho: &VectorF64, eta: &VectorF64) -> (Value, usize) {
    let mut idx = 0;
    let ret = unsafe {
        sys::gsl_multifit_linear_lcorner2(rho.unwrap_shared(), eta.unwrap_shared(), &mut idx)
    };
    (Value::from(ret), idx)
}

pub fn linear_Lk(p: usize, k: usize, L: &mut MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_multifit_linear_Lk(p, k, L.unwrap_unique()) })
}
