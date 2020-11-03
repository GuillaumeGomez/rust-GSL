//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi::FFI;

/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
pub fn covar(J: &::MatrixF64, epsrel: f64, covar: &mut ::MatrixF64) -> ::Value {
    ::Value::from(unsafe {
        sys::gsl_multifit_covar(J.unwrap_shared(), epsrel, covar.unwrap_unique())
    })
}

pub fn test_delta(dx: &::VectorF64, x: &::VectorF64, epsabs: f64, epsrel: f64) -> ::Value {
    ::Value::from(unsafe {
        sys::gsl_multifit_test_delta(dx.unwrap_shared(), x.unwrap_shared(), epsabs, epsrel)
    })
}

pub fn gradient(J: &::MatrixF64, f: &::VectorF64, g: &mut ::VectorF64) -> ::Value {
    ::Value::from(unsafe {
        sys::gsl_multifit_gradient(J.unwrap_shared(), f.unwrap_shared(), g.unwrap_unique())
    })
}
