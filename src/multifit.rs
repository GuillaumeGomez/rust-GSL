//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;

/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
pub fn covar(J: &::MatrixF64, epsrel: f64, covar: &mut ::MatrixF64) -> ::Value {
    unsafe { ffi::gsl_multifit_covar(ffi::FFI::unwrap(J), epsrel, ffi::FFI::unwrap(covar)) }
}

pub fn test_delta(dx: &::VectorF64, x: &::VectorF64, epsabs: f64, epsrel: f64) -> ::Value {
    unsafe { ffi::gsl_multifit_test_delta(ffi::FFI::unwrap(dx),
        ffi::FFI::unwrap(x), epsabs, epsrel) }
}

pub fn gradient(J: &::MatrixF64, f: &::VectorF64, g: &mut ::VectorF64) -> ::Value {
    unsafe { ffi::gsl_multifit_gradient(ffi::FFI::unwrap(J),
        ffi::FFI::unwrap(f), ffi::FFI::unwrap(g)) }
}