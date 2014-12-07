//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;

/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
pub fn covar(J: &::MatrixF64, epsrel: f64, covar: &mut ::MatrixF64) -> ::Value {
    unsafe { ffi::gsl_multifit_covar(ffi::FFI::unwrap(J) as *const ffi::gsl_matrix, epsrel, ffi::FFI::unwrap(covar)) }
}

pub fn test_delta(dx: &::VectorF64, x: &::VectorF64, epsabs: f64, epsrel: f64) -> ::Value {
    unsafe { ffi::gsl_multifit_test_delta(ffi::FFI::unwrap(dx) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(x) as *const ffi::gsl_vector, epsabs, epsrel) }
}

pub fn gradient(J: &::MatrixF64, f: &::VectorF64, g: &mut ::VectorF64) -> ::Value {
    unsafe { ffi::gsl_multifit_gradient(ffi::FFI::unwrap(J) as *const ffi::gsl_matrix,
        ffi::FFI::unwrap(f) as *const ffi::gsl_vector, ffi::FFI::unwrap(g)) }
}