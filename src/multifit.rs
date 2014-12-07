//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;

/// Compute the covariance matrix cov = inv (J^T J) by QRP^T decomposition of J
pub fn covar(J: &::MatrixF64, epsrel: f64, covar: &mut ::MatrixF64) -> ::Value {
	unsafe { ffi::gsl_multifit_covar(ffi::FFI::unwrap(J) as *const ffi::gsl_matrix, epsrel, ffi::FFI::unwrap(covar)) }
}