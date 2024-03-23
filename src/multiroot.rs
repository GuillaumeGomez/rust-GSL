//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Multiroot test algorithms, See `rgsl::types::multiroot` for solvers.

use crate::ffi::FFI;
use crate::Value;

#[doc(alias = "gsl_multiroot_test_delta")]
pub fn test_delta(dx: &crate::VectorF64, x: &crate::VectorF64, epsabs: f64, epsrel: f64) -> Value {
    Value::from(unsafe {
        sys::gsl_multiroot_test_delta(dx.unwrap_shared(), x.unwrap_shared(), epsabs, epsrel)
    })
}

#[doc(alias = "gsl_multiroot_test_residual")]
pub fn test_residual(f: &crate::VectorF64, epsabs: f64) -> Value {
    Value::from(unsafe { sys::gsl_multiroot_test_residual(f.unwrap_shared(), epsabs) })
}
