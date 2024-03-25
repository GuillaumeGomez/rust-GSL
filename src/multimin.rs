//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use crate::Value;

/// This function tests the minimizer specific characteristic size (if applicable to the used minimizer) against absolute tolerance `epsabs`.
/// The test returns `crate::Value::Success` if the size is smaller than tolerance, otherwise crate::Value::Continue is returned.
#[doc(alias = "gsl_multimin_test_size")]
pub fn test_size(size: f64, epsabs: f64) -> Value {
    Value::from(unsafe { sys::gsl_multimin_test_size(size, epsabs) })
}

/// This function tests the norm of the gradient `g` against the absolute tolerance `epsabs`.
/// The gradient of a multidimensional function goes to zero at a minimum. The test returns `crate::Value::Success` if the following condition is achieved, |g| < epsabs
/// and returns `crate::Value::Continue` otherwise. A suitable choice of `epsabs` can be made from the desired accuracy in the function for small variations in `x`.
/// The relationship between these quantities is given by \delta{f} = g\,\delta{x}.
#[doc(alias = "gsl_multimin_test_gradient")]
pub fn test_gradient(g: &crate::VectorF64, epsabs: f64) -> Value {
    Value::from(unsafe { sys::gsl_multimin_test_gradient(g.unwrap_shared(), epsabs) })
}
