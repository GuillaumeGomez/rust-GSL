//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::Value;

/// This function tests the minimizer specific characteristic size (if applicable to the used minimizer) against absolute tolerance `epsabs`.
/// The test returns `crate::Value::Success` if the size is smaller than tolerance, otherwise crate::Value::Continue is returned.

#[doc(alias = "gsl_multimin_test_size")]
pub fn test_size(size: f64, epsabs: f64) -> Value {
    Value::from(unsafe { sys::gsl_multimin_test_size(size, epsabs) })
}
