//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::Value;

/// This function tests for the convergence of the interval [x_lower, x_upper] with absolute error epsabs and relative error epsrel. The
/// test returns ::Value::Success if the following condition is achieved,
///
/// ```text
/// |a - b| < epsabs + epsrel min(|a|,|b|)
/// ```
///
/// when the interval x = [a,b] does not include the origin. If the interval includes the origin then \min(|a|,|b|) is replaced by zero (
/// which is the minimum value of |x| over the interval). This ensures that the relative error is accurately estimated for minima close to
/// the origin.
///
/// This condition on the interval also implies that any estimate of the minimum x_m in the interval satisfies the same condition with
/// respect to the true minimum x_m^*,
///
/// ```text
/// |x_m - x_m^*| < epsabs + epsrel x_m^*
/// ```
///
/// assuming that the true minimum x_m^* is contained within the interval.
pub fn test_interval(x_lower: f64, x_upper: f64, epsabs: f64, epsrel: f64) -> Value {
    Value::from(unsafe { sys::gsl_min_test_interval(x_lower, x_upper, epsabs, epsrel) })
}
