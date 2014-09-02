//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::intrinsics::fabsf64;

pub fn subinterval_too_small(a1: f64, a2: f64, b2: f64) -> bool {
    let e = ::DBL_EPSILON;
    let u = ::DBL_MIN;

    let tmp = unsafe { (1f64 + 100f64 * e) * (fabsf64(a2) + 1000f64 * u) };

    unsafe { fabsf64(a1) <= tmp && fabsf64(b2) <= tmp }
}