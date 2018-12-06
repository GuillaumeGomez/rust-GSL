//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The following functions are equivalent to the function gsl_pow_int (see [1Small integer powers1](Pow.html)) with an error estimate.

use std::mem::zeroed;
use enums;
use ffi;

/// This routine computes the power x^n for integer n. The power is computed using the minimum number of multiplications.
/// For example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
/// For reasons of efficiency, these functions do not check for overflow or underflow conditions.
/// 
/// ```Rust
/// use rgsl::power;
/// /* compute 3.0**12 */
/// double y = gsl_sf_pow_int(3.0, 12);
/// ```
pub fn pow_int(x: f64, n: i32) -> f64 {
    unsafe { ffi::gsl_sf_pow_int(x, n) }
}

/// This routine computes the power x^n for integer n. The power is computed using the minimum number of multiplications.
/// For example, x^8 is computed as ((x^2)^2)^2, requiring only 3 multiplications.
/// For reasons of efficiency, these functions do not check for overflow or underflow conditions.
/// 
/// ```Rust
/// use rgsl::power;
/// /* compute 3.0**12 */
/// double y = gsl_sf_pow_int(3.0, 12);
/// ```
pub fn pow_int_e(x: f64, n: i32) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_pow_int_e(x, n, &mut result) };

    (enums::Value::from(ret), ::types::Result{val: result.val, err: result.err})
}
