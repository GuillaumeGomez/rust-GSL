//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
Lambert’s W functions, W(x), are defined to be solutions of the equation W(x) \exp(W(x)) = x. This function has multiple branches for x < 0; however, it has only two real-valued branches.
We define W_0(x) to be the principal branch, where W > -1 for x < 0, and W_{-1}(x) to be the other real branch, where W < -1 for x < 0.
!*/

use enums;
use ffi;
use std::mem::zeroed;

/// This computes the principal branch of the Lambert W function, W_0(x).
pub fn lambert_W0(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_lambert_W0(x) }
}

/// This computes the principal branch of the Lambert W function, W_0(x).
pub fn lambert_W0_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_lambert_W0_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}

/// This computes the secondary real-valued branch of the Lambert W function, W_{-1}(x).
pub fn lambert_Wm1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_lambert_Wm1(x) }
}

/// This computes the secondary real-valued branch of the Lambert W function, W_{-1}(x).
pub fn lambert_Wm1_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_lambert_Wm1_e(x, &mut result) };

    (
        enums::Value::from(ret),
        ::types::Result {
            val: result.val,
            err: result.err,
        },
    )
}
