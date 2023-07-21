//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{MatrixF64, Value, VectorF64};
use ffi::FFI;

#[doc(alias = "gsl_multilarge_linear_L_decomp")]
pub fn linear_L_decomp(L: &mut MatrixF64, tau: &mut VectorF64) -> Result<(), Value> {
    let ret =
        unsafe { sys::gsl_multilarge_linear_L_decomp(L.unwrap_unique(), tau.unwrap_unique()) };
    result_handler!(ret, ())
}
