//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use crate::{Error, MatrixF64, VectorF64};

#[doc(alias = "gsl_multilarge_linear_L_decomp")]
pub fn linear_L_decomp(L: &mut MatrixF64, tau: &mut VectorF64) -> Result<(), Error> {
    let ret =
        unsafe { sys::gsl_multilarge_linear_L_decomp(L.unwrap_unique(), tau.unwrap_unique()) };
    Error::handle(ret, ())
}
