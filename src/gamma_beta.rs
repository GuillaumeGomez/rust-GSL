//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
This following routines compute the gamma and beta functions in their full and incomplete forms, as well as various kinds of factorials.
!*/

/// The Gamma function is defined by the following integral,
/// 
/// \Gamma(x) = \int_0^\infty dt  t^{x-1} \exp(-t)
/// 
/// It is related to the factorial function by \Gamma(n)=(n-1)! for positive integer n.
/// Further information on the Gamma function can be found in Abramowitz & Stegun, Chapter 6.
pub mod gamma {
    use ffi;
    use std::mem::zeroed;
    use enums;

    /// These routines compute the Gamma function \Gamma(x), subject to x not being a negative integer or zero. The function is computed using the real Lanczos method.
    /// The maximum value of x such that \Gamma(x) is not considered an overflow is given by the macro GSL_SF_GAMMA_XMAX and is 171.0.
    pub fn gamma(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_gamma(x) }
    }

    /// This routine provides an exponential function \exp(x) using GSL semantics and error checking.
    pub fn gamma_e(x: f64) -> (enums::GslValue, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_gamma_e(x, &mut result) };

        (ret, ::ffi::FFI::wrap(&result))
    }
}