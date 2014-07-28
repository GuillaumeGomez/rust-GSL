/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

pub mod Airy {
    use ffi;
    use types::*;
    use std::mem::zeroed;

    /// These routines compute the Airy function Ai(x) with an accuracy specified by mode.
    pub fn Ai(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Ai(x, mode as u32) }
    }

    pub fn Ai_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Ai_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the Airy function Bi(x) with an accuracy specified by mode.
    pub fn Bi(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Bi(x, mode as u32) }
    }

    pub fn Bi_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Bi_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute a scaled version of the Airy function S_A(x) Ai(x). For x>0 the scaling factor S_A(x) is \exp(+(2/3) x^(3/2)), and is 1 for x<0.
    pub fn Ai_scaled(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Ai_scaled(x, mode as u32) }
    }

    pub fn Ai_scaled_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Ai_scaled_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute a scaled version of the Airy function S_B(x) Bi(x). For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
    pub fn Bi_scaled(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Bi_scaled(x, mode as u32) }
    }

    pub fn Bi_scaled_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Bi_scaled_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }
}