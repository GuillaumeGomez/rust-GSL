/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

pub mod Airy {
    use ffi;
    use types::*;
    use enums::*;
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

    /// These routines compute the Airy function derivative Ai'(x) with an accuracy specified by mode.
    pub fn Ai_deriv(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Ai_deriv(x, mode as u32) }
    }

    pub fn Ai_deriv_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Ai_deriv_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the Airy function derivative Bi'(x) with an accuracy specified by mode.
    pub fn Bi_deriv(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Bi_deriv(x, mode as u32) }
    }

    pub fn Bi_deriv_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Bi_deriv_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the scaled Airy function derivative S_A(x) Ai'(x). For x>0 the scaling factor S_A(x) is \exp(+(2/3) x^(3/2)), and is 1 for x<0.
    pub fn Ai_deriv_scaled(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Ai_deriv_scaled(x, mode as u32) }
    }

    pub fn Ai_deriv_scaled_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Ai_deriv_scaled_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the scaled Airy function derivative S_B(x) Bi'(x). For x>0 the scaling factor S_B(x) is exp(-(2/3) x^(3/2)), and is 1 for x<0.
    pub fn Bi_deriv_scaled(x: f64, mode: Gsl::Mode) -> f64 {
        unsafe { ffi::gsl_sf_airy_Bi_deriv_scaled(x, mode as u32) }
    }

    pub fn Bi_deriv_scaled_e(x: f64, mode: Gsl::Mode) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_Bi_deriv_scaled_e(x, mode as u32, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the location of the s-th zero of the Airy function Ai(x).
    pub fn zero_Ai(s: u32) -> f64 {
        unsafe { ffi::gsl_sf_airy_zero_Ai(s) }
    }

    pub fn zero_Ai_e(s: u32) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_zero_Ai_e(s, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the location of the s-th zero of the Airy function Bi(x).
    pub fn zero_Bi(s: u32) -> f64 {
        unsafe { ffi::gsl_sf_airy_zero_Bi(s) }
    }

    pub fn zero_Bi_e(s: u32) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_zero_Bi_e(s, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the location of the s-th zero of the Airy function derivative Ai'(x).
    pub fn zero_Ai_deriv(s: u32) -> f64 {
        unsafe { ffi::gsl_sf_airy_zero_Ai_deriv(s) }
    }

    pub fn zero_Ai_deriv_e(s: u32) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_zero_Ai_deriv_e(s, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }

    /// These routines compute the location of the s-th zero of the Airy function derivative Bi'(x).
    pub fn zero_Bi_deriv(s: u32) -> f64 {
        unsafe { ffi::gsl_sf_airy_zero_Bi_deriv(s) }
    }

    pub fn zero_Bi_deriv_e(s: u32) -> (i32, GslResult) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_airy_zero_Bi_deriv_e(s, &mut result) };

        (ret, GslResult{val: result.val, err: result.err})
    }
}