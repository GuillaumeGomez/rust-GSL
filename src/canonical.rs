//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Conical Functions P^\mu_{-(1/2)+i\lambda}(x) and Q^\mu_{-(1/2)+i\lambda} are described in Abramowitz & Stegun, Section 8.12.
!*/

use ffi;
use gsl;
use enums;
use std::mem::zeroed;

/// This routine computes the irregular Spherical Conical Function P^{1/2}_{-1/2 + i \lambda}(x) for x > -1.
pub fn half(lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_conicalP_half(lambda, x) }
}

/// This routine computes the irregular Spherical Conical Function P^{1/2}_{-1/2 + i \lambda}(x) for x > -1.
pub fn half_e(lambda: f64, x: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_conicalP_half_e(lambda, x, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}

/// This routine computes the regular Spherical Conical Function P^{-1/2}_{-1/2 + i \lambda}(x) for x > -1.
pub fn mhalf(lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_conicalP_mhalf(lambda, x) }
}

/// This routine computes the regular Spherical Conical Function P^{-1/2}_{-1/2 + i \lambda}(x) for x > -1.
pub fn mhalf_e(lambda: f64, x: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_conicalP_mhalf_e(lambda, x, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}

/// This routine computes the conical function P^0_{-1/2 + i \lambda}(x) for x > -1.
pub fn _0(lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_conicalP_0(lambda, x) }
}

/// This routine computes the conical function P^0_{-1/2 + i \lambda}(x) for x > -1.
pub fn _0_e(lambda: f64, x: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_conicalP_0_e(lambda, x, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}

/// This routine computes the conical function P^1_{-1/2 + i \lambda}(x) for x > -1.
pub fn _1(lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_conicalP_1(lambda, x) }
}

/// This routine computes the conical function P^1_{-1/2 + i \lambda}(x) for x > -1.
pub fn _1_e(lambda: f64, x: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_conicalP_1_e(lambda, x, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}

/// This routine computes the Regular Spherical Conical Function P^{-1/2-l}_{-1/2 + i \lambda}(x) for x > -1, l >= -1.
pub fn sph_reg(l: i32, lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_conicalP_sph_reg(l, lambda, x) }
}

/// This routine computes the Regular Spherical Conical Function P^{-1/2-l}_{-1/2 + i \lambda}(x) for x > -1, l >= -1.
pub fn sph_reg_e(l: i32, lambda: f64, x: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_conicalP_sph_reg_e(l, lambda, x, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}

/// This routine computes the Regular Cylindrical Conical Function P^{-m}_{-1/2 + i \lambda}(x) for x > -1, m >= -1.
pub fn cyl_reg(m: i32, lambda: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_conicalP_cyl_reg(m, lambda, x) }
}

/// This routine computes the Regular Cylindrical Conical Function P^{-m}_{-1/2 + i \lambda}(x) for x > -1, m >= -1.
pub fn cyl_reg_e(m: i32, lambda: f64, x: f64) -> (enums::GslValue, gsl::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_conicalP_cyl_reg_e(m, lambda, x, &mut result) };

    (ret, gsl::Result{val: result.val, err: result.err})
}