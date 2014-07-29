/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

use libc::{c_double, c_int, c_uint};
use types;

extern "C" {
    // Airy functions
    pub fn gsl_sf_airy_Ai(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Ai_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Ai_deriv(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi_deriv(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Ai_deriv_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi_deriv_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_zero_Ai(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Ai_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_zero_Bi(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Bi_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_zero_Ai_deriv(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Ai_deriv_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_zero_Bi_deriv(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Bi_deriv_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;

    // Bessel functions
    // Regular Modified Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_I0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_I1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_In_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_I0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_I1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_In_scaled_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In_scaled_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Regular Modified Spherical Bessel Functions
    // The regular modified spherical Bessel functions i_l(x) are related to the modified Bessel functions of fractional order, i_l(x) = \sqrt{\pi/(2x)} I_{l+1/2}(x)
    pub fn gsl_sf_bessel_i0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_i1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_i2_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i2_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_il_scaled(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_il_scaled_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_il_scaled_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_Inu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Inu_scaled(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_scaled_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Regular Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_J0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_J1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Jn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Jn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Regular Spherical Bessel Functions
    pub fn gsl_sf_bessel_j0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_j1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_j2(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_jl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_jl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_jl_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_jl_steed_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Regular Bessel Function—Fractional Order
    pub fn gsl_sf_bessel_Jnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jnu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_sequence_Jnu_e(nu: c_double, mode: types::gsl_mode_t, size: i64, v: *mut c_double) -> c_int;
    // Irregular Modified Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_K0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_K1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_K0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_K1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_scaled_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_scaled_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;

}

pub struct gsl_sf_result {
    pub val: c_double,
    pub err: c_double
}