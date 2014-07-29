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
}

pub struct gsl_sf_result {
    pub val: c_double,
    pub err: c_double
}