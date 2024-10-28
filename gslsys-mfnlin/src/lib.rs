//
// FFI binding for the GSL library
//

/*
#![cfg_attr(feature = "dox", feature(doc_cfg))]
#![allow(improper_ctypes)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::redundant_static_lifetimes)]
#![allow(clippy::upper_case_acronyms)]
*/

use gsl_sys::gsl_multifit_nlinear_fdtype;
use gsl_sys::gsl_multifit_nlinear_scale;
use gsl_sys::gsl_multifit_nlinear_solver;
use gsl_sys::gsl_multifit_nlinear_trs;
use gsl_sys::gsl_multifit_nlinear_type;
use gsl_sys::gsl_multifit_nlinear_workspace;

pub extern crate libc;

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct gsl_multifit_nlinear_parameters {
    trs: *const gsl_multifit_nlinear_trs,
    scale: *const gsl_multifit_nlinear_scale,
    solver: *const gsl_multifit_nlinear_solver,
    fdtype: gsl_multifit_nlinear_fdtype,
    factor_up: f64,
    factor_down: f64,
    avmax: f64,
    h_df: f64,
    h_fvv: f64
}

extern "C" {
    pub fn gsl_multifit_nlinear_default_parameters() -> gsl_multifit_nlinear_parameters;
}
extern "C" {
    pub fn gsl_multifit_nlinear_alloc(
        T: *const gsl_multifit_nlinear_type,
        params: *const gsl_multifit_nlinear_parameters,
        n: usize,
        p: usize,
    ) -> *mut gsl_multifit_nlinear_workspace;
}
extern "C" {
    pub fn gsl_multifit_nlinear_free(w: *mut gsl_multifit_nlinear_workspace);
}
