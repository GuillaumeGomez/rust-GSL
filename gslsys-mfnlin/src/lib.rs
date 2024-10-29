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

use std::ptr;

use gsl_sys::GSL_SUCCESS;

use gsl_sys::gsl_blas_ddot;
use gsl_sys::gsl_matrix;
use gsl_sys::gsl_matrix_alloc;
use gsl_sys::gsl_multifit_nlinear_free;
use gsl_sys::gsl_multifit_nlinear_fdtype;
use gsl_sys::gsl_multifit_nlinear_iterate;
use gsl_sys::gsl_multifit_nlinear_rcond;
use gsl_sys::gsl_multifit_nlinear_residual;
use gsl_sys::gsl_multifit_nlinear_scale;
use gsl_sys::gsl_multifit_nlinear_solver;
use gsl_sys::gsl_multifit_nlinear_test;
use gsl_sys::gsl_multifit_nlinear_trs;
use gsl_sys::gsl_multifit_nlinear_trust;
use gsl_sys::gsl_multifit_nlinear_type;
use gsl_sys::gsl_multifit_nlinear_workspace;
use gsl_sys::gsl_vector;
use gsl_sys::gsl_vector_get;
use gsl_sys::gsl_vector_set;
use gsl_sys::gsl_vector_view;
use gsl_sys::gsl_vector_view_array;


pub extern crate libc;

#[repr(C)]
#[derive(Debug, Clone)]
struct multifit_nlinear_data {
  func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
  func_dfs: Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
  params_len: usize,
  ts: Vec<f64>,
  ys: Vec<f64>,
  args: Vec<f64>
}

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

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct gsl_multifit_nlinear_fdf {
  f: fn(*const gsl_vector, *mut multifit_nlinear_data, *mut gsl_vector) -> i32,
  df: *const fn(*const gsl_vector, *mut multifit_nlinear_data, *mut gsl_matrix) -> i32,
  fvv: *const fn(*const gsl_vector, *const gsl_vector, *mut multifit_nlinear_data, *mut gsl_vector) -> i32,
  n: usize,
  p: usize,
  params: *mut multifit_nlinear_data,
  nevalf: usize,
  nevaldf: usize,
  nevalfvv: usize
}

extern "C" {
    pub fn gsl_multifit_nlinear_default_parameters() -> self::gsl_multifit_nlinear_parameters;
}
extern "C" {
    pub fn gsl_multifit_nlinear_alloc(
        T: *const gsl_multifit_nlinear_type,
        params: *const self::gsl_multifit_nlinear_parameters,
        n: usize,
        p: usize,
    ) -> *mut gsl_multifit_nlinear_workspace;
}

extern "C" {
    pub fn gsl_multifit_nlinear_init(
        x: *const gsl_vector,
        fdf: *mut self::gsl_multifit_nlinear_fdf,
        w: *mut gsl_multifit_nlinear_workspace,
    ) -> ::std::os::raw::c_int;
}

fn call_f(x: *const gsl_vector, data: *mut multifit_nlinear_data, f: *mut gsl_vector) -> i32 {

    let n = unsafe { (*data).ts.len() };

    let ts: &Vec<f64> = unsafe { &(*data).ts };
    let ys: &Vec<f64> = unsafe { &(*data).ys };
    let args: &Vec<f64> = unsafe { &(*data).args };

    let mut params: Vec<f64> = Vec::new();
    let params_len = unsafe { (*data).params_len };

    for i in 0..params_len {
        let param = unsafe { gsl_vector_get(x, i) };
        params.push(param);
    }

    for i in 0..n {
        let func_f = unsafe { (*data).func_f };
        let yi = func_f(params.clone(), ts[i], args.clone());
        unsafe { gsl_vector_set(f, i, yi - ys[i]); }
    }

    GSL_SUCCESS
}

fn call_df(x: *const gsl_vector, params: *mut multifit_nlinear_data, df: *mut gsl_matrix) -> i32 {
    GSL_SUCCESS
}

fn call_fvv(x: *const gsl_vector, v: *const gsl_vector, params: *mut multifit_nlinear_data, fvv: *mut gsl_vector) -> i32 {
    GSL_SUCCESS
}

pub unsafe fn gsl_multifit_nlinear_basic_df(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
    params_in: Vec<f64>,
    ts: &Vec<f64>,
    ys: &Vec<f64>,
    args: &Vec<f64>,
    max_iters: u64
) -> (Vec<f64>, Vec<f64>) {

    if ts.len() != ys.len() {
        eprintln!("Time length does not match Ys length!");
        return (vec![], vec![]);
    }

    let mut params: Vec<f64> = params_in.clone();
    let mut covars: Vec<f64> = Vec::with_capacity(params_in.len());

    covars.resize(params_in.len(), 0.0);

    run_gsl_multifit_nlinear_df(
        func_f,
        func_dfs,
        params.as_mut_ptr(),
        covars.as_mut_ptr(),
        params.len(),
        ts,
        ys,
        args,
        max_iters
    );

    (params, covars)
}

fn run_gsl_multifit_nlinear_df(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
    params: *mut f64,
    covars: *mut f64,
    params_len: usize,
    ts: &Vec<f64>,
    ys: &Vec<f64>,
    args: &Vec<f64>,
    max_iters: u64
) {

    let vars_len = ts.len();

    let trust: *const gsl_multifit_nlinear_type = unsafe { gsl_multifit_nlinear_trust };
    let fdf_params: self::gsl_multifit_nlinear_parameters = unsafe { gsl_multifit_nlinear_default_parameters() };

    let f: *mut gsl_vector;
    let jacobian: *mut gsl_matrix;
    let covar: *mut gsl_matrix;

    let data_struct = multifit_nlinear_data {
        func_f: func_f,
        func_dfs: func_dfs.to_vec(),
        params_len: params_len,
        ts: ts.to_vec(),
        ys: ys.to_vec(),
        args: args.to_vec()
    };
    let data_ptr: *mut multifit_nlinear_data = Box::into_raw(Box::new(data_struct));
    let x = unsafe { gsl_vector_view_array(params, params_len) };
    let x_ptr: *const gsl_vector = Box::into_raw(Box::new(x.vector));
    let dof: f64 = (vars_len - params_len) as f64;

    let mut chisq: f64 = 0.0;
    let mut chisq0: f64 = 0.0;
    let mut status: i32 = 0;
    let mut info: i32 = 0;

    const xtol: f64 = 1e-8;
    const gtol: f64 = 1e-8;
    const ftol: f64 = 1e-8;

    let fdf = self::gsl_multifit_nlinear_fdf {
        f: call_f,
        df: ptr::null(),
        fvv: ptr::null(),
        n: vars_len,
        p: params_len,
        params: data_ptr,
        nevalf: 0,
        nevaldf: 0,
        nevalfvv: 0
    };
    let fdf_ptr: *mut self::gsl_multifit_nlinear_fdf = Box::into_raw(Box::new(fdf));

    /* allocate workspace with default parameters */
    let w: *mut gsl_multifit_nlinear_workspace = unsafe {
        gsl_multifit_nlinear_alloc(
            trust,
            &fdf_params as *const self::gsl_multifit_nlinear_parameters,
            vars_len,
            params_len
        )
    };

    /* initialize solver with starting point */
    unsafe {
        gsl_multifit_nlinear_init (
            x_ptr,
            fdf_ptr,
            w
        );
    }
    /* compute initial cost function */
    let f = unsafe { gsl_multifit_nlinear_residual(w) };
    unsafe {
        gsl_blas_ddot(f, f, &mut chisq0);
    }

    /* solve the system within a maximum of max_iters iterations */
    for i in 0..max_iters {

        let mut rcond: f64 = 0.0;

        unsafe { gsl_multifit_nlinear_iterate(w); }
        unsafe { gsl_multifit_nlinear_test(xtol, gtol, ftol, &mut status, w); }

        if status != 0 {

            if status == 1 {
                println!("Reason for stopping: Small Step Size");
            } else if status == 2 {
                println!("Reason for stopping: Small Gradient");
            }

            break;
        }

        unsafe { gsl_multifit_nlinear_rcond(&mut rcond, w) };
        println!("{}", rcond);
    }

    let res = unsafe { gsl_multifit_nlinear_residual(w) };
    for i in 0..params_len {
        unsafe { params.wrapping_add(i).write(gsl_vector_get(res, i)); }
    }

    unsafe {
        gsl_multifit_nlinear_free(w);
    }
}
