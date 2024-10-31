//
// FFI binding for the GSL library
//

#![cfg_attr(feature = "dox", feature(doc_cfg))]
#![allow(improper_ctypes)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(clippy::redundant_static_lifetimes)]
#![allow(clippy::upper_case_acronyms)]


pub extern crate libc;

mod auto;

pub use auto::*;

use std::ptr;


macro_rules! result_handler {
    ($ret:ident, $value:expr) => {{
        if $ret == GSL_SUCCESS {
            Ok($value)
        } else {
            Err($ret)
        }
    }};
}


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

fn call_df(x: *const gsl_vector, data: *mut multifit_nlinear_data, jacobian: *mut gsl_matrix) -> i32 {

    let n = unsafe { (*data).ts.len() };

    let ts: &Vec<f64> = unsafe { &(*data).ts };
    let func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64> = unsafe { &(*data).func_dfs };
    let args: &Vec<f64> = unsafe { &(*data).args };

    let mut params: Vec<f64> = Vec::new();
    let params_len = unsafe { (*data).params_len };

    for i in 0..params_len {
        let param = unsafe { gsl_vector_get(x, i) };
        params.push(param);
    }

    for i in 0..n {
        for j in 0..params_len {
            unsafe { gsl_matrix_set(jacobian, i, j, func_dfs[j](params.clone(), ts[i], args.clone())); }
        }   
    }

    GSL_SUCCESS
}

fn call_fvv(x: *const gsl_vector, v: *const gsl_vector, params: *mut multifit_nlinear_data, fvv: *mut gsl_vector) -> i32 {
    GSL_SUCCESS
}

pub unsafe fn gsl_multifit_nlinear_basic(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    params_in: Vec<f64>,
    ts: &Vec<f64>,
    ys: &Vec<f64>,
    args: &Vec<f64>,
    max_iters: u64
) -> Result<(Vec<f64>, Vec<f64>, i32), i32> {

    gsl_multifit_nlinear_basic_df(
        func_f,
        &Vec::new(),
        params_in,
        ts,
        ys,
        args,
        max_iters
    )
}

pub unsafe fn gsl_multifit_nlinear_basic_df(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
    params_in: Vec<f64>,
    ts: &Vec<f64>,
    ys: &Vec<f64>,
    args: &Vec<f64>,
    max_iters: u64
) -> Result<(Vec<f64>, Vec<f64>, i32), i32> {

    let mut params: Vec<f64> = params_in.clone();
    let mut parerr: Vec<f64> = Vec::with_capacity(params_in.len());
    let mut status: i32 = 0;
    let mut err_code: i32 = GSL_FAILURE;

    parerr.resize(params_in.len(), 0.0);

    if ts.len() != ys.len() {
        eprintln!("Time length does not match Ys length!");
    } else if params_in.len() > ts.len() {
        eprintln!("Time length is shorter than parameter length!");
    } else {
        run_gsl_multifit_nlinear_df(
            func_f,
            func_dfs,
            params.as_mut_ptr(),
            parerr.as_mut_ptr(),
            params.len(),
            ts,
            ys,
            args,
            max_iters,
            &mut status
        );
        err_code = GSL_SUCCESS;
    }

    result_handler!(err_code, (params, parerr, status))
}

fn run_gsl_multifit_nlinear_df(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
    params: *mut f64,
    parerr: *mut f64,
    params_len: usize,
    ts: &Vec<f64>,
    ys: &Vec<f64>,
    args: &Vec<f64>,
    max_iters: u64,
    status: &mut i32
) {

    let vars_len = ts.len();

    let trust: *const gsl_multifit_nlinear_type = unsafe { gsl_multifit_nlinear_trust };
    let fdf_params: self::gsl_multifit_nlinear_parameters = unsafe { self::gsl_multifit_nlinear_default_parameters() };

    let f: *mut gsl_vector;
    let jacobian: *mut gsl_matrix;
    let covar: *mut gsl_matrix = unsafe { gsl_matrix_alloc(params_len, params_len) };

    let mut data_struct = multifit_nlinear_data {
        func_f: func_f,
        func_dfs: func_dfs.to_vec(),
        params_len: params_len,
        ts: ts.to_vec(),
        ys: ys.to_vec(),
        args: args.to_vec()
    };
    let data_ptr: *mut multifit_nlinear_data = &mut data_struct as *mut multifit_nlinear_data;
    let x = unsafe { gsl_vector_view_array(params, params_len) };
    let x_ptr: *const gsl_vector = &x.vector as *const gsl_vector;
    let dof: f64 = (vars_len - params_len) as f64;

    let mut chisq: f64 = 0.0;
    let mut chisq0: f64 = 0.0;
    let mut info: i32 = 0;

    const xtol: f64 = 1e-8;
    const gtol: f64 = 1e-8;
    const ftol: f64 = 1e-8;

    let call_df_ptr = if func_dfs.len() == params_len {
        call_df as *const fn(*const gsl_vector, *mut multifit_nlinear_data, *mut gsl_matrix) -> i32
    } else {
        ptr::null()
    };

    let mut fdf = self::gsl_multifit_nlinear_fdf {
        f: call_f,
        df: call_df_ptr,
        fvv: ptr::null(),
        n: vars_len,
        p: params_len,
        params: data_ptr,
        nevalf: 0,
        nevaldf: 0,
        nevalfvv: 0
    };
    let fdf_ptr: *mut self::gsl_multifit_nlinear_fdf = &mut fdf as *mut self::gsl_multifit_nlinear_fdf;

    /* allocate workspace with default parameters */
    let w: *mut gsl_multifit_nlinear_workspace = unsafe {
        self::gsl_multifit_nlinear_alloc(
            trust,
            &fdf_params as *const self::gsl_multifit_nlinear_parameters,
            vars_len,
            params_len
        )
    };

    /* initialize solver with starting point */
    unsafe {
        self::gsl_multifit_nlinear_init (
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
        unsafe { gsl_multifit_nlinear_test(xtol, gtol, ftol, &mut *status, w); }

        if *status != 0 {

            if *status == 1 {
                println!("Reason for stopping: Small Step Size");
            } else if *status == 2 {
                println!("Reason for stopping: Small Gradient");
            }

            break;
        }

        unsafe { gsl_multifit_nlinear_rcond(&mut rcond, w); }
        if rcond.is_nan() && i != 0 {
            println!("Reason for stopping: Invalid Status");
            *status = -1;

            break;
        }
    }

    /* compute covariance of best fit parameters */
    let jacobian = unsafe { gsl_multifit_nlinear_jac(w) };
    unsafe { gsl_multifit_nlinear_covar(jacobian, 0.0, covar); }

    /* compute final cost */
    unsafe { gsl_blas_ddot(f, f, &mut chisq); }
    let chisq_dof = (chisq / dof).sqrt();

    let residuals = unsafe { gsl_multifit_nlinear_residual(w) };
    for i in 0..params_len {
        unsafe { params.wrapping_add(i).write(gsl_vector_get(residuals, i)); }
        unsafe { parerr.wrapping_add(i).write(chisq_dof * gsl_matrix_get(covar, i, i).sqrt()); }
    }

    unsafe { gsl_multifit_nlinear_free(w); }
    unsafe { gsl_matrix_free(covar); }
}
