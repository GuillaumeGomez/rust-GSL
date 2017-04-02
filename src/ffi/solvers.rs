use libc::{c_double, c_char, c_void, size_t};

use super::{gsl_vector, gsl_matrix};

use enums;

extern "C" {
    pub static gsl_multifit_fdfsolver_lmder: *mut gsl_multifit_fdfsolver_type;
    pub static gsl_multifit_fdfsolver_lmsder: *mut gsl_multifit_fdfsolver_type;

    pub static gsl_root_fsolver_bisection: *mut gsl_root_fsolver_type;
    pub static gsl_root_fsolver_brent: *mut gsl_root_fsolver_type;
    pub static gsl_root_fsolver_falsepos: *mut gsl_root_fsolver_type;

    pub static gsl_root_fdfsolver_newton: *mut gsl_root_fdfsolver_type;
    pub static gsl_root_fdfsolver_secant: *mut gsl_root_fdfsolver_type;
    pub static gsl_root_fdfsolver_steffenson: *mut gsl_root_fdfsolver_type;

    // multifit
    pub fn gsl_multifit_covar(j: *const gsl_matrix,
                              epsrel: c_double,
                              covar: *mut gsl_matrix)
                              -> enums::Value;

    pub fn gsl_multifit_fdfsolver_alloc(T: *const gsl_multifit_fdfsolver_type,
                                        n: size_t,
                                        p: size_t)
                                        -> *mut gsl_multifit_fdfsolver;
    pub fn gsl_multifit_fdfsolver_set(s: *mut gsl_multifit_fdfsolver,
                                      fdf: *mut gsl_multifit_function_fdf,
                                      x: *const gsl_vector)
                                      -> enums::Value;
    pub fn gsl_multifit_fdfsolver_iterate(s: *mut gsl_multifit_fdfsolver) -> enums::Value;
    pub fn gsl_multifit_fdfsolver_free(s: *mut gsl_multifit_fdfsolver);
    pub fn gsl_multifit_fdfsolver_name(s: *const gsl_multifit_fdfsolver) -> *const c_char;
    pub fn gsl_multifit_fdfsolver_position(s: *const gsl_multifit_fdfsolver) -> *mut gsl_vector;

    pub fn gsl_multifit_fsolver_alloc(T: *const gsl_multifit_fsolver_type,
                                      n: size_t,
                                      p: size_t)
                                      -> *mut gsl_multifit_fsolver;
    pub fn gsl_multifit_fsolver_free(s: *mut gsl_multifit_fsolver);
    pub fn gsl_multifit_fsolver_set(s: *mut gsl_multifit_fsolver,
                                    f: *mut ::MultiFitFunction,
                                    x: *const gsl_vector)
                                    -> enums::Value;
    pub fn gsl_multifit_fsolver_iterate(s: *mut gsl_multifit_fsolver) -> enums::Value;
    pub fn gsl_multifit_fsolver_name(s: *const gsl_multifit_fsolver) -> *const c_char;
    pub fn gsl_multifit_fsolver_position(s: *const gsl_multifit_fsolver) -> *mut gsl_vector;

    //pub fn gsl_multifit_robust_weight();
    //pub fn gsl_multifit_robust_residuals();

    pub fn gsl_multifit_gradient(j: *const gsl_matrix,
                                 f: *const gsl_vector,
                                 g: *mut gsl_vector)
                                 -> enums::Value;

    pub fn gsl_multifit_test_delta(dx: *const gsl_vector,
                                   x: *const gsl_vector,
                                   epsabs: c_double,
                                   epsrel: c_double)
                                   -> enums::Value;
    //pub fn gsl_multifit_test_gradient();

    // one-dimensional root
    pub fn gsl_root_fsolver_alloc(T: *const gsl_root_fsolver_type) -> *mut gsl_root_fsolver;
    pub fn gsl_root_fsolver_free(s: *mut gsl_root_fsolver);
    pub fn gsl_root_fsolver_set(s: *mut gsl_root_fsolver,
                                f: *mut gsl_function,
                                x_lower: c_double,
                                x_upper: c_double)
                                -> enums::Value;
    pub fn gsl_root_fsolver_iterate(s: *mut gsl_root_fsolver) -> enums::Value;
    pub fn gsl_root_fsolver_name(s: *const gsl_root_fsolver) -> *const c_char;
    pub fn gsl_root_fsolver_root(s: *const gsl_root_fsolver) -> c_double;
    pub fn gsl_root_fsolver_x_lower(s: *const gsl_root_fsolver) -> c_double;
    pub fn gsl_root_fsolver_x_upper(s: *const gsl_root_fsolver) -> c_double;

    pub fn gsl_root_fdfsolver_alloc(T: *const gsl_root_fdfsolver_type) -> *mut gsl_root_fdfsolver;
    pub fn gsl_root_fdfsolver_free(s: *mut gsl_root_fdfsolver);
    pub fn gsl_root_fdfsolver_set(s: *mut gsl_root_fdfsolver,
                                  fdf: *mut gsl_function_fdf,
                                  root: c_double)
                                  -> enums::Value;
    pub fn gsl_root_fdfsolver_iterate(s: *mut gsl_root_fdfsolver) -> enums::Value;
    pub fn gsl_root_fdfsolver_name(s: *const gsl_root_fdfsolver) -> *const c_char;
    pub fn gsl_root_fdfsolver_root(s: *const gsl_root_fdfsolver) -> c_double;

    pub fn gsl_root_test_interval(x_lower: c_double,
                                  x_upper: c_double,
                                  epsabs: c_double,
                                  epsrel: c_double)
                                  -> enums::Value;
    pub fn gsl_root_test_residual(f: c_double, epsabs: c_double) -> enums::Value;
    pub fn gsl_root_test_delta(x1: c_double,
                               x0: c_double,
                               epsabs: c_double,
                               epsrel: c_double)
                               -> enums::Value;
}

// multifit fsolver/fdfsolver types:

#[repr(C)]
pub struct gsl_multifit_fsolver_type {
    name: *const c_char,
    pub size: size_t,
    pub alloc: Option<extern "C" fn(state: *mut c_void, n: size_t, p: size_t) -> enums::Value>,
    pub set: Option<extern "C" fn(state: *mut c_void,
                                  function: *mut ::MultiFitFunction,
                                  x: *mut gsl_vector,
                                  f: *mut gsl_vector,
                                  dx: *mut gsl_vector)
                                  -> enums::Value>,
    pub iterate: Option<extern "C" fn(state: *mut c_void,
                                      function: *mut ::MultiFitFunction,
                                      x: *mut gsl_vector,
                                      f: *mut gsl_vector,
                                      dx: *mut gsl_vector)
                                      -> enums::Value>,
    pub free: Option<extern "C" fn(state: *mut c_void)>,
}

#[repr(C)]
pub struct gsl_multifit_fsolver {
    pub type_: *const gsl_multifit_fsolver_type,
    pub function: *mut ::MultiFitFunction,
    pub x: *mut gsl_vector,
    pub f: *mut gsl_vector,
    pub dx: *mut gsl_vector,
    pub state: *mut c_void,
}

#[repr(C)]
pub struct gsl_multifit_function_fdf {
    pub f: Option<extern "C" fn(x: *mut gsl_vector, params: *mut c_void, f: *mut gsl_vector)
                                -> enums::Value>,
    pub df: Option<extern "C" fn(x: *mut gsl_vector, params: *mut c_void, df: *mut gsl_matrix)
                                 -> enums::Value>,
    pub fdf: Option<extern "C" fn(x: *mut gsl_vector,
                                  params: *mut c_void,
                                  f: *mut gsl_vector,
                                  df: *mut gsl_matrix)
                                  -> enums::Value>,
    pub n: size_t,
    pub p: size_t,
    pub params: *mut c_void,
}

#[repr(C)]
pub struct gsl_multifit_fdfsolver_type {
    pub name: *const c_char,
    pub size: size_t,
    pub alloc: Option<extern "C" fn(state: *mut c_void, n: size_t, p: size_t) -> enums::Value>,
    pub set: Option<extern "C" fn(state: *mut c_void,
                                  fdf: *mut gsl_multifit_function_fdf,
                                  x: *mut gsl_vector,
                                  f: *mut gsl_vector,
                                  J: *mut gsl_matrix,
                                  dx: *mut gsl_vector)
                                  -> enums::Value>,
    pub iterate: Option<extern "C" fn(state: *mut c_void,
                                      fdf: *mut gsl_multifit_function_fdf,
                                      x: *mut gsl_vector,
                                      f: *mut gsl_vector,
                                      J: *mut gsl_matrix,
                                      dx: *mut gsl_vector)
                                      -> enums::Value>,
    pub free: Option<extern "C" fn(state: *mut c_void)>,
}

#[repr(C)]
pub struct gsl_multifit_fdfsolver {
    pub type_: *const gsl_multifit_fdfsolver_type,
    pub fdf: *mut gsl_multifit_function_fdf,
    pub x: *mut gsl_vector,
    pub f: *mut gsl_vector,
    pub J: *mut gsl_matrix,
    pub dx: *mut gsl_vector,
    pub state: *mut c_void,
}

// one dimensional root solvers:

#[repr(C)]
pub struct gsl_function {
    pub function: Option<extern "C" fn(x: c_double, params: *mut c_void) -> c_double>,
    pub params: *mut c_void,
}

#[repr(C)]
pub struct gsl_root_fsolver_type {
    pub name: c_char,
    pub size: size_t,
    pub set: Option<extern "C" fn(state: *mut c_void,
                                  f: *mut gsl_function,
                                  root: c_double,
                                  x_lower: c_double,
                                  x_upper: c_double)
                                  -> enums::Value>,
    pub iterate: Option<extern "C" fn(state: *mut c_void,
                                      f: *mut gsl_function,
                                      root: c_double,
                                      x_lower: c_double,
                                      x_upper: c_double)
                                      -> enums::Value>,
}

#[repr(C)]
pub struct gsl_root_fsolver {
    pub type_: *const gsl_root_fsolver_type,
    pub gsl_function: *mut gsl_function,
    pub root: c_double,
    pub x_lower: c_double,
    pub x_upper: c_double,
    pub state: *mut c_void,
}

#[repr(C)]
pub struct gsl_function_fdf {
    pub f: Option<extern "C" fn(x: c_double, params: *mut c_void) -> c_double>,
    pub df: Option<extern "C" fn(x: c_double, params: *mut c_void) -> c_double>,
    pub fdf: Option<extern "C" fn(x: c_double,
                                  params: *mut c_void,
                                  y: &mut c_double,
                                  dy: &mut c_double)>,
    pub params: *mut c_void,
}

#[repr(C)]
pub struct gsl_root_fdfsolver_type {
    pub name: c_char,
    pub size: size_t,
    pub set: Option<extern "C" fn(state: *mut c_void, f: *mut gsl_function_fdf, root: c_double)
                                  -> enums::Value>,
    pub iterate: Option<extern "C" fn(state: *mut c_void,
                                      f: *mut gsl_function_fdf,
                                      root: c_double)
                                      -> enums::Value>,
}

#[repr(C)]
pub struct gsl_root_fdfsolver {
    pub type_: *const gsl_root_fdfsolver_type,
    pub gsl_function: *mut gsl_function_fdf,
    pub root: c_double,
    pub state: *mut c_void,
}
