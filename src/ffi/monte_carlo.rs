use libc::{c_double, c_int, c_uint, c_void, size_t, FILE};

use super::{gsl_rng};
use enums;

extern "C" {
    // PLAIN Monte Carlo
    pub fn gsl_monte_plain_alloc(dim: size_t) -> *mut gsl_monte_plain_state;
    pub fn gsl_monte_plain_init(s: *mut gsl_monte_plain_state) -> enums::Value;
    pub fn gsl_monte_plain_free(s: *mut gsl_monte_plain_state);
    pub fn gsl_monte_plain_integrate(f: *mut c_void,
                                     xl: *const c_double,
                                     xu: *const c_double,
                                     dim: size_t,
                                     calls: size_t,
                                     r: *mut gsl_rng,
                                     s: *mut gsl_monte_plain_state,
                                     result: *mut c_double,
                                     abserr: *mut c_double)
                                     -> enums::Value;
    // MISER
    pub fn gsl_monte_miser_alloc(dim: size_t) -> *mut gsl_monte_miser_state;
    pub fn gsl_monte_miser_init(s: *mut gsl_monte_miser_state) -> enums::Value;
    pub fn gsl_monte_miser_free(s: *mut gsl_monte_miser_state);
    pub fn gsl_monte_miser_integrate(f: *mut c_void,
                                     xl: *const c_double,
                                     xu: *const c_double,
                                     dim: size_t,
                                     calls: size_t,
                                     r: *mut gsl_rng,
                                     s: *mut gsl_monte_miser_state,
                                     result: *mut c_double,
                                     abserr: *mut c_double)
                                     -> enums::Value;
    pub fn gsl_monte_miser_params_get(s: *mut gsl_monte_miser_state, m: *mut ::MiserParams);
    pub fn gsl_monte_miser_params_set(s: *mut gsl_monte_miser_state, m: *const ::MiserParams);
    // VEGAS
    pub fn gsl_monte_vegas_alloc(dim: size_t) -> *mut gsl_monte_vegas_state;
    pub fn gsl_monte_vegas_init(s: *mut gsl_monte_vegas_state) -> enums::Value;
    pub fn gsl_monte_vegas_free(s: *mut gsl_monte_vegas_state);
    pub fn gsl_monte_vegas_integrate(f: *mut c_void,
                                     xl: *const c_double,
                                     xu: *const c_double,
                                     dim: size_t,
                                     calls: size_t,
                                     r: *mut gsl_rng,
                                     s: *mut gsl_monte_vegas_state,
                                     result: *mut c_double,
                                     abserr: *mut c_double)
                                     -> enums::Value;
    pub fn gsl_monte_vegas_chisq(s: *const gsl_monte_vegas_state) -> c_double;
    pub fn gsl_monte_vegas_runval(s: *const gsl_monte_vegas_state,
                                  result: *mut c_double,
                                  sigma: *mut c_double);
    pub fn gsl_monte_vegas_params_get(s: *const gsl_monte_vegas_state,
                                      params: *mut gsl_monte_vegas_params);
    pub fn gsl_monte_vegas_params_set(s: *mut gsl_monte_vegas_state,
                                      params: *const gsl_monte_vegas_params);
}


#[repr(C)]
pub struct gsl_monte_plain_state {
    pub dim: size_t,
    pub x: *mut c_double,
}

#[repr(C)]
pub struct gsl_monte_function {
    pub f: *mut c_void,
    pub dim: size_t,
    pub params: *mut c_void,
}

#[repr(C)]
pub struct gsl_monte_miser_state {
    pub min_calls: size_t,
    pub min_calls_per_bisection: size_t,
    pub dither: c_double,
    pub estimate_frac: c_double,
    pub alpha: c_double,
    pub dim: size_t,
    pub estimate_style: c_int,
    pub depth: c_int,
    pub verbose: c_int,
    pub x: *mut c_double,
    pub xmid: *mut c_double,
    pub sigma_l: *mut c_double,
    pub sigma_r: *mut c_double,
    pub fmax_l: *mut c_double,
    pub fmax_r: *mut c_double,
    pub fmin_l: *mut c_double,
    pub fmin_r: *mut c_double,
    pub fsum_l: *mut c_double,
    pub fsum_r: *mut c_double,
    pub fsum2_l: *mut c_double,
    pub fsum2_r: *mut c_double,
    pub hits_l: *mut size_t,
    pub hits_r: *mut size_t,
}

#[repr(C)]
pub struct gsl_monte_vegas_state {
    /* grid */
    pub dim: size_t,
    pub bins_max: size_t,
    pub bins: c_uint,
    pub boxes: c_uint, /* these are both counted along the axes */
    pub xi: *mut c_double,
    pub xin: *mut c_double,
    pub delx: *mut c_double,
    pub weight: *mut c_double,
    pub vol: c_double,

    pub x: *mut c_double,
    pub bin: *mut c_int,
    pub box_: *mut c_int,

    /* distribution */
    pub d: *mut c_double,

    /* control variables */
    pub alpha: c_double,
    pub mode: ::VegasMode,
    pub verbose: c_int,
    pub iterations: c_uint,
    pub stage: c_int,

    /* scratch variables preserved between calls to vegas1/2/3  */
    pub jac: c_double,
    pub wtd_int_sum: c_double,
    pub sum_wgts: c_double,
    pub chi_sum: c_double,
    pub chisq: c_double,

    pub result: c_double,
    pub sigma: c_double,

    pub it_start: c_uint,
    pub it_num: c_uint,
    pub samples: c_uint,
    pub calls_per_box: c_uint,

    pub ostream: *mut FILE,
}

#[repr(C)]
pub struct gsl_monte_vegas_params {
    pub alpha: c_double,
    pub iterations: size_t,
    pub stage: c_int,
    pub mode: ::VegasMode,
    pub verbose: c_int,
    pub ostream: *mut FILE,
}
