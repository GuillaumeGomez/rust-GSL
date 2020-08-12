use libc::{c_double, c_int, c_uint, c_void, size_t, FILE};

use super::gsl_rng;

extern "C" {
    // PLAIN Monte Carlo
    pub fn gsl_monte_plain_alloc(dim: size_t) -> *mut gsl_monte_plain_state;
    pub fn gsl_monte_plain_init(s: *mut gsl_monte_plain_state) -> c_int;
    pub fn gsl_monte_plain_free(s: *mut gsl_monte_plain_state);
    pub fn gsl_monte_plain_integrate(
        f: *mut c_void,
        xl: *const c_double,
        xu: *const c_double,
        dim: size_t,
        calls: size_t,
        r: *mut gsl_rng,
        s: *mut gsl_monte_plain_state,
        result: *mut c_double,
        abserr: *mut c_double,
    ) -> c_int;
    // MISER
    pub fn gsl_monte_miser_alloc(dim: size_t) -> *mut gsl_monte_miser_state;
    pub fn gsl_monte_miser_init(s: *mut gsl_monte_miser_state) -> c_int;
    pub fn gsl_monte_miser_free(s: *mut gsl_monte_miser_state);
    pub fn gsl_monte_miser_integrate(
        f: *mut c_void,
        xl: *const c_double,
        xu: *const c_double,
        dim: size_t,
        calls: size_t,
        r: *mut gsl_rng,
        s: *mut gsl_monte_miser_state,
        result: *mut c_double,
        abserr: *mut c_double,
    ) -> c_int;
    pub fn gsl_monte_miser_params_get(
        s: *mut gsl_monte_miser_state,
        m: *mut gsl_monte_miser_params,
    );
    pub fn gsl_monte_miser_params_set(
        s: *mut gsl_monte_miser_state,
        m: *const gsl_monte_miser_params,
    );
    // VEGAS
    pub fn gsl_monte_vegas_alloc(dim: size_t) -> *mut gsl_monte_vegas_state;
    pub fn gsl_monte_vegas_init(s: *mut gsl_monte_vegas_state) -> c_int;
    pub fn gsl_monte_vegas_free(s: *mut gsl_monte_vegas_state);
    pub fn gsl_monte_vegas_integrate(
        f: *mut c_void,
        xl: *const c_double,
        xu: *const c_double,
        dim: size_t,
        calls: size_t,
        r: *mut gsl_rng,
        s: *mut gsl_monte_vegas_state,
        result: *mut c_double,
        abserr: *mut c_double,
    ) -> c_int;
    pub fn gsl_monte_vegas_chisq(s: *const gsl_monte_vegas_state) -> c_double;
    pub fn gsl_monte_vegas_runval(
        s: *const gsl_monte_vegas_state,
        result: *mut c_double,
        sigma: *mut c_double,
    );
    pub fn gsl_monte_vegas_params_get(
        s: *const gsl_monte_vegas_state,
        params: *mut gsl_monte_vegas_params,
    );
    pub fn gsl_monte_vegas_params_set(
        s: *mut gsl_monte_vegas_state,
        params: *const gsl_monte_vegas_params,
    );
}

#[derive(Debug, Clone)]
#[repr(C)]
pub struct gsl_monte_miser_params {
    /// This parameter specifies the fraction of the currently available number of function calls which
    /// are allocated to estimating the variance at each recursive step. The default value is 0.1.
    pub estimate_frac: c_double,
    /// This parameter specifies the minimum number of function calls required for each estimate of the
    /// variance. If the number of function calls allocated to the estimate using estimate_frac falls
    /// below min_calls then min_calls are used instead. This ensures that each estimate maintains a
    /// reasonable level of accuracy. The default value of min_calls is 16 * dim.
    pub min_calls: size_t,
    /// This parameter specifies the minimum number of function calls required to proceed with a bisection
    /// step. When a recursive step has fewer calls available than min_calls_per_bisection it performs
    /// a plain Monte Carlo estimate of the current sub-region and terminates its branch of the recursion.
    /// The default value of this parameter is 32 * min_calls.
    pub min_calls_per_bisection: size_t,
    /// This parameter controls how the estimated variances for the two sub-regions of a bisection are
    /// combined when allocating points. With recursive sampling the overall variance should scale better
    /// than 1/N, since the values from the sub-regions will be obtained using a procedure which explicitly
    /// minimizes their variance. To accommodate this behavior the MISER algorithm allows the total variance
    /// to depend on a scaling parameter \alpha,
    ///
    /// \Var(f) = {\sigma_a \over N_a^\alpha} + {\sigma_b \over N_b^\alpha}.
    ///
    /// The authors of the original paper describing MISER recommend the value \alpha = 2 as a good choice,
    /// obtained from numerical experiments, and this is used as the default value in this implementation.
    pub alpha: c_double,
    /// This parameter introduces a random fractional variation of size dither into each bisection, which
    /// can be used to break the symmetry of integrands which are concentrated near the exact center of
    /// the hypercubic integration region. The default value of dither is zero, so no variation is introduced.
    /// If needed, a typical value of dither is 0.1.
    pub dither: c_double,
}

#[derive(Debug)]
#[repr(C)]
pub struct gsl_monte_plain_state {
    pub dim: size_t,
    pub x: *mut c_double,
}

#[derive(Debug)]
#[repr(C)]
pub struct gsl_monte_function {
    pub f: *mut c_void,
    pub dim: size_t,
    pub params: *mut c_void,
}

#[derive(Debug)]
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

#[derive(Debug)]
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
    pub mode: c_int,
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

#[derive(Debug)]
#[repr(C)]
pub struct gsl_monte_vegas_params {
    pub alpha: c_double,
    pub iterations: size_t,
    pub stage: c_int,
    pub mode: c_int,
    pub verbose: c_int,
    pub ostream: *mut FILE,
}
