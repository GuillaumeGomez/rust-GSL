//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#One dimensional Minimization

This chapter describes routines for finding minima of arbitrary one-dimensional functions. The library provides low level components for a 
variety of iterative minimizers and convergence tests. These can be combined by the user to achieve the desired solution, with full access 
to the intermediate steps of the algorithms. Each class of methods uses the same framework, so that you can switch between minimizers at 
runtime without needing to recompile your program. Each instance of a minimizer keeps track of its own state, allowing the minimizers to be 
used in multi-threaded programs.

##Overview

The minimization algorithms begin with a bounded region known to contain a minimum. The region is described by a lower bound a and an upper 
bound b, with an estimate of the location of the minimum x.

The value of the function at x must be less than the value of the function at the ends of the interval,

f(a) > f(x) < f(b)
This condition guarantees that a minimum is contained somewhere within the interval. On each iteration a new point x' is selected using one 
of the available algorithms. If the new point is a better estimate of the minimum, i.e. where f(x') < f(x), then the current estimate of 
the minimum x is updated. The new point also allows the size of the bounded interval to be reduced, by choosing the most compact set of 
points which satisfies the constraint f(a) > f(x) < f(b). The interval is reduced until it encloses the true minimum to a desired tolerance. 
This provides a best estimate of the location of the minimum and a rigorous error estimate.

Several bracketing algorithms are available within a single framework. The user provides a high-level driver for the algorithm, and the library 
provides the individual functions necessary for each of the steps. There are three main phases of the iteration. The steps are,

 * initialize minimizer state, s, for algorithm T
 * update s using the iteration T
 * test s for convergence, and repeat iteration if necessary

The state for the minimizers is held in a gsl_min_fminimizer struct. The updating procedure uses only function evaluations (not derivatives).

##Caveats

Note that minimization functions can only search for one minimum at a time. When there are several minima in the search area, the first minimum 
to be found will be returned; however it is difficult to predict which of the minima this will be. In most cases, no error will be reported if 
you try to find a minimum in an area where there is more than one.

With all minimization algorithms it can be difficult to determine the location of the minimum to full numerical precision. The behavior of the 
function in the region of the minimum x^* can be approximated by a Taylor expansion,

y = f(x^*) + (1/2) f''(x^*) (x - x^*)^2

and the second term of this expansion can be lost when added to the first term at finite precision. This magnifies the error in locating x^*, 
making it proportional to \sqrt \epsilon (where \epsilon is the relative accuracy of the floating point numbers). For functions with higher 
order minima, such as x^4, the magnification of the error is correspondingly worse. The best that can be achieved is to converge to the limit 
of numerical accuracy in the function values, rather than the location of the minimum itself.

##Providing the function to minimize

You must provide a continuous function of one variable for the minimizers to operate on. In order to allow for general parameters the functions 
are defined by a gsl_function data type (see Providing the function to solve).

##Iteration

The following functions drive the iteration of each algorithm. Each function performs one iteration to update the state of any minimizer of 
the corresponding type. The same functions work for all minimizers so that different methods can be substituted at runtime without modifications 
to the code.
!*/

use ffi;
use libc::{c_void};
use enums;
use std::intrinsics::{fabsf64, sqrtf64};

static REL_ERR_VAL    : f64 = 1.0e-06f64;
static ABS_ERR_VAL    : f64 = 1.0e-10f64;
/* (3 - sqrt(5))/2 */
static GOLDEN_MEAN    : f64 = 0.3819660112501052f64;
/* (1 + sqrt(5))/2 */
static GOLDEN_RATIO   : f64 = 1.6180339887498950f64;

struct InternMinimizer<T> {
    arg: *mut c_void,
    func: Option<::function<T>>
}

pub struct Minimizer<T> {
    t: *mut ffi::gsl_min_fminimizer,
    intern: InternMinimizer<T>
}

impl<T> Minimizer<T> {
    /// This function returns a pointer to a newly allocated instance of a minimizer of type T. For example, the following code creates an
    /// instance of a golden section minimizer,
    /// 
    /// ```C
    /// const gsl_min_fminimizer_type * T 
    ///   = gsl_min_fminimizer_goldensection;
    /// gsl_min_fminimizer * s 
    ///   = gsl_min_fminimizer_alloc (T);
    /// ```
    /// 
    /// If there is insufficient memory to create the minimizer then the function returns a null pointer and the error handler is invoked
    /// with an error code of enums::NoMem.
    pub fn new<T>(t: &MinimizerType<T>) -> Option<Minimizer<T>> {
        let tmp_pointer = unsafe { ffi::gsl_min_fminimizer_alloc(t.t as *const ffi::gsl_min_fminimizer_type) };

        if tmp_pointer.is_null() {
            None
        } else {
            Some(Minimizer {
                t: tmp_pointer,
                intern: InternMinimizer {
                    arg: ::std::ptr::mut_null(),
                    func: None
                }
            })
        }
    }

    /// This function sets, or resets, an existing minimizer s to use the function f and the initial search interval [x_lower, x_upper], with
    /// a guess for the location of the minimum x_minimum.
    /// 
    /// If the interval given does not contain a minimum, then the function returns an error code of enums::Inval.
    pub fn set<T>(&mut self, f: ::function<T>, arg: &mut T, x_minimum: f64, x_lower: f64, x_upper: f64) -> enums::Value {
        unsafe {
            self.intern.arg = ::std::mem::transmute(arg);
            self.intern.func = Some(::std::mem::transmute(f));

            ffi::gsl_min_fminimizer_set(self.t, &mut self.g_f, x_minimum, x_lower, x_upper)
        }
    }

    /// This function is equivalent to gsl_min_fminimizer_set but uses the values f_minimum, f_lower and f_upper instead of computing
    /// f(x_minimum), f(x_lower) and f(x_upper).
    pub fn set_with_values<T>(&mut self, f: ::function<T>, arg: &mut T, x_minimum: f64, f_minimum: f64, x_lower: f64, f_lower: f64,
        x_upper: f64, f_upper: f64) -> enums::Value {
        unsafe {
            self.intern.arg = ::std::mem::transmute(arg);
            self.intern.func = Some(::std::mem::transmute(f));

            ffi::gsl_min_fminimizer_set_with_values(self.t, &mut self.g_f, x_minimum, f_minimum, x_lower, f_lower, x_upper, f_upper)
        }
    }

    /// This function returns a pointer to the name of the minimizer. For example,
    /// 
    /// ```C
    /// printf("s is a '%s' minimizer\n", gsl_min_fminimizer_name (s));
    /// ```
    /// 
    /// would print something like s is a 'brent' minimizer.
    pub fn name(&self) -> Option<String> {
        let tmp_pointer = unsafe { ffi::gsl_min_fminimizer_name(self.t as *const ffi::gsl_min_fminimizer) };

        if tmp_pointer.is_null() {
            None
        } else {
            unsafe { Some(::std::string::raw::from_buf(tmp_pointer as *const u8)) }
        }
    }

    /// This function performs a single iteration of the minimizer s. If the iteration encounters an unexpected problem then an error code
    /// will be returned,
    /// 
    /// enums::BadFunc
    /// the iteration encountered a singular point where the function evaluated to Inf or NaN.
    /// 
    /// enums::Failure
    /// the algorithm could not improve the current best approximation or bounding interval.
    /// 
    /// The minimizer maintains a current best estimate of the position of the minimum at all times, and the current interval bounding the
    /// minimum. This information can be accessed with the following auxiliary functions,
    pub fn iterate(&self) -> enums::Value {

    }
}

impl<T> Drop for Minimizer<T> {
    fn drop(&mut self) {
        unsafe { ffi::gsl_min_fminimizer_free(self.t) };
        self.t = ::std::ptr::mut_null();
    }
}

pub struct MinimizerType<T> {
    pub name: String,
    size: u64,
    set: fn(state: *mut c_void, f: ::function<T>, arg: &mut T, x_minimum: f64, f_minimum: f64, x_lower: f64, f_lower: f64, x_upper: f64,
        f_upper: f64) -> enums::Value,
    iterate: fn(state: *mut c_void, f: ::function<T>, arg: &mut T, x_minimum: &mut f64, f_minimum: &mut f64, x_lower: &mut f64,
        f_lower: &mut f64, x_upper: &mut f64, f_upper: &mut f64) -> enums::Value
}

impl<T> MinimizerType<T> {
    pub fn golden_section<T>() -> MinimizerType<T> {
        MinimizerType {
            name: "goldensection".to_string(),
            size: ::std::mem::size_of::<goldensection_state_t>() as u64,
            set: goldensection_init,
            iterate: goldensection_iterate
        }
    }

    pub fn brent<T>() -> MinimizerType<T> {
        MinimizerType {
            name: "brent".to_string(),
            size: ::std::mem::size_of::<brent_state_t>() as u64,
            set: brent_init,
            iterate: brent_iterate
        }
    }

    pub fn quad_golden<T>() -> MinimizerType<T> {
        MinimizerType {
            name: "quad-golden",
            size: ::std::mem::size_of::<quad_golden_state_t>() as u64,
            set: quad_golden_init,
            iterate: quad_golden_iterate
        }
    }
}

struct goldensection_state_t {
    dummy: f64
}

fn goldensection_init<T>(vstate: *mut c_void, f: ::function<T>, arg: &mut T, x_minimum: f64, f_minimum: f64, x_lower: f64, f_lower: f64,
    x_upper: f64, f_upper: f64) -> enums::Value {
    let state : &mut goldensection_state_t = ::std::mem::transmute(vstate);

    state.dummy = 0f64;
    enums::Success
}

fn goldensection_iterate<T>(vstate: *mut c_void, f: ::function<f64>, arg: &mut T, x_minimum: &mut f64, f_minimum: &mut f64, x_lower: &mut f64,
    f_lower: &mut f64, x_upper: &mut f64, f_upper: &mut f64) -> enums::Value {
    let x_center = *x_minimum;
    let x_left = *x_lower;
    let x_right = *x_upper;

    let f_min = *f_minimum;

    /* golden = (3 - sqrt(5))/2 */
    let golden = 0.3819660f64;

    let w_lower = (x_center - x_left);
    let w_upper = (x_right - x_center);

    let mut f_new = 0f64;

    let x_new = x_center + golden * if w_upper > w_lower { w_upper } else { -w_lower };

    f(x_new, &mut f_new);

    if f_new < f_min {
        *x_minimum = x_new;
        *f_minimum = f_new;
        enums::Success
    } else if x_new < x_center && f_new > f_min {
        *x_lower = x_new;
        *f_lower = f_new;
        enums::Success
    } else if x_new > x_center && f_new > f_min {
        *x_upper = x_new;
        *f_upper = f_new;
        enums::Success
    } else {
        enums::Failure
    }
}

struct brent_state_t {
    d: f64,
    e: f64,
    v: f64,
    w: f64,
    f_v: f64,
    f_w: f64
}

fn brent_init<T>(vstate: *mut c_void, f: ::function<f64>, arg: &mut T, x_minimum: f64, f_minimum: f64, x_lower: f64, f_lower: f64, x_upper: f64,
    f_upper: f64) -> enums::Value {
    let state : &mut brent_state_t = ::std::mem::transmute(vstate);

    /* golden = (3 - sqrt(5))/2 */
    let golden = 0.3819660f64;

    state.v = x_lower + golden * (x_upper - x_lower);
    state.w = state.v;

    let mut f_vw = 0f64;

    state.d = 0f64;
    state.e = 0f64;

    f(state.v, &mut f_vw);

    state.f_v = f_vw;
    state.f_w = f_vw;

    enums::Success
}

fn brent_iterate<T>(vstate: *mut c_void, f: ::function<f64>, arg: &mut T, x_minimum: &mut f64, f_minimum: &mut f64, x_lower: &mut f64,
    f_lower: &mut f64, x_upper: &mut f64, f_upper: &mut f64) -> enums::Value {
    let state : &mut brent_state_t = ::std::mem::transmute(vstate);

    let x_left = *x_lower;
    let x_right = *x_upper;

    let z = *x_minimum;
    let mut d = state.e;
    let mut e = state.d;
    let mut u = 0f64;
    let mut f_u = 0f64;
    let v = state.v;
    let w = state.w;
    let f_v = state.f_v;
    let f_w = state.f_w;
    let f_z = *f_minimum;

    /* golden = (3 - sqrt(5))/2 */
    let golden = 0.3819660f64;

    let w_lower = (z - x_left);
    let w_upper = (x_right - z);

    let tolerance = ::SQRT_DBL_EPSILON * fabsf64(z);

    let mut p = 0f64;
    let mut q = 0f64;
    let mut r = 0f64;

    let midpoint = 0.5f64 * (x_left + x_right);

    if fabsf64(e) > tolerance {
        /* fit parabola */

        r = (z - w) * (f_z - f_v);
        q = (z - v) * (f_z - f_w);
        p = (z - v) * q - (z - w) * r;
        q = 2f64 * (q - r);

        if q > 0f64 {
            p = -p;
        } else {
            q = -q;
        }

        r = e;
        e = d;
    }

    if fabsf64(p) < fabsf64(0.5f64 * q * r) && p < q * w_lower && p < q * w_upper {
        let t2 = 2f64 * tolerance;

        d = p / q;
        u = z + d;

        if (u - x_left) < t2 || (x_right - u) < t2 {
            d = if z < midpoint { tolerance } else { -tolerance };
        }
    } else {
        e = if z < midpoint { x_right - z } else { -(z - x_left) };
        d = golden * e;
    }


    if fabsf64(d) >= tolerance {
        u = z + d;
    } else {
        u = z + if d > 0f64 { tolerance } else { -tolerance };
    }

    state.e = e;
    state.d = d;

    f(u, &mut f_u);

    if f_u <= f_z {
        if u < z {
            *x_upper = z;
            *f_upper = f_z;
        } else {
            *x_lower = z;
            *f_lower = f_z;
        }

        state.v = w;
        state.f_v = f_w;
        state.w = z;
        state.f_w = f_z;
        *x_minimum = u;
        *f_minimum = f_u;

        enums::Success
    } else {
        if u < z {
            *x_lower = u;
            *f_lower = f_u;
        } else {
            *x_upper = u;
            *f_upper = f_u;
        }

        if f_u <= f_w || w == z {
            state.v = w;
            state.f_v = f_w;
            state.w = u;
            state.f_w = f_u;

            enums::Success
        } else if f_u <= f_v || v == z || v == w {
            state.v = u;
            state.f_v = f_u;

            enums::Success
        } else {
            enums::Success
        }
    }
}

struct quad_golden_state_t {
    step_size: f64,
    stored_step: f64,
    prev_stored_step: f64,
    x_prev_small: f64,
    f_prev_small: f64,
    x_small: f64,
    f_small: f64,
    num_iter: i32
}

fn quad_golden_init<T>(vstate: *mut c_void, f: ::function<T>, arg: &mut T, x_minimum: f64, f_minimum: f64, x_lower: f64, f_lower: f64, x_upper: f64,
    f_upper: f64) -> enums::Value {
    let state : &mut quad_golden_state_t = ::std::mem::transmute(vstate);

    /* For the original behavior, the first value for x_minimum_minimum
     passed in by the user should be a golden section step but we
     don't enforce this here. */

    state.x_prev_small = x_minimum;
    state.x_small = x_minimum;

    state.f_prev_small = f_minimum;
    state.f_small = f_minimum;

    state.step_size = 0f64;
    state.stored_step = 0f64;
    state.prev_stored_step = 0f64;
    state.num_iter = 0;

    enums::Success
}

fn quad_golden_iterate<T>(vstate: *mut c_void, f: ::function<f64>, arg: &mut T, x_minimum: &mut f64, f_minimum: &mut f64, x_lower: &mut f64,
    f_lower: &mut f64, x_upper: &mut f64, f_upper: &mut f64) -> enums::Value {
    let state : &mut quad_golden_state_t = ::std::mem::transmute(vstate);

    let x_m = *x_minimum;
    let f_m = *f_minimum;

    let x_l = *x_lower;
    let x_u = *x_upper;

    let x_small = state.x_small;
    let f_small = state.f_small;

    let x_prev_small = state.x_prev_small;
    let f_prev_small = state.f_prev_small;

    /* update on exit */
    let mut stored_step = state.stored_step;
    /* update on exit */
    let mut prev_stored_step = state.prev_stored_step;
    /* update on exit */
    let mut step_size = state.step_size;

    let mut quad_step_size = prev_stored_step;

    let mut f_eval = 0f64;

    let x_midpoint = 0.5f64 * (x_l + x_u);
    /* total error tolerance */
    let tol = REL_ERR_VAL * fabsf64(x_m) + ABS_ERR_VAL;

    if fabsf64(stored_step) - tol > -2.0f64 * ::DBL_EPSILON {
        /* Fit quadratic */
        let c3 = (x_m - x_small) * (f_m - f_prev_small);
        let c2 = (x_m - x_prev_small) * (f_m - f_small);
        let c1 = (x_m - x_prev_small) * c2 - (x_m - x_small) * c3;

        c2 = 2f64 * (c2 - c3);

        /* if( c2 != 0 ) */
        if fabsf64(c2) > ::DBL_EPSILON {
            if c2 > 0f64 {
                c1 = -c1;
            }

            c2 = fabsf64(c2);

            quad_step_size = c1 / c2;
        } else {
            /* Handle case where c2 ~=~ 0  */
            /* Insure that the line search will NOT take a quadratic
                interpolation step in this iteration */
            quad_step_size = stored_step;
        }

        prev_stored_step = stored_step;
        stored_step = step_size;
    }

    let mut x_trial = x_m + quad_step_size;

    if fabsf64(quad_step_size) < fabsf64(0.5f64 * prev_stored_step) && x_trial > x_l && x_trial < x_u {
        /* Take quadratic interpolation step */
        step_size = quad_step_size;

        /* Do not evaluate function too close to x_l or x_u */
        if (x_trial - x_l) < 2.0 * tol || (x_u - x_trial) < 2f64 * tol {
            step_size = if x_midpoint >= x_m { 1f64 } else { -1f64 } * fabsf64(tol);
        }

        // This line is supposed to do nothing
        //DEBUG_PRINTF(("quadratic step: %g\n", step_size));
    }
    else if (x_small != x_prev_small && x_small < x_m && x_prev_small < x_m) ||
        (x_small != x_prev_small && x_small > x_m && x_prev_small > x_m) {
        /* Take safeguarded function comparison step */
        let mut outside_interval = 0f64;
        let mut inside_interval = 0f64;

        if x_small < x_m {
            outside_interval = x_l - x_m;
            inside_interval = x_u - x_m;
        } else {
            outside_interval = x_u - x_m;
            inside_interval = x_l - x_m;
        }

        if fabsf64(inside_interval) <= tol {
            /* Swap inside and outside intervals */
            let tmp = outside_interval;

            outside_interval = inside_interval;
            inside_interval = tmp;
        }

        {
            let step = inside_interval;
            let scale_factor;

            if fabsf64(outside_interval) < fabsf64(inside_interval) {
                scale_factor = 0.5f64 * sqrtf64(-outside_interval / inside_interval);
            } else {
                scale_factor = (5f64 / 11f64) * (0.1f64 - inside_interval / outside_interval);
            }

            state.stored_step = step;
            step_size = scale_factor * step;
        }

        // This line is supposed to do nothing
        //DEBUG_PRINTF(("safeguard step: %g\n", step_size));
    } else {
        /* Take golden section step */
        let step = if x_m < x_midpoint {
            x_u - x_m
        } else {
            x_l - x_m
        };

        state.stored_step = step;
        step_size = GOLDEN_MEAN * step;

        // This line is supposed to do nothing
        //DEBUG_PRINTF(("golden step: %g\n", step_size));
    }

    /* Do not evaluate function too close to x_minimum */
    let mut x_eval = if fabsf64(step_size) > tol {
        x_m + step_size
    } else {
        x_m + if step_size >= 0f64 { 1f64 } else { -1f64 } * fabsf64(tol)
    };

    /* Evaluate function at the new point x_eval */
    f(x_eval, &mut f_eval);

    /* Update {x,f}_lower, {x,f}_upper, {x,f}_prev_small, {x,f}_small, and {x,f}_minimum */
    if f_eval <= f_m {
        if x_eval < x_m {
            *x_upper = x_m;
            *f_upper = f_m;
        } else {
            *x_lower = x_m;
            *f_upper = f_m;
        }

      state.x_prev_small = x_small;
      state.f_prev_small = f_small;

      state.x_small = x_m;
      state.f_small = f_m;

      *x_minimum = x_eval;
      *f_minimum = f_eval;
    } else {
        if x_eval < x_m {
            *x_lower = x_eval;
            *f_lower = f_eval;
        } else {    
            *x_upper = x_eval;
            *f_upper = f_eval;
        }

        if f_eval <= f_small || fabsf64(x_small - x_m) < 2f64 * ::DBL_EPSILON {
            state.x_prev_small = x_small;
            state.f_prev_small = f_small;

            state.x_small = x_eval;
            state.f_small = f_eval;
        } else if f_eval <= f_prev_small || fabsf64(x_prev_small - x_m) < 2f64 * ::DBL_EPSILON ||
           fabsf64(x_prev_small - x_small) < 2f64 * ::DBL_EPSILON {
            state.x_prev_small = x_eval;
            state.f_prev_small = f_eval;
        }
    }

    /* Update stored values for next iteration */
    state.stored_step = stored_step;
    state.prev_stored_step = prev_stored_step;
    state.step_size = step_size;
    state.num_iter += 1;


    // This line is supposed to do nothing
    //DEBUG_PRINTF(("[%d] Final State: %g  %g  %g\n", state->num_iter, x_l, x_m, x_u));

    enums::Success
}