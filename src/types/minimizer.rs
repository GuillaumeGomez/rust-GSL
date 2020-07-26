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
are defined by a gsl_function data type (see [Providing the function to solve](http://www.gnu.org/software/gsl/manual/html_node/Providing-the-function-to-solve.html#Providing-the-function-to-solve)).

##Iteration

The following functions drive the iteration of each algorithm. Each function performs one iteration to update the state of any minimizer of
the corresponding type. The same functions work for all minimizers so that different methods can be substituted at runtime without modifications
to the code.

##Stopping Parameters

A minimization procedure should stop when one of the following conditions is true:

 * A minimum has been found to within the user-specified precision.
 * A user-specified maximum number of iterations has been reached.
 * An error has occurred.

The handling of these conditions is under user control. The function below allows the user to test the precision of the current result.

##Minimization Algorithms

The minimization algorithms described in this section require an initial interval which is guaranteed to contain a minimum—if a and b are the
endpoints of the interval and x is an estimate of the minimum then f(a) > f(x) < f(b). This ensures that the function has at least one minimum
somewhere in the interval. If a valid initial interval is used then these algorithm cannot fail, provided the function is well-behaved.
!*/

use ffi;
use libc::{c_void, free, malloc};

static REL_ERR_VAL: f64 = 1.0e-06f64;
static ABS_ERR_VAL: f64 = 1.0e-10f64;
/* (3 - sqrt(5))/2 */
static GOLDEN_MEAN: f64 = 0.3819660112501052f64;
/* (1 + sqrt(5))/2 */
//static GOLDEN_RATIO   : f64 = 1.6180339887498950f64;

fn safe_func_call<T>(f: ::function<T>, arg: &mut T, x: f64, yp: &mut f64) {
    *yp = f(x, arg);
    if !yp.is_finite() {
        rgsl_error!(
            "computed function value is infinite or NaN",
            ::Value::BadFunction
        );
    }
}

fn compute_f_values<T>(
    f: ::function<T>,
    arg: &mut T,
    x_minimum: f64,
    f_minimum: &mut f64,
    x_lower: f64,
    f_lower: &mut f64,
    x_upper: f64,
    f_upper: &mut f64,
) -> ::Value {
    safe_func_call(f, arg, x_lower, f_lower);
    safe_func_call(f, arg, x_upper, f_upper);
    safe_func_call(f, arg, x_minimum, f_minimum);

    ::Value::Success
}

pub struct Minimizer<T> {
    type_: MinimizerType<T>,
    function: Option<::function<T>>,
    arg: Option<*mut c_void>,
    x_minimum: f64,
    x_lower: f64,
    x_upper: f64,
    f_minimum: f64,
    f_lower: f64,
    f_upper: f64,
    state: *mut c_void,
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
    /// with an error code of ::NoMem.
    pub fn new(t: &MinimizerType<T>) -> Option<Minimizer<T>> {
        let state = unsafe { malloc(t.size) };

        if state.is_null() {
            None
        } else {
            Some(Minimizer {
                type_: t.clone(),
                function: None,
                arg: None,
                x_minimum: 0f64,
                x_lower: 0f64,
                x_upper: 0f64,
                f_minimum: 0f64,
                f_lower: 0f64,
                f_upper: 0f64,
                state: state,
            })
        }
    }

    /// This function sets, or resets, an existing minimizer s to use the function f and the initial search interval [x_lower, x_upper], with
    /// a guess for the location of the minimum x_minimum.
    ///
    /// If the interval given does not contain a minimum, then the function returns an error code of ::Value::Invalid.
    pub fn set(
        &mut self,
        f: ::function<T>,
        arg: &mut T,
        x_minimum: f64,
        x_lower: f64,
        x_upper: f64,
    ) -> ::Value {
        let mut f_minimum = 0f64;
        let mut f_lower = 0f64;
        let mut f_upper = 0f64;

        let status = ::Value::from(compute_f_values(
            f,
            arg,
            x_minimum,
            &mut f_minimum,
            x_lower,
            &mut f_lower,
            x_upper,
            &mut f_upper,
        ));

        if status != ::Value::Success {
            status
        } else {
            self.set_with_values(
                f, arg, x_minimum, f_minimum, x_lower, f_lower, x_upper, f_upper,
            )
        }
    }

    /// This function is equivalent to gsl_min_fminimizer_set but uses the values f_minimum, f_lower and f_upper instead of computing
    /// f(x_minimum), f(x_lower) and f(x_upper).
    pub fn set_with_values(
        &mut self,
        f: ::function<T>,
        arg: &mut T,
        x_minimum: f64,
        f_minimum: f64,
        x_lower: f64,
        f_lower: f64,
        x_upper: f64,
        f_upper: f64,
    ) -> ::Value {
        self.function = Some(f);
        self.arg = unsafe { Some(::std::mem::transmute(arg)) };
        self.x_minimum = x_minimum;
        self.x_lower = x_lower;
        self.x_upper = x_upper;

        if x_lower > x_upper {
            rgsl_error!("invalid interval (lower > upper)", ::Value::Invalid);
        }

        if x_minimum >= x_upper || x_minimum <= x_lower {
            rgsl_error!(
                "x_minimum must lie inside interval (lower < x < upper)",
                ::Value::Invalid
            );
        }

        self.f_upper = f_upper;
        self.f_minimum = f_minimum;
        self.f_lower = f_lower;

        if f_minimum >= f_lower || f_minimum >= f_upper {
            rgsl_error!("endpoints do not enclose a minimum", ::Value::Invalid);
        }

        ::Value::from(unsafe {
            (self.type_.set)(
                self.state,
                self.function.unwrap(),
                ::std::mem::transmute(self.arg.unwrap()),
                x_minimum,
                f_minimum,
                x_lower,
                f_lower,
                x_upper,
                f_upper,
            )
        })
    }

    /// This function returns a pointer to the name of the minimizer. For example,
    ///
    /// ```C
    /// printf("s is a '%s' minimizer\n", gsl_min_fminimizer_name (s));
    /// ```
    ///
    /// would print something like s is a 'brent' minimizer.
    pub fn name(&self) -> String {
        self.type_.name.clone()
    }

    /// This function returns the current estimate of the position of the minimum for the minimizer s.
    pub fn x_minimum(&self) -> f64 {
        self.x_minimum
    }

    /// This function returns the current upper and lower bound of the interval for the minimizer s.
    pub fn x_lower(&self) -> f64 {
        self.x_lower
    }

    /// /// This function returns the current upper and lower bound of the interval for the minimizer s.
    pub fn x_upper(&self) -> f64 {
        self.x_upper
    }

    /// This function returns the value of the function at the current estimate of the minimum and at the upper and lower bounds of the
    /// interval for the minimizer s.
    pub fn f_minimum(&self) -> f64 {
        self.f_minimum
    }

    /// This function returns the value of the function at the current estimate of the minimum and at the upper and lower bounds of the
    /// interval for the minimizer s.
    pub fn f_lower(&self) -> f64 {
        self.f_lower
    }

    /// This function returns the value of the function at the current estimate of the minimum and at the upper and lower bounds of the
    /// interval for the minimizer s.
    pub fn f_upper(&self) -> f64 {
        self.f_upper
    }

    /// This function performs a single iteration of the minimizer s. If the iteration encounters an unexpected problem then an error code
    /// will be returned,
    ///
    /// ::Value::BadFunc
    /// the iteration encountered a singular point where the function evaluated to Inf or NaN.
    ///
    /// ::Value::Failure
    /// the algorithm could not improve the current best approximation or bounding interval.
    ///
    /// The minimizer maintains a current best estimate of the position of the minimum at all times, and the current interval bounding the
    /// minimum. This information can be accessed with the following auxiliary functions,
    pub fn iterate(&mut self) -> ::Value {
        ::Value::from(unsafe {
            (self.type_.iterate)(
                self.state,
                self.function.unwrap(),
                ::std::mem::transmute(self.arg.unwrap()),
                &mut self.x_minimum,
                &mut self.f_minimum,
                &mut self.x_lower,
                &mut self.f_lower,
                &mut self.x_upper,
                &mut self.f_upper,
            )
        })
    }
}

impl<T> Drop for Minimizer<T> {
    fn drop(&mut self) {
        unsafe { free(self.state) };
        self.state = ::std::ptr::null_mut();
    }
}

pub struct MinimizerType<T> {
    pub name: String,
    size: usize,
    set: fn(
        state: *mut c_void,
        f: ::function<T>,
        arg: &mut T,
        x_minimum: f64,
        f_minimum: f64,
        x_lower: f64,
        f_lower: f64,
        x_upper: f64,
        f_upper: f64,
    ) -> ::Value,
    iterate: fn(
        state: *mut c_void,
        f: ::function<T>,
        arg: &mut T,
        x_minimum: &mut f64,
        f_minimum: &mut f64,
        x_lower: &mut f64,
        f_lower: &mut f64,
        x_upper: &mut f64,
        f_upper: &mut f64,
    ) -> ::Value,
}

impl<T> MinimizerType<T> {
    /// The golden section algorithm is the simplest method of bracketing the minimum of a function. It is the slowest algorithm provided
    /// by the library, with linear convergence.
    ///
    /// On each iteration, the algorithm first compares the subintervals from the endpoints to the current minimum. The larger subinterval
    /// is divided in a golden section (using the famous ratio (3-\sqrt 5)/2 = 0.3189660…) and the value of the function at this new point
    /// is calculated. The new value is used with the constraint f(a') > f(x') < f(b') to a select new interval containing the minimum, by
    /// discarding the least useful point. This procedure can be continued indefinitely until the interval is sufficiently small. Choosing
    /// the golden section as the bisection ratio can be shown to provide the fastest convergence for this type of algorithm.
    pub fn golden_section() -> MinimizerType<T> {
        MinimizerType {
            name: "goldensection".to_string(),
            size: ::std::mem::size_of::<goldensection_state_t>() as usize,
            set: goldensection_init,
            iterate: goldensection_iterate,
        }
    }

    /// The Brent minimization algorithm combines a parabolic interpolation with the golden section algorithm. This produces a fast algorithm
    /// which is still robust.
    ///
    /// The outline of the algorithm can be summarized as follows: on each iteration Brent’s method approximates the function using an
    /// interpolating parabola through three existing points. The minimum of the parabola is taken as a guess for the minimum. If it lies
    /// within the bounds of the current interval then the interpolating point is accepted, and used to generate a smaller interval. If the
    /// interpolating point is not accepted then the algorithm falls back to an ordinary golden section step. The full details of Brent’s
    /// method include some additional checks to improve convergence.
    pub fn brent() -> MinimizerType<T> {
        MinimizerType {
            name: "brent".to_string(),
            size: ::std::mem::size_of::<brent_state_t>() as usize,
            set: brent_init,
            iterate: brent_iterate,
        }
    }

    /// This is a variant of Brent’s algorithm which uses the safeguarded step-length algorithm of Gill and Murray.
    pub fn quad_golden() -> MinimizerType<T> {
        MinimizerType {
            name: "quad-golden".to_string(),
            size: ::std::mem::size_of::<quad_golden_state_t>() as usize,
            set: quad_golden_init,
            iterate: quad_golden_iterate,
        }
    }
}

impl<T> Clone for MinimizerType<T> {
    fn clone(&self) -> MinimizerType<T> {
        MinimizerType {
            name: self.name.clone(),
            size: self.size,
            set: self.set,
            iterate: self.iterate,
        }
    }
}

struct goldensection_state_t {
    dummy: f64,
}

#[allow(unused_variables)]
fn goldensection_init<T>(
    vstate: *mut c_void,
    f: ::function<T>,
    arg: &mut T,
    x_minimum: f64,
    f_minimum: f64,
    x_lower: f64,
    f_lower: f64,
    x_upper: f64,
    f_upper: f64,
) -> ::Value {
    let state: &mut goldensection_state_t = unsafe { ::std::mem::transmute(vstate) };

    state.dummy = 0f64;
    ::Value::Success
}

#[allow(unused_variables)]
fn goldensection_iterate<T>(
    vstate: *mut c_void,
    f: ::function<T>,
    arg: &mut T,
    x_minimum: &mut f64,
    f_minimum: &mut f64,
    x_lower: &mut f64,
    f_lower: &mut f64,
    x_upper: &mut f64,
    f_upper: &mut f64,
) -> ::Value {
    let x_center = *x_minimum;
    let x_left = *x_lower;
    let x_right = *x_upper;

    let f_min = *f_minimum;

    /* golden = (3 - sqrt(5))/2 */
    let golden = 0.3819660f64;

    let w_lower = x_center - x_left;
    let w_upper = x_right - x_center;

    let mut f_new = 0f64;

    let x_new = x_center + golden * if w_upper > w_lower { w_upper } else { -w_lower };

    safe_func_call(f, arg, x_new, &mut f_new);

    if f_new < f_min {
        *x_minimum = x_new;
        *f_minimum = f_new;
        ::Value::Success
    } else if x_new < x_center && f_new > f_min {
        *x_lower = x_new;
        *f_lower = f_new;
        ::Value::Success
    } else if x_new > x_center && f_new > f_min {
        *x_upper = x_new;
        *f_upper = f_new;
        ::Value::Success
    } else {
        ::Value::Failure
    }
}

struct brent_state_t {
    d: f64,
    e: f64,
    v: f64,
    w: f64,
    f_v: f64,
    f_w: f64,
}

#[allow(unused_variables)]
fn brent_init<T>(
    vstate: *mut c_void,
    f: ::function<T>,
    arg: &mut T,
    x_minimum: f64,
    f_minimum: f64,
    x_lower: f64,
    f_lower: f64,
    x_upper: f64,
    f_upper: f64,
) -> ::Value {
    let state: &mut brent_state_t = unsafe { ::std::mem::transmute(vstate) };

    /* golden = (3 - sqrt(5))/2 */
    let golden = 0.3819660f64;

    state.v = x_lower + golden * (x_upper - x_lower);
    state.w = state.v;

    let mut f_vw = 0f64;

    state.d = 0f64;
    state.e = 0f64;

    safe_func_call(f, arg, state.v, &mut f_vw);

    state.f_v = f_vw;
    state.f_w = f_vw;

    ::Value::Success
}

#[allow(unused_assignments)]
fn brent_iterate<T>(
    vstate: *mut c_void,
    f: ::function<T>,
    arg: &mut T,
    x_minimum: &mut f64,
    f_minimum: &mut f64,
    x_lower: &mut f64,
    f_lower: &mut f64,
    x_upper: &mut f64,
    f_upper: &mut f64,
) -> ::Value {
    unsafe {
        let state: &mut brent_state_t = ::std::mem::transmute(vstate);

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

        let w_lower = z - x_left;
        let w_upper = x_right - z;

        let tolerance = ::SQRT_DBL_EPSILON * z.abs();

        let mut p = 0f64;
        let mut q = 0f64;
        let mut r = 0f64;

        let midpoint = 0.5f64 * (x_left + x_right);

        if e.abs() > tolerance {
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

        if p.abs() < (0.5f64 * q * r).abs() && p < q * w_lower && p < q * w_upper {
            let t2 = 2f64 * tolerance;

            d = p / q;
            u = z + d;

            if (u - x_left) < t2 || (x_right - u) < t2 {
                d = if z < midpoint { tolerance } else { -tolerance };
            }
        } else {
            e = if z < midpoint {
                x_right - z
            } else {
                -(z - x_left)
            };
            d = golden * e;
        }

        if d.abs() >= tolerance {
            u = z + d;
        } else {
            u = z + if d > 0f64 { tolerance } else { -tolerance };
        }

        state.e = e;
        state.d = d;

        safe_func_call(f, arg, u, &mut f_u);

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

            ::Value::Success
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

                ::Value::Success
            } else if f_u <= f_v || v == z || v == w {
                state.v = u;
                state.f_v = f_u;

                ::Value::Success
            } else {
                ::Value::Success
            }
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
    num_iter: i32,
}

#[allow(unused_variables)]
fn quad_golden_init<T>(
    vstate: *mut c_void,
    f: ::function<T>,
    arg: &mut T,
    x_minimum: f64,
    f_minimum: f64,
    x_lower: f64,
    f_lower: f64,
    x_upper: f64,
    f_upper: f64,
) -> ::Value {
    let state: &mut quad_golden_state_t = unsafe { ::std::mem::transmute(vstate) };

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

    ::Value::Success
}

#[allow(unused_assignments)]
fn quad_golden_iterate<T>(
    vstate: *mut c_void,
    f: ::function<T>,
    arg: &mut T,
    x_minimum: &mut f64,
    f_minimum: &mut f64,
    x_lower: &mut f64,
    f_lower: &mut f64,
    x_upper: &mut f64,
    f_upper: &mut f64,
) -> ::Value {
    unsafe {
        let state: &mut quad_golden_state_t = ::std::mem::transmute(vstate);

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
        let tol = REL_ERR_VAL * x_m.abs() + ABS_ERR_VAL;

        if stored_step.abs() - tol > -2.0f64 * ::DBL_EPSILON {
            /* Fit quadratic */
            let c3 = (x_m - x_small) * (f_m - f_prev_small);
            let mut c2 = (x_m - x_prev_small) * (f_m - f_small);
            let mut c1 = (x_m - x_prev_small) * c2 - (x_m - x_small) * c3;

            c2 = 2f64 * (c2 - c3);

            /* if( c2 != 0 ) */
            if c2.abs() > ::DBL_EPSILON {
                if c2 > 0f64 {
                    c1 = -c1;
                }

                c2 = c2.abs();

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

        let x_trial = x_m + quad_step_size;

        if quad_step_size.abs() < (0.5f64 * prev_stored_step).abs()
            && x_trial > x_l
            && x_trial < x_u
        {
            /* Take quadratic interpolation step */
            step_size = quad_step_size;

            /* Do not evaluate function too close to x_l or x_u */
            if (x_trial - x_l) < 2.0 * tol || (x_u - x_trial) < 2f64 * tol {
                step_size = if x_midpoint >= x_m { 1f64 } else { -1f64 } * tol.abs();
            }

        // This line is supposed to do nothing
        //DEBUG_PRINTF(("quadratic step: %g\n", step_size));
        } else if (x_small != x_prev_small && x_small < x_m && x_prev_small < x_m)
            || (x_small != x_prev_small && x_small > x_m && x_prev_small > x_m)
        {
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

            if inside_interval.abs() <= tol {
                /* Swap inside and outside intervals */
                let tmp = outside_interval;

                outside_interval = inside_interval;
                inside_interval = tmp;
            }

            {
                let step = inside_interval;
                let scale_factor;

                if outside_interval.abs() < inside_interval.abs() {
                    scale_factor = 0.5f64 * (-outside_interval / inside_interval).sqrt();
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
        let x_eval = if step_size.abs() > tol {
            x_m + step_size
        } else {
            x_m + if step_size >= 0f64 { 1f64 } else { -1f64 } * tol.abs()
        };

        /* Evaluate function at the new point x_eval */
        safe_func_call(f, arg, x_eval, &mut f_eval);

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

            if f_eval <= f_small || (x_small - x_m).abs() < 2f64 * ::DBL_EPSILON {
                state.x_prev_small = x_small;
                state.f_prev_small = f_small;

                state.x_small = x_eval;
                state.f_small = f_eval;
            } else if f_eval <= f_prev_small
                || (x_prev_small - x_m).abs() < 2f64 * ::DBL_EPSILON
                || (x_prev_small - x_small).abs() < 2f64 * ::DBL_EPSILON
            {
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

        ::Value::Success
    }
}
