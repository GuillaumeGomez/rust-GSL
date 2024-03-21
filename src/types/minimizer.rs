//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# One dimensional Minimization

This chapter describes routines for finding minima of arbitrary one-dimensional functions. The
library provides low level components for a variety of iterative minimizers and convergence tests.
These can be combined by the user to achieve the desired solution, with full access to the
intermediate steps of the algorithms. Each class of methods uses the same framework, so that you can
switch between minimizers at runtime without needing to recompile your program. Each instance of a
minimizer keeps track of its own state, allowing the minimizers to be used in multi-threaded
programs.

## Overview

The minimization algorithms begin with a bounded region known to contain a minimum. The region is
described by a lower bound a and an upper bound b, with an estimate of the location of the minimum
x.

The value of the function at x must be less than the value of the function at the ends of the
interval,

f(a) > f(x) < f(b)

This condition guarantees that a minimum is contained somewhere within the interval. On each
iteration a new point x' is selected using one of the available algorithms. If the new point is a
better estimate of the minimum, i.e. where f(x') < f(x), then the current estimate of the minimum x
is updated. The new point also allows the size of the bounded interval to be reduced, by choosing
the most compact set of points which satisfies the constraint f(a) > f(x) < f(b). The interval is
reduced until it encloses the true minimum to a desired tolerance. This provides a best estimate of
the location of the minimum and a rigorous error estimate.

Several bracketing algorithms are available within a single framework. The user provides a
high-level driver for the algorithm, and the library provides the individual functions necessary for
each of the steps. There are three main phases of the iteration. The steps are,

 * initialize minimizer state, s, for algorithm T
 * update s using the iteration T
 * test s for convergence, and repeat iteration if necessary

The state for the minimizers is held in a gsl_min_fminimizer struct. The updating procedure uses
only function evaluations (not derivatives).

## Caveats

Note that minimization functions can only search for one minimum at a time. When there are several
minima in the search area, the first minimum to be found will be returned; however it is difficult
to predict which of the minima this will be. In most cases, no error will be reported if you try to
find a minimum in an area where there is more than one.

With all minimization algorithms it can be difficult to determine the location of the minimum to
full numerical precision. The behavior of the function in the region of the minimum x^* can be
approximated by a Taylor expansion,

y = f(x^*) + (1/2) f''(x^*) (x - x^*)^2

and the second term of this expansion can be lost when added to the first term at finite precision.
This magnifies the error in locating x^*, making it proportional to sqrt epsilon (where epsilon
is the relative accuracy of the floating point numbers). For functions with higher order minima,
such as x^4, the magnification of the error is correspondingly worse. The best that can be achieved
is to converge to the limit of numerical accuracy in the function values, rather than the location
of the minimum itself.

## Providing the function to minimize

You must provide a continuous function of one variable for the minimizers to operate on. In order to
allow for general parameters the functions are defined by a gsl_function data type (see
[Providing the function to solve](http://www.gnu.org/software/gsl/manual/html_node/Providing-the-function-to-solve.html#Providing-the-function-to-solve)).

## Iteration

The following functions drive the iteration of each algorithm. Each function performs one iteration
to update the state of any minimizer of the corresponding type. The same functions work for all
minimizers so that different methods can be substituted at runtime without modifications to the
code.

## Stopping Parameters

A minimization procedure should stop when one of the following conditions is true:

 * A minimum has been found to within the user-specified precision.
 * A user-specified maximum number of iterations has been reached.
 * An error has occurred.

The handling of these conditions is under user control. The function below allows the user to test
the precision of the current result.

## Minimization Algorithms

The minimization algorithms described in this section require an initial interval which is
guaranteed to contain a minimum—if a and b are the endpoints of the interval and x is an estimate of
the minimum then f(a) > f(x) < f(b). This ensures that the function has at least one minimum
somewhere in the interval. If a valid initial interval is used then these algorithm cannot fail,
provided the function is well-behaved.
!*/

use crate::Value;
use ffi::FFI;

ffi_wrapper!(
    Minimizer<'a>,
    *mut sys::gsl_min_fminimizer,
    gsl_min_fminimizer_free
    ;inner_call: sys::gsl_function_struct => sys::gsl_function_struct { function: None, params: std::ptr::null_mut() };
    ;inner_closure: Option<Box<dyn Fn(f64) -> f64 + 'a>> => None;
);

impl<'a> Minimizer<'a> {
    /// This function returns a pointer to a newly allocated instance of a minimizer of type T. For
    /// example, the following code creates an instance of a golden section minimizer,
    ///
    /// ```C
    /// const gsl_min_fminimizer_type * T
    ///   = gsl_min_fminimizer_goldensection;
    /// gsl_min_fminimizer * s
    ///   = gsl_min_fminimizer_alloc (T);
    /// ```
    ///
    /// If there is insufficient memory to create the minimizer then the function returns a null
    /// pointer and the error handler is invoked with an error code of ::NoMem.
    #[doc(alias = "gsl_min_fminimizer_alloc")]
    pub fn new(t: MinimizerType) -> Option<Minimizer<'a>> {
        let ptr = unsafe { sys::gsl_min_fminimizer_alloc(t.unwrap_shared()) };

        if ptr.is_null() {
            None
        } else {
            Some(Self::wrap(ptr))
        }
    }

    /// This function sets, or resets, an existing minimizer s to use the function f and the initial
    /// search interval [x_lower, x_upper], with a guess for the location of the minimum x_minimum.
    ///
    /// If the interval given does not contain a minimum, then the function returns an error code of
    /// `Value::Invalid`.
    #[doc(alias = "gsl_min_fminimizer_set")]
    pub fn set<F: Fn(f64) -> f64 + 'a>(
        &mut self,
        f: F,
        x_minimum: f64,
        x_lower: f64,
        x_upper: f64,
    ) -> Result<(), Value> {
        self.inner_call = wrap_callback!(f, F + 'a);
        self.inner_closure = Some(Box::new(f));

        let ret = unsafe {
            sys::gsl_min_fminimizer_set(
                self.unwrap_unique(),
                &mut self.inner_call,
                x_minimum,
                x_lower,
                x_upper,
            )
        };
        result_handler!(ret, ())
    }

    /// This function is equivalent to gsl_min_fminimizer_set but uses the values f_minimum, f_lower
    /// and f_upper instead of computing f(x_minimum), f(x_lower) and f(x_upper).
    #[doc(alias = "gsl_min_fminimizer_set_with_values")]
    pub fn set_with_values<F: Fn(f64) -> f64 + 'a>(
        &mut self,
        f: F,
        x_minimum: f64,
        f_minimum: f64,
        x_lower: f64,
        f_lower: f64,
        x_upper: f64,
        f_upper: f64,
    ) -> Result<(), Value> {
        self.inner_call = wrap_callback!(f, F + 'a);
        self.inner_closure = Some(Box::new(f));

        let ret = unsafe {
            sys::gsl_min_fminimizer_set_with_values(
                self.unwrap_unique(),
                &mut self.inner_call,
                x_minimum,
                f_minimum,
                x_lower,
                f_lower,
                x_upper,
                f_upper,
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_min_fminimizer_name")]
    pub fn name(&self) -> Option<String> {
        let n = unsafe { sys::gsl_min_fminimizer_name(self.unwrap_shared()) };
        if n.is_null() {
            return None;
        }
        let mut len = 0;
        loop {
            if unsafe { *n.offset(len) } == 0 {
                break;
            }
            len += 1;
        }
        let slice = unsafe { std::slice::from_raw_parts(n as _, len as _) };
        std::str::from_utf8(slice).ok().map(|x| x.to_owned())
    }

    #[doc(alias = "gsl_min_fminimizer_x_minimum")]
    pub fn x_minimum(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_x_minimum(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_min_fminimizer_x_lower")]
    pub fn x_lower(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_x_lower(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_min_fminimizer_x_upper")]
    pub fn x_upper(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_x_upper(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_min_fminimizer_f_minimum")]
    pub fn f_minimum(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_f_minimum(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_min_fminimizer_f_lower")]
    pub fn f_lower(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_f_lower(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_min_fminimizer_f_upper")]
    pub fn f_upper(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_f_upper(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_min_fminimizer_minimum")]
    pub fn minimum(&self) -> f64 {
        unsafe { sys::gsl_min_fminimizer_minimum(self.unwrap_shared()) }
    }

    /// This function performs a single iteration of the minimizer s. If the iteration encounters an
    /// unexpected problem then an error code will be returned,
    ///
    /// `Value::BadFunc`
    /// the iteration encountered a singular point where the function evaluated to Inf or NaN.
    ///
    /// `Value::Failure`
    /// the algorithm could not improve the current best approximation or bounding interval.
    ///
    /// The minimizer maintains a current best estimate of the position of the minimum at all times,
    /// and the current interval bounding the minimum. This information can be accessed with the
    /// following auxiliary functions,
    #[doc(alias = "gsl_min_fminimizer_iterate")]
    pub fn iterate(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_min_fminimizer_iterate(self.unwrap_unique()) };
        result_handler!(ret, ())
    }
}

ffi_wrapper!(MinimizerType, *const sys::gsl_min_fminimizer_type);

impl MinimizerType {
    #[doc(alias = "gsl_min_fminimizer_goldensection")]
    pub fn goldensection() -> Self {
        ffi_wrap!(gsl_min_fminimizer_goldensection)
    }

    #[doc(alias = "gsl_min_fminimizer_brent")]
    pub fn brent() -> Self {
        ffi_wrap!(gsl_min_fminimizer_brent)
    }

    #[doc(alias = "gsl_min_fminimizer_quad_golden")]
    pub fn quad_golden() -> Self {
        ffi_wrap!(gsl_min_fminimizer_quad_golden)
    }
}

#[cfg(any(test, doctest))]
mod test {
    /// This doc block will be used to ensure that the closure can't be set everywhere!
    ///
    /// ```compile_fail
    /// use rgsl::*;
    /// use rgsl::minimizer::test_interval;
    ///
    /// fn set(min: &mut Minimizer) {
    ///     let y = "lalal".to_owned();
    ///     min.set(|x| {
    ///         println!("==> {:?}", y);
    ///         x * x - 5.
    ///     }, 1.0, -5.0, 5.0);
    /// }
    ///
    /// let mut min = Minimizer::new(MinimizerType::brent()).unwrap();
    /// set(&mut min);
    /// let status = min.iterate();
    /// ```
    ///
    /// Same but a working version:
    ///
    /// ```
    /// use rgsl::*;
    /// use rgsl::minimizer::test_interval;
    ///
    /// fn set(min: &mut Minimizer) {
    ///     min.set(|x| x * x - 5., 1.0, -5.0, 5.0);
    /// }
    ///
    /// let mut min = Minimizer::new(MinimizerType::brent()).unwrap();
    /// set(&mut min);
    /// let status = min.iterate();
    /// ```
    use super::*;
    use minimizer::test_interval;

    fn quadratic_test_fn(x: f64) -> f64 {
        x.powf(2.0) - 5.0
    }

    #[test]
    fn test_min() {
        let mut min = Minimizer::new(MinimizerType::brent()).unwrap();
        min.set(quadratic_test_fn, 1.0, -5.0, 5.0).unwrap();

        let max_iter = 5_usize;
        let eps_abs = 0.0001;
        let eps_rel = 0.0000001;

        let mut status = Value::Continue;
        let mut iter = 0_usize;

        while matches!(status, Value::Continue) && iter < max_iter {
            // iterate for next value
            min.iterate().unwrap(); // fails here w/ segfault

            // test for convergence
            let r = min.minimum();
            let x_lo = min.x_lower();
            let x_hi = min.x_upper();

            status = test_interval(x_lo, x_hi, eps_abs, eps_rel);

            // check if iteration succeeded
            if matches!(status, Value::Success) {
                println!("Converged");
            }

            // print current iteration
            println!("{} [{}, {}] {} {}", iter, x_lo, x_hi, r, x_hi - x_lo);

            iter += 1;
        }
    }
}
