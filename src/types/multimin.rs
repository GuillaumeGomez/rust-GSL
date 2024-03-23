//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!

# Multidimensional Minimization

This chapter describes routines for finding minima of arbitrary multidimensional functions. The
library provides low level components for a variety of iterative minimizers and convergence tests.
These can be combined by the user to achieve the desired solution, while providing full access
to the intermediate steps of the algorithms. Each class of methods uses the same framework,
so that you can switch between minimizers at runtime without needing to recompile your program.
Each instance of a minimizer keeps track of its own state, allowing the minimizers to be used
in multi-threaded programs. The minimization algorithms can be used to maximize a function by
inverting its sign.

## Overview

The problem of multidimensional minimization requires finding a point x such that the scalar function,

```text
f(x_1, ... , x_n)
```

takes a value which is lower than at any neighboring point. For smooth functions the gradient g = \nabla f
vanishes at the minimum. In general there are no bracketing methods available for the minimization of
n-dimensional functions. The algorithms proceed from an initial guess using a search algorithm which
attempts to move in a downhill direction.

Algorithms making use of the gradient of the function perform a one-dimensional line minimisation along this
direction until the lowest point is found to a suitable tolerance. The search direction is then
updated with local information from the function and its derivatives, and the whole process repeated
until the true n-dimensional minimum is found.

Algorithms which do not require the gradient of the function use different strategies. For example,
the Nelder-Mead Simplex algorithm maintains n+1 trial parameter vectors as the vertices of a n-dimensional
simplex. On each iteration it tries to improve the worst vertex of the simplex by geometrical transformations.
The iterations are continued until the overall size of the simplex has decreased sufficiently.

Both types of algorithms use a standard framework. The user provides a high-level driver for the algorithms,
and the library provides the individual functions necessary for each of the steps. There are three
main phases of the iteration. The steps are,

 * initialize minimizer state, s, for algorithm T
 * update s using the iteration T
 * test s for convergence, and repeat iteration if necessary

Each iteration step consists either of an improvement to the line-minimisation in the current
direction or an update to the search direction itself. The state for the minimizers is held
in a gsl_multimin_fdfminimizer struct or a gsl_multimin_fminimizer struct.

## Caveats

Note that the minimization algorithms can only search for one local minimum at a time. When there are several
local minima in the search area, the first minimum to be found will be returned; however it is difficult
to predict which of the minima this will be. In most cases, no error will be reported if you try to
find a local minimum in an area where there is more than one.

It is also important to note that the minimization algorithms find local minima; there is no way
to determine whether a minimum is a global minimum of the function in question.

## Providing a function to minimize

You must provide a parametric function of n variables for the minimizers to operate on.

*/

use crate::ffi::FFI;
use crate::{Value, VectorF64, View};
use sys::libc::c_void;

ffi_wrapper!(
    Minimizer<'a>,
    *mut sys::gsl_multimin_fminimizer,
    gsl_multimin_fminimizer_free
    ;inner_call: sys::gsl_multimin_function_struct => sys::gsl_multimin_function_struct { f: None, n: 0_usize, params: std::ptr::null_mut() };
    ;inner_closure: Option<Box<dyn Fn(&VectorF64) -> f64 + 'a>> => None;
);

impl<'a> Minimizer<'a> {
    /// Creates a minimizer of type `t` for an n-dimensional function.
    /// If there is insufficient memory to create the minimizer then
    /// the function returns a null pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    #[doc(alias = "gsl_multimin_fminimizer_alloc")]
    pub fn new(t: MinimizerType, n: usize) -> Option<Minimizer<'a>> {
        let ptr = unsafe { sys::gsl_multimin_fminimizer_alloc(t.unwrap_shared(), n) };

        if ptr.is_null() {
            None
        } else {
            Some(Self::wrap(ptr))
        }
    }

    /// This function initializes the minimizer `s` to minimize the function `f`, starting from the initial point `x`.
    /// The size of the initial trial steps is given in vector `step_size`. The precise meaning of this
    /// parameter depends on the method used.
    #[doc(alias = "gsl_multimin_fminimizer_set")]
    pub fn set<F: Fn(&VectorF64) -> f64 + 'a>(
        &mut self,
        f: F,
        x: &VectorF64,
        step_size: &VectorF64,
    ) -> Result<(), Value> {
        unsafe extern "C" fn inner_f<F: Fn(&VectorF64) -> f64>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
        ) -> f64 {
            let f: &F = &*(params as *const F);
            let x_new = VectorF64::soft_wrap(x as *const _ as *mut _);
            f(&x_new)
        }
        self.inner_call = sys::gsl_multimin_function_struct {
            f: Some(inner_f::<F>),
            n: x.len(),
            params: &f as *const _ as *mut _,
        };

        self.inner_closure = Some(Box::new(f));

        let ret = unsafe {
            sys::gsl_multimin_fminimizer_set(
                self.unwrap_unique(),
                &mut self.inner_call,
                x.unwrap_shared(),
                step_size.unwrap_shared(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function returns a pointer to the name of the minimizer.
    #[doc(alias = "gsl_multimin_fminimizer_name")]
    pub fn name(&self) -> Option<String> {
        let n = unsafe { sys::gsl_multimin_fminimizer_name(self.unwrap_shared()) };
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

    #[doc(alias = "gsl_multimin_fminimizer_x")]
    pub fn x(&self) -> View<'_, VectorF64> {
        unsafe { View::new(sys::gsl_multimin_fminimizer_x(self.unwrap_shared())) }
    }

    #[doc(alias = "gsl_multimin_fminimizer_minimum")]
    pub fn minimum(&self) -> f64 {
        unsafe { sys::gsl_multimin_fminimizer_minimum(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_multimin_fminimizer_size")]
    pub fn size(&self) -> f64 {
        unsafe { sys::gsl_multimin_fminimizer_size(self.unwrap_shared()) }
    }

    /// This function performs a single iteration of the minimizer s. If the iteration encounters
    /// an unexpected problem then an error code will be returned. The error code `Value::NoProgress`
    /// signifies that the minimizer is unable to improve on its current estimate, either due
    /// to numerical difficulty or because a genuine local minimum has been reached.
    #[doc(alias = "gsl_multimin_fminimizer_iterate")]
    pub fn iterate(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_multimin_fminimizer_iterate(self.unwrap_unique()) };
        result_handler!(ret, ())
    }
}

ffi_wrapper!(MinimizerType, *const sys::gsl_multimin_fminimizer_type);

impl MinimizerType {
    #[doc(alias = "gsl_multimin_fminimizer_nmsimplex2")]
    pub fn nm_simplex2() -> Self {
        ffi_wrap!(gsl_multimin_fminimizer_nmsimplex2)
    }

    #[doc(alias = "gsl_multimin_fminimizer_nmsimplex")]
    pub fn nm_simplex() -> Self {
        ffi_wrap!(gsl_multimin_fminimizer_nmsimplex)
    }

    #[doc(alias = "gsl_multimin_fminimizer_nmsimplex2rand")]
    pub fn nm_simplex2_rand() -> Self {
        ffi_wrap!(gsl_multimin_fminimizer_nmsimplex2rand)
    }
}

#[cfg(any(test, doctest))]
mod test {
    /// This doc block will be used to ensure that the closure can't be set everywhere!
    ///
    /// ```compile_fail
    /// extern crate rgsl;
    /// use crate::rgsl::types::multimin::{Minimizer,MinimizerType};
    ///
    /// fn set(m: &mut Minimizer) {
    ///     let dummy = "lalal".to_owned();
    ///     m.set(|x| {
    ///     println!("==> {:?}", dummy);
    ///     x.get(0) + x.get(1)},
    ///     &rgsl::VectorF64::from_slice(&[-10.0, 1.0]).unwrap(),
    ///     &rgsl::VectorF64::from_slice(&[1.0, 1.0]).unwrap()).unwrap();
    ///
    ///     let mut mint = Minimizer::new(MinimizerType::nm_simplex(), 2).unwrap();
    ///     set(&mut mint);
    ///     let _status = mint.iterate();
    /// ```
    ///
    /// Same but a working version:
    ///
    /// ```
    /// extern crate rgsl;
    /// use crate::rgsl::types::multimin::{Minimizer,MinimizerType};
    ///
    /// fn set(m: &mut Minimizer) {
    ///     m.set(|x| {
    ///     x.get(0) + x.get(1)},
    ///     &rgsl::VectorF64::from_slice(&[-10.0, 1.0]).unwrap(),
    ///     &rgsl::VectorF64::from_slice(&[1.0, 1.0]).unwrap()).unwrap();
    ///
    ///     let mut mint = Minimizer::new(MinimizerType::nm_simplex(), 2).unwrap();
    ///     set(&mut mint);
    ///     let _status = mint.iterate();
    /// }
    /// ```
    use super::*;
    use crate::multimin::test_size;

    fn print_state(min: &Minimizer, iter: usize) {
        let f = min.minimum();
        let x = min.x();
        println!(
            "iter: {}, f = {:+.2e}, x = [{:+.5}, {:+.5}]",
            iter,
            f,
            x.get(0),
            x.get(1),
        )
    }

    #[test]
    fn test_multi_min() {
        let mut min = Minimizer::new(MinimizerType::nm_simplex2(), 2).unwrap();
        let guess_value = VectorF64::from_slice(&[5.0, 7.0]).unwrap();
        let step_size = VectorF64::from_slice(&[1.0, 1.0]).unwrap();

        let paraboloid = |v: &VectorF64| -> f64 {
            let center = (1.0, 2.0);
            let scale = (10.0, 20.0);
            let minimum = 30.0;

            let x = v.get(0);
            let y = v.get(1);
            let result =
                scale.0 * (x - center.0).powf(2.0) + scale.1 * (y - center.1).powf(2.0) + minimum;
            result
        };

        min.set(paraboloid, &guess_value, &step_size).unwrap();

        let max_iter = 100_usize;
        let eps_abs = 0.01;

        let mut status = Value::Continue;
        let mut iter = 0_usize;

        while matches!(status, Value::Continue) && iter < max_iter {
            // iterate for next value
            min.iterate().unwrap(); // fails here w/ segfault

            // test for convergence
            let size = min.size();

            status = test_size(size, eps_abs);

            // check if iteration succeeded
            if matches!(status, Value::Success) {
                println!("Converged");
            }

            // print current iteration
            print_state(&min, iter);

            iter += 1;
        }
    }
}
