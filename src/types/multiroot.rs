//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Multi dimensional Root-Finding

This chapter describes functions for multidimensional root-finding (solving nonlinear systems with
n equations in n unknowns). The library provides low level components for a variety of iterative
solvers and convergence tests. These can be combined by the user to achieve the desired solution,
with full access to the intermediate steps of the iteration. Each class of methods uses the same
framework, so that you can switch between solvers at runtime without needing to recompile your
program. Each instance of a solver keeps track of its own state, allowing the solvers to be used
in multi-threaded programs. The solvers are based on the original Fortran library MINPACK.
The header file `gsl_multiroots.h` contains prototypes for the multidimensional root finding functions
and related declarations.

## Overview
The problem of multidimensional root finding requires the simultaneous solution of n equations,
f_i, in n variables, x_i,

```text
f_i (x_1, \dots, x_n) = 0 \qquad\hbox{for}~i = 1 \dots n.
```

In general there are no bracketing methods available for n dimensional systems, and no way of
knowing whether any solutions exist. All algorithms proceed from an initial guess using a
variant of the Newton iteration,

```text
x \to x' = x - J^{-1} f(x)
```
where x, f are vector quantities and J is the Jacobian matrix J_{ij} = \partial f_i / \partial x_j.
Additional strategies can be used to enlarge the region of convergence. These include requiring a
decrease in the norm |f| on each step proposed by Newton’s method, or taking steepest-descent
steps in the direction of the negative gradient of |f|.

Several root-finding algorithms are available within a single framework. The user provides a
high-level driver for the algorithms, and the library provides the individual functions necessary
for each of the steps. There are three main phases of the iteration. The steps are,

- initialize solver state, `s`, for algorithm `T`
- update `s` using the iteration `T`
- test `s` for convergence, and repeat iteration if necessary

The evaluation of the Jacobian matrix can be problematic, either because programming the derivatives
is intractable or because computation of the n^2 terms of the matrix becomes too expensive.
For these reasons the algorithms provided by the library are divided into two classes according to
whether the derivatives are available or not.

The state for solvers with an analytic Jacobian matrix is held in a `gsl_multiroot_fdfsolver` struct.
The updating procedure requires both the function and its derivatives to be supplied by the user.

The state for solvers which do not use an analytic Jacobian matrix is held in a
`gsl_multiroot_fsolver` struct. The updating procedure uses only function evaluations (not derivatives).
The algorithms estimate the matrix J or J^{-1} by approximate methods.
!*/

use crate::{Value, VectorF64, View};
use ffi::FFI;
use sys::libc::{c_int, c_void};

ffi_wrapper!(
    MultiRootFSolverType,
    *const sys::gsl_multiroot_fsolver_type,
    "The multiroot algorithms described in this section do not require any derivative information to be
    supplied by the user. Any derivatives needed are approximated by finite differences.
    Note that if the finite-differencing step size chosen by these routines is inappropriate,
    an explicit user-supplied numerical derivative can always be used with
    derivative-based algorithms."
);

impl MultiRootFSolverType {
    ///This is a version of the Hybrid algorithm which replaces calls to the Jacobian function by
    /// its finite difference approximation. The finite difference approximation is computed
    /// using `gsl_multiroots_fdjac()` with a relative step size of `GSL_SQRT_DBL_EPSILON`.
    /// Note that this step size will not be suitable for all problems.
    #[doc(alias = "gsl_multiroot_fsolver_hybrids")]
    pub fn hybrids() -> MultiRootFSolverType {
        ffi_wrap!(gsl_multiroot_fsolver_hybrids)
    }

    ///This is a finite difference version of the Hybrid algorithm without internal scaling.
    #[doc(alias = "gsl_multiroot_fsolver_hybrid")]
    pub fn hybrid() -> MultiRootFSolverType {
        ffi_wrap!(gsl_multiroot_fsolver_hybrid)
    }

    /// The discrete Newton algorithm is the simplest method of solving a multidimensional
    /// system. It uses the Newton iteration
    ///
    ///```text
    ///x \to x - J^{-1} f(x)
    ///```
    ///
    /// where the Jacobian matrix J is approximated by taking finite differences of the function f.
    /// The approximation scheme used by this implementation is,
    ///
    ///```text
    ///J_{ij} = (f_i(x + \delta_j) - f_i(x)) / \delta_j
    ///```
    ///
    /// where \delta_j is a step of size \sqrt\epsilon |x_j| with \epsilon being the machine
    /// precision (\epsilon \approx 2.22 \times 10^{-16}). The order of convergence of Newton’s
    /// algorithm is quadratic, but the finite differences require n^2 function evaluations on
    /// each iteration. The algorithm may become unstable if the finite differences are not a
    /// good approximation to the true derivatives.
    #[doc(alias = "gsl_multiroot_fsolver_dnewton")]
    pub fn dnewton() -> MultiRootFSolverType {
        ffi_wrap!(gsl_multiroot_fsolver_dnewton)
    }

    /// The Broyden algorithm is a version of the discrete Newton algorithm which attempts to
    /// avoids the expensive update of the Jacobian matrix on each iteration. The changes to
    /// the Jacobian are also approximated, using a rank-1 update,
    ///
    ///```text
    ///J^{-1} \to J^{-1} - (J^{-1} df - dx) dx^T J^{-1} / dx^T J^{-1} df
    ///```
    ///
    /// where the vectors dx and df are the changes in x and f. On the first iteration the inverse
    /// Jacobian is estimated using finite differences, as in the discrete Newton algorithm.
    ///
    /// This approximation gives a fast update but is unreliable if the changes are not small, and
    /// the estimate of the inverse Jacobian becomes worse as time passes. The algorithm has a
    /// tendency to become unstable unless it starts close to the root. The Jacobian is refreshed
    /// if this instability is detected (consult the source for details).
    ///
    /// This algorithm is included only for demonstration purposes, and is not recommended for
    /// serious use.
    #[doc(alias = "gsl_multiroot_fsolver_broyden")]
    pub fn broyden() -> MultiRootFSolverType {
        ffi_wrap!(gsl_multiroot_fsolver_broyden)
    }
}

ffi_wrapper!(
    MultiRootFSolver<'a>,
    *mut sys::gsl_multiroot_fsolver,
    gsl_multiroot_fsolver_free
    ;inner_call: sys::gsl_multiroot_function_struct => sys::gsl_multiroot_function_struct{ f: None, n: 0, params: std::ptr::null_mut() };
    ;inner_closure: Option<Box<dyn Fn(&VectorF64, &mut VectorF64) -> Value + 'a>> => None;,
    "This is a workspace for multidimensional root-finding without derivatives."
);

impl<'a> MultiRootFSolver<'a> {
    /// This function returns a pointer to a newly allocated instance of a solver of type `T` with
    /// `n` unknowns.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    #[doc(alias = "gsl_multiroot_fsolver_alloc")]
    pub fn new(t: &MultiRootFSolverType, n: usize) -> Option<MultiRootFSolver<'a>> {
        let ptr = unsafe { sys::gsl_multiroot_fsolver_alloc(t.unwrap_shared(), n) };

        if ptr.is_null() {
            None
        } else {
            Some(MultiRootFSolver::wrap(ptr))
        }
    }

    /// This function initializes, or reinitializes, an existing solver `s` to use the multi
    /// function `f` with `n` unknowns.
    #[doc(alias = "gsl_multiroot_fsolver_set")]
    pub fn set<F: Fn(&VectorF64, &mut VectorF64) -> Value + 'a>(
        &mut self,
        f: F,
        n: usize,
        x: &VectorF64,
    ) -> Result<(), Value> {
        unsafe extern "C" fn inner_f<A: Fn(&VectorF64, &mut VectorF64) -> Value>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            f: *mut sys::gsl_vector,
        ) -> c_int {
            let g: &A = &*(params as *const A);
            let x_new = VectorF64::soft_wrap(x as *const _ as *mut _);
            Value::into(g(&x_new, &mut VectorF64::soft_wrap(f)))
        }

        self.inner_call = sys::gsl_multiroot_function_struct {
            f: Some(inner_f::<F>),
            n,
            params: &f as *const _ as *mut _,
        };
        self.inner_closure = Some(Box::new(f));

        let ret = unsafe {
            sys::gsl_multiroot_fsolver_set(
                self.unwrap_unique(),
                &mut self.inner_call,
                x.unwrap_shared(),
            )
        };
        result_handler!(ret, ())
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
    #[doc(alias = "gsl_multiroot_fsolver_iterate")]
    pub fn iterate(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_multiroot_fsolver_iterate(self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function returns the current estimate of the root for the solver `s`, given by `s->x`.
    #[doc(alias = "gsl_multiroot_fsolver_root")]
    pub fn root(&self) -> View<'_, VectorF64> {
        unsafe { View::new(sys::gsl_multiroot_fsolver_root(self.unwrap_shared())) }
    }

    /// This function returns the last step `dx` taken by the solver `s`, given by `s->dx`.
    #[doc(alias = "gsl_multiroot_fsolver_dx")]
    pub fn dx(&self) -> View<'_, VectorF64> {
        unsafe { View::new(sys::gsl_multiroot_fsolver_dx(self.unwrap_shared())) }
    }

    /// This function returns the function value `f(x)` at the current estimate of the root for
    /// the solver `s`, given by `s->f`.
    #[doc(alias = "gsl_multiroot_fsolver_f")]
    pub fn f(&self) -> View<'_, VectorF64> {
        unsafe { View::new(sys::gsl_multiroot_fsolver_f(self.unwrap_shared())) }
    }
}

#[cfg(any(test, doctest))]
mod tests {
    /// This doc block will be used to ensure that the closure can't be set everywhere!
    ///
    /// ```compile_fail
    /// use rgsl::*;
    /// use rgsl::types::multiroot::{MultiRootFSolver, MultiRootFSolverType};
    ///
    /// fn set(root: &mut MultiRootFSolver) {
    ///     let dummy = "lalal".to_owned();
    ///     root.set(|x, y| {
    ///         println!("==> {:?}", dummy);
    ///         y.set(0, 1.0 - x.get(0));
    ///         y.set(1, x.get(0) - x.get(1));
    ///         rgsl::Value::Success}, 2, &rgsl::VectorF64::from_slice(&[-10.0, 1.0]).unwrap());
    /// }
    ///
    /// let mut root = MultiRootFSolver::new(&MultiRootFSolverType::hybrid(), 2).unwrap();
    /// set(&mut root);
    /// let status = root.iterate();
    /// ```
    ///
    /// Same but a working version:
    ///
    /// ```
    /// use rgsl::types::multiroot::{MultiRootFSolver, MultiRootFSolverType};
    /// use rgsl::*;
    ///
    /// fn set(root: &mut MultiRootFSolver) {
    ///     root.set(|x, y| {
    ///         y.set(0, 1.0 - x.get(0));
    ///         y.set(1, x.get(0) - x.get(1));
    ///         rgsl::Value::Success}, 2, &rgsl::VectorF64::from_slice(&[-10.0, 1.0]).unwrap());
    /// }
    ///
    /// let mut root = MultiRootFSolver::new(&MultiRootFSolverType::hybrid(), 2).unwrap();
    /// set(&mut root);
    /// let status = root.iterate();
    /// ```
    ///
    use super::*;
    use multiroot::test_residual;

    /// checking a test function
    /// must return a success criteria (or failure)
    fn rosenbrock_f(x: &VectorF64, f: &mut VectorF64) -> Value {
        f.set(0, 1.0 - x.get(0));
        f.set(1, x.get(0) - x.get(1).powf(2.0));
        Value::Success
    }

    fn print_state(solver: &mut MultiRootFSolver, iteration: usize) {
        let f = solver.f();
        let x = solver.root();
        println!(
            "iter: {}, f = [{:+.2e}, {:+.2e}], x = [{:+.5}, {:+.5}]",
            iteration,
            f.get(0),
            f.get(1),
            x.get(0),
            x.get(1)
        )
    }

    #[test]
    fn test_multiroot_fsolver() {
        // setup workspace
        let mut multi_root = MultiRootFSolver::new(&MultiRootFSolverType::hybrid(), 2).unwrap();
        let array_size: usize = 2;
        let guess_value = VectorF64::from_slice(&[-10.0, -5.0]).unwrap();
        multi_root
            .set(rosenbrock_f, array_size, &guess_value)
            .unwrap();

        // iteration counters
        let max_iter: usize = 100;
        let mut iter = 0;

        // convergence checks
        let mut status = crate::Value::Continue;
        let epsabs = 1e-6;

        print_state(&mut multi_root, 0);

        while matches!(status, crate::Value::Continue) && iter < max_iter {
            // iterate solver
            multi_root.iterate().unwrap();

            // print current iteration
            print_state(&mut multi_root, iter);

            // test for convergence
            let f_value = multi_root.f();
            status = test_residual(&f_value, epsabs);

            // check if iteration succeeded
            if matches!(status, crate::Value::Success) {
                println!("Converged");
            }

            iter += 1;
        }
        assert!(matches!(status, crate::Value::Success))
    }
}
