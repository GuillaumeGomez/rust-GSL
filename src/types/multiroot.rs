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

use crate::ffi::FFI;
use crate::{Error, MatrixF64, VectorF64, View};
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
    /// This is a version of the Hybrid algorithm which replaces calls to the Jacobian function by
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
    ;inner_closure: Option<Box<dyn Fn(&VectorF64, &mut VectorF64) -> Result<(), Error> + 'a>> => None;,
    "This is a workspace for multidimensional root-finding without derivatives."
);

impl<'a> MultiRootFSolver<'a> {
    /// This function returns a pointer to a newly allocated instance of a solver of type `T` with
    /// `n` unknowns.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Error::NoMemory`.
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
    pub fn set<F: Fn(&VectorF64, &mut VectorF64) -> Result<(), Error> + 'a>(
        &mut self,
        f: F,
        n: usize,
        x: &VectorF64,
    ) -> Result<(), Error> {
        unsafe extern "C" fn inner_f<A>(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            f: *mut sys::gsl_vector,
        ) -> c_int
        where
            A: Fn(&VectorF64, &mut VectorF64) -> Result<(), Error>,
        {
            let g: &A = &*(params as *const A);
            let x_new = VectorF64::soft_wrap(x as *const _ as *mut _);
            Error::to_c(g(&x_new, &mut VectorF64::soft_wrap(f)))
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
        Error::handle(ret, ())
    }

    /// This function performs a single iteration of the minimizer s. If the iteration encounters an
    /// unexpected problem then an error code will be returned,
    ///
    /// `Error::BadFunc`
    /// the iteration encountered a singular point where the function evaluated to Inf or NaN.
    ///
    /// `Error::Failure`
    /// the algorithm could not improve the current best approximation or bounding interval.
    ///
    /// The minimizer maintains a current best estimate of the position of the minimum at all times,
    /// and the current interval bounding the minimum. This information can be accessed with the
    /// following auxiliary functions,
    #[doc(alias = "gsl_multiroot_fsolver_iterate")]
    pub fn iterate(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_multiroot_fsolver_iterate(self.unwrap_unique()) };
        Error::handle(ret, ())
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

ffi_wrapper!(
    MultiRootFdfSolverType,
    *const sys::gsl_multiroot_fdfsolver_type
);

impl MultiRootFdfSolverType {
    /// This is a modified version of Powell’s Hybrid method as implemented in the HYBRJ algorithm in MINPACK.
    /// Minpack was written by Jorge J. Moré, Burton S. Garbow and Kenneth E. Hillstrom.
    /// The Hybrid algorithm retains the fast convergence of Newton’s method but will also reduce the
    /// residual when Newton’s method is unreliable.
    ///
    /// The algorithm uses a generalized trust region to keep each step under control.
    /// In order to be accepted a proposed new position x' must satisfy the condition |D (x' - x)| < \delta,
    /// where D is a diagonal scaling matrix and \delta is the size of the trust region.
    /// The components of D are computed internally, using the column norms of the Jacobian to estimate the
    /// sensitivity of the residual to each component of x. This improves the behavior of the algorithm for badly scaled functions.
    ///
    /// On each iteration the algorithm first determines the standard Newton step by solving the system
    /// J dx = - f. If this step falls inside the trust region it is used as a trial step in the
    /// next stage. If not, the algorithm uses the linear combination of the Newton and gradient directions
    /// which is predicted to minimize the norm of the function while staying inside the trust region,
    /// dx = - \alpha J^{-1} f(x) - \beta \nabla |f(x)|^2
    /// This combination of Newton and gradient directions is referred to as a dogleg step.
    ///
    /// The proposed step is now tested by evaluating the function at the resulting point, x'.
    /// If the step reduces the norm of the function sufficiently then it is accepted and size of
    /// the trust region is increased. If the proposed step fails to improve the solution then the
    /// size of the trust region is decreased and another trial step is computed.
    ///
    /// The speed of the algorithm is increased by computing the changes to the Jacobian approximately,
    /// using a rank-1 update. If two successive attempts fail to reduce the residual then the full
    /// Jacobian is recomputed. The algorithm also monitors the progress of the solution and returns an error if several steps fail to make any improvement,
    /// `crate::Error::NoProgress` the iteration is not making any progress, preventing the algorithm from continuing.
    /// `crate::Error::NoProgressJacobian` re-evaluations of the Jacobian indicate that the iteration is not making any progress, preventing the algorithm from continuing.
    #[doc(alias = "gsl_multiroot_fdfsolver_hybridsj")]
    pub fn hybridsj() -> Self {
        ffi_wrap!(gsl_multiroot_fdfsolver_hybridsj)
    }

    /// This algorithm is an unscaled version of `hybridsj`. The steps are controlled by a spherical trust region
    /// |x' - x| < \delta, instead of a generalized region. This can be useful if the generalized region estimated
    /// by `hybridsj` is inappropriate.
    #[doc(alias = "gsl_multiroot_fdfsolver_hybridj")]
    pub fn hybridj() -> Self {
        ffi_wrap!(gsl_multiroot_fdfsolver_hybridj)
    }

    /// Newton’s Method is the standard root-polishing algorithm. The algorithm begins with
    /// an initial guess for the location of the solution. On each iteration a linear approximation to
    /// the function F is used to estimate the step which will zero all the components of the residual.
    /// The iteration is defined by the following sequence, x \to x' = x - J^{-1} f(x)
    /// where the Jacobian matrix J is computed from the derivative functions provided by f.
    /// The step dx is obtained by solving the linear system, J dx = - f(x)
    /// using LU decomposition. If the Jacobian matrix is singular, an error code of
    /// `crate::Error::Domain` is returned.
    #[doc(alias = "gsl_multiroot_fdfsolver_newton")]
    pub fn newton() -> Self {
        ffi_wrap!(gsl_multiroot_fdfsolver_newton)
    }

    /// This is a modified version of Newton’s method which attempts to improve global convergence
    /// by requiring every step to reduce the Euclidean norm of the residual, |f(x)|.
    /// If the Newton step leads to an increase in the norm then a reduced step of relative size,
    /// t = (\sqrt{1 + 6 r} - 1) / (3 r)
    /// is proposed, with r being the ratio of norms |f(x')|^2/|f(x)|^2.
    /// This procedure is repeated until a suitable step size is found.
    #[doc(alias = "gsl_multiroot_fdfsolver_gnewton")]
    pub fn gnewton() -> Self {
        ffi_wrap!(gsl_multiroot_fdfsolver_gnewton)
    }
}

pub struct MultiRootFdfSolverFunction<'a> {
    pub f: Box<dyn Fn(&VectorF64, &mut VectorF64) -> Result<(), Error> + 'a>,
    pub df: Box<dyn Fn(&VectorF64, &mut MatrixF64) -> Result<(), Error> + 'a>,
    pub fdf: Box<dyn Fn(&VectorF64, &mut VectorF64, &mut MatrixF64) -> Result<(), Error> + 'a>,
    pub n: usize,
    intern: sys::gsl_multiroot_function_fdf,
}

impl<'a> MultiRootFdfSolverFunction<'a> {
    #[doc(alias = "gsl_multiroot_function_fdf")]
    pub fn new<
        F: Fn(&VectorF64, &mut VectorF64) -> Result<(), Error> + 'a,
        DF: Fn(&VectorF64, &mut MatrixF64) -> Result<(), Error> + 'a,
        FDF: Fn(&VectorF64, &mut VectorF64, &mut MatrixF64) -> Result<(), Error> + 'a,
    >(
        f: F,
        df: DF,
        fdf: FDF,
        n: usize,
    ) -> MultiRootFdfSolverFunction<'a> {
        unsafe extern "C" fn inner_f(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            f: *mut sys::gsl_vector,
        ) -> i32 {
            let t = &*(params as *mut MultiRootFdfSolverFunction);
            let i_f = &t.f;
            Error::to_c(i_f(
                &VectorF64::soft_wrap(x as *const _ as *mut _),
                &mut VectorF64::soft_wrap(f as *const _ as *mut _),
            ))
        }

        unsafe extern "C" fn inner_df(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            J: *mut sys::gsl_matrix,
        ) -> i32 {
            let t = &*(params as *mut MultiRootFdfSolverFunction);
            let i_df = &t.df;
            Error::to_c(i_df(
                &VectorF64::soft_wrap(x as *const _ as *mut _),
                &mut MatrixF64::soft_wrap(J as *const _ as *mut _),
            ))
        }

        unsafe extern "C" fn inner_fdf(
            x: *const sys::gsl_vector,
            params: *mut c_void,
            f: *mut sys::gsl_vector,
            J: *mut sys::gsl_matrix,
        ) -> i32 {
            let t = &*(params as *mut MultiRootFdfSolverFunction);
            let i_fdf = &t.fdf;
            Error::to_c(i_fdf(
                &VectorF64::soft_wrap(x as *const _ as *mut _),
                &mut VectorF64::soft_wrap(f as *const _ as *mut _),
                &mut MatrixF64::soft_wrap(J as *const _ as *mut _),
            ))
        }

        MultiRootFdfSolverFunction {
            f: Box::new(f),
            df: Box::new(df),
            fdf: Box::new(fdf),
            n,
            intern: sys::gsl_multiroot_function_fdf {
                f: Some(inner_f),
                df: Some(inner_df),
                fdf: Some(inner_fdf),
                n,
                params: std::ptr::null_mut(),
            },
        }
    }

    #[allow(clippy::wrong_self_convention)]
    fn to_raw(&mut self) -> *mut sys::gsl_multiroot_function_fdf {
        self.intern.n = self.n;
        self.intern.params = self as *mut MultiRootFdfSolverFunction as *mut c_void;
        &mut self.intern
    }
}

ffi_wrapper!(
    MultiRootFdfSolver,
    *mut sys::gsl_multiroot_fdfsolver,
    gsl_multiroot_fdfsolver_free
);

impl MultiRootFdfSolver {
    /// This function returns a pointer to a newly allocated instance of a derivative solver of type T
    /// for a system of `n` dimensions.
    #[doc(alias = "gsl_multiroot_fdfsolver_alloc")]
    pub fn new(t: MultiRootFdfSolverType, n: usize) -> Option<MultiRootFdfSolver> {
        let ptr = unsafe { sys::gsl_multiroot_fdfsolver_alloc(t.unwrap_shared(), n) };

        if ptr.is_null() {
            None
        } else {
            Some(Self::wrap(ptr))
        }
    }

    /// These functions set, or reset, an existing solver to use the functions `f`, and the initial guess `x.
    #[doc(alias = "gsl_multiroot_fdfsolver_set")]
    pub fn set(&mut self, f: &mut MultiRootFdfSolverFunction, x: &VectorF64) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_multiroot_fdfsolver_set(self.unwrap_unique(), f.to_raw(), x.unwrap_shared())
        };
        Error::handle(ret, ())
    }

    /// Return the name of the solver.
    #[doc(alias = "gsl_multiroot_fdfsolver_name")]
    pub fn name(&self) -> Option<String> {
        let n = unsafe { sys::gsl_multiroot_fdfsolver_name(self.unwrap_shared()) };
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

    /// Perform a single iteration of the solver. If the iteration
    /// encounters an unexpected problem then an error code will be
    /// returned,
    ///
    /// * `crate::Error::BadFunc` the iteration encountered a singular
    ///   point where the function or its derivative evaluated to Inf or NaN.
    ///
    /// * `crate::Error::NoProgress` the iteration is not making any progress,
    ///    preventing the algorithm from continuing.
    ///
    /// The solver maintains a current best estimate of the root and
    /// its function value at all times.  This information can be
    /// accessed with `root`, `f`, and `dx` functions.
    #[doc(alias = "gsl_multiroot_fdfsolver_iterate")]
    pub fn iterate(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_multiroot_fdfsolver_iterate(self.unwrap_unique()) };
        Error::handle(ret, ())
    }

    /// Returns the current estimate of the root for the solver.
    #[doc(alias = "gsl_multiroot_fdfsolver_root")]
    pub fn root(&self) -> View<'_, VectorF64> {
        unsafe { View::new(sys::gsl_multiroot_fdfsolver_root(self.unwrap_shared())) }
    }

    /// Returns the last step taken by the solver.
    #[doc(alias = "gsl_multiroot_fdfsolver_dx")]
    pub fn dx(&self) -> View<'_, VectorF64> {
        unsafe { View::new(sys::gsl_multiroot_fdfsolver_dx(self.unwrap_shared())) }
    }

    /// Returns the function value f(x) at the current estimate of the root for the solver.
    #[doc(alias = "gsl_multiroot_fdfsolver_f")]
    pub fn f(&self) -> View<'_, VectorF64> {
        unsafe { View::new(sys::gsl_multiroot_fdfsolver_f(self.unwrap_shared())) }
    }
}

#[cfg(any(test, doctest))]
mod tests {
    /// This doc block will be used to ensure that the closure can't be set everywhere!
    ///
    /// ```compile_fail
    /// use crate::rgsl::*;
    /// use crate::rgsl::types::multiroot::{MultiRootFSolver, MultiRootFSolverType};
    ///
    /// fn set(root: &mut MultiRootFSolver) {
    ///     let dummy = "lalal".to_owned();
    ///     root.set(|x, y| {
    ///         println!("==> {:?}", dummy);
    ///         y.set(0, 1.0 - x.get(0));
    ///         y.set(1, x.get(0) - x.get(1));
    ///         rgsl::Error::Success}, 2, &rgsl::VectorF64::from_slice(&[-10.0, 1.0]).unwrap());
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
    /// use crate::rgsl::types::multiroot::{MultiRootFSolver, MultiRootFSolverType};
    /// use crate::rgsl::*;
    ///
    /// fn set(root: &mut MultiRootFSolver) {
    ///     root.set(|x, y| {
    ///         y.set(0, 1.0 - x.get(0));
    ///         y.set(1, x.get(0) - x.get(1));
    ///         Ok(())
    ///     },
    ///     2, &rgsl::VectorF64::from_slice(&[-10.0, 1.0]).unwrap());
    /// }
    ///
    /// let mut root = MultiRootFSolver::new(&MultiRootFSolverType::hybrid(), 2).unwrap();
    /// set(&mut root);
    /// let status = root.iterate();
    /// ```
    ///
    use super::*;
    use crate::multiroot::test_residual;

    const RPARAMS: (f64, f64) = (1.0, 10.0);

    /// checking a test function
    /// must return a success criteria (or failure)
    fn rosenbrock_f(x: &VectorF64, f: &mut VectorF64) -> Result<(), Error> {
        f.set(0, RPARAMS.0 * (1.0 - x.get(0)));
        f.set(1, RPARAMS.1 * (x.get(1) - x.get(0).powf(2.0)));
        Ok(())
    }

    fn rosenbrock_df(x: &VectorF64, J: &mut MatrixF64) -> Result<(), Error> {
        J.set(0, 0, -RPARAMS.0);
        J.set(0, 1, 0f64);
        J.set(1, 0, -2.0 * RPARAMS.1 * x.get(0));
        J.set(1, 1, RPARAMS.1);
        Ok(())
    }

    fn rosenbrock_fdf(x: &VectorF64, f: &mut VectorF64, J: &mut MatrixF64) -> Result<(), Error> {
        rosenbrock_f(x, f)?;
        rosenbrock_df(x, J)?;
        Ok(())
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

    fn print_fdf_state(solver: &mut MultiRootFdfSolver, iteration: usize) {
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
        let mut status = Err(Error::Continue);
        let epsabs = 1e-6;

        print_state(&mut multi_root, 0);

        while matches!(status, Err(Error::Continue)) && iter < max_iter {
            // iterate solver
            multi_root.iterate().unwrap();

            // print current iteration
            print_state(&mut multi_root, iter);

            // test for convergence
            let f_value = multi_root.f();
            status = test_residual(&f_value, epsabs);
            if status.is_ok() {
                println!("Converged");
            }

            iter += 1;
        }
        assert!(status.is_ok())
    }

    #[test]
    fn test_multiroot_fdf_fsolver() {
        // setup workspace
        let mut multi_root = MultiRootFdfSolver::new(MultiRootFdfSolverType::gnewton(), 2).unwrap();
        let array_size: usize = 2;
        let guess_value = VectorF64::from_slice(&[-10.0, -5.0]).unwrap();
        let mut fs = MultiRootFdfSolverFunction::new(
            rosenbrock_f,
            rosenbrock_df,
            rosenbrock_fdf,
            array_size,
        );
        multi_root.set(&mut fs, &guess_value).unwrap();

        // iteration counters
        let max_iter: usize = 100;
        let mut iter = 0;

        // convergence checks
        let mut status = Err(Error::Continue);
        let epsabs = 1e-6;

        print_fdf_state(&mut multi_root, 0);

        while matches!(status, Err(Error::Continue)) && iter < max_iter {
            // iterate solver
            multi_root.iterate().unwrap();

            // print current iteration
            print_fdf_state(&mut multi_root, iter);

            // test for convergence
            let f_value = multi_root.f();
            status = test_residual(&f_value, epsabs);
            if status.is_ok() {
                println!("Converged");
            }

            iter += 1;
        }
        assert!(status.is_ok())
    }
}
