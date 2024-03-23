//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# One dimensional Root-Finding

This chapter describes routines for finding roots of arbitrary one-dimensional functions.
The library provides low level components for a variety of iterative solvers and convergence
tests. These can be combined by the user to achieve the desired solution, with full access to
the intermediate steps of the iteration. Each class of methods uses the same framework, so
that you can switch between solvers at runtime without needing to recompile your program.
Each instance of a solver keeps track of its own state, allowing the solvers to be used in
multi-threaded programs.

## Overview

One-dimensional root finding algorithms can be divided into two classes, root bracketing and
root polishing. Algorithms which proceed by bracketing a root are guaranteed to converge.
Bracketing algorithms begin with a bounded region known to contain a root. The size of
this bounded region is reduced, iteratively, until it encloses the root to a desired tolerance.
This provides a rigorous error estimate for the location of the root.

The technique of root polishing attempts to improve an initial guess to the root. These
algorithms converge only if started “close enough” to a root, and sacrifice a rigorous error
bound for speed. By approximating the behavior of a function in the vicinity of a root they
attempt to find a higher order improvement of an initial guess. When the behavior of the
function is compatible with the algorithm and a good initial guess is available a polishing
algorithm can provide rapid convergence.

In GSL both types of algorithm are available in similar frameworks. The user provides
a high-level driver for the algorithms, and the library provides the individual functions
necessary for each of the steps. There are three main phases of the iteration. The steps are,
• initialize solver state, s, for algorithm T
• update s using the iteration T
• test s for convergence, and repeat iteration if necessary

The state for bracketing solvers is held in a gsl_root_fsolver struct. The updating
procedure uses only function evaluations (not derivatives). The state for root polishing
solvers is held in a gsl_root_fdfsolver struct. The updates require both the function and
its derivative (hence the name fdf) to be supplied by the user.
!*/

use crate::ffi::FFI;
use crate::Value;
use sys::libc::{c_double, c_void};

ffi_wrapper!(
    RootFSolverType,
    *const sys::gsl_root_fsolver_type,
    "The root bracketing algorithms described in this section require an initial interval which is
guaranteed to contain a root—if a and b are the endpoints of the interval then f (a) must
differ in sign from f (b). This ensures that the function crosses zero at least once in the
interval. If a valid initial interval is used then these algorithm cannot fail, provided the
function is well-behaved.
Note that a bracketing algorithm cannot find roots of even degree, since these do not
cross the x-axis."
);

impl RootFSolverType {
    /// The bisection algorithm is the simplest method of bracketing the roots of a function.
    /// It is the slowest algorithm provided by the library, with linear convergence.
    /// On each iteration, the interval is bisected and the value of the function at the midpoint
    /// is calculated. The sign of this value is used to determine which half of the interval does
    /// not contain a root. That half is discarded to give a new, smaller interval containing
    /// the root. This procedure can be continued indefinitely until the interval is sufficiently
    /// small.
    ///
    /// At any time the current estimate of the root is taken as the midpoint of the interval.
    #[doc(alias = "gsl_root_fsolver_bisection")]
    pub fn bisection() -> RootFSolverType {
        ffi_wrap!(gsl_root_fsolver_bisection)
    }

    /// The false position algorithm is a method of finding roots based on linear interpolation.
    /// Its convergence is linear, but it is usually faster than bisection.
    ///
    /// On each iteration a line is drawn between the endpoints (a, f (a)) and (b, f (b)) and
    /// the point where this line crosses the x-axis taken as a “midpoint”. The value of the
    /// function at this point is calculated and its sign is used to determine which side of the
    /// interval does not contain a root. That side is discarded to give a new, smaller interval
    /// containing the root. This procedure can be continued indefinitely until the interval
    /// is sufficiently small.
    ///
    /// The best estimate of the root is taken from the linear interpolation of the interval on
    /// the current iteration.
    #[doc(alias = "gsl_root_fsolver_brent")]
    pub fn brent() -> RootFSolverType {
        ffi_wrap!(gsl_root_fsolver_brent)
    }

    /// The Brent-Dekker method (referred to here as Brent’s method) combines an interpo-
    /// lation strategy with the bisection algorithm. This produces a fast algorithm which is
    /// still robust.

    /// On each iteration Brent’s method approximates the function using an interpolating
    /// curve. On the first iteration this is a linear interpolation of the two endpoints. For
    /// subsequent iterations the algorithm uses an inverse quadratic fit to the last three
    /// points, for higher accuracy. The intercept of the interpolating curve with the x-axis
    /// is taken as a guess for the root. If it lies within the bounds of the current interval
    /// then the interpolating point is accepted, and used to generate a smaller interval. If
    /// the interpolating point is not accepted then the algorithm falls back to an ordinary
    /// bisection step.
    ///
    /// The best estimate of the root is taken from the most recent interpolation or bisection.
    #[doc(alias = "gsl_root_fsolver_falsepos")]
    pub fn falsepos() -> RootFSolverType {
        ffi_wrap!(gsl_root_fsolver_falsepos)
    }
}

ffi_wrapper!(
    RootFSolver<'a>,
    *mut sys::gsl_root_fsolver,
    gsl_root_fsolver_free
    ;inner_call: sys::gsl_function_struct => sys::gsl_function_struct { function: None, params: std::ptr::null_mut() };
    ;inner_closure: Option<Box<dyn Fn(f64) -> f64 + 'a>> => None;,
    "This is a workspace for finding roots using methods which do not require derivatives."
);

impl<'a> RootFSolver<'a> {
    /// This function returns a pointer to a newly allocated instance of a solver of type T.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    #[doc(alias = "gsl_root_fsolver_alloc")]
    pub fn new(t: RootFSolverType) -> Option<RootFSolver<'a>> {
        let tmp = unsafe { sys::gsl_root_fsolver_alloc(t.unwrap_shared()) };

        if tmp.is_null() {
            None
        } else {
            Some(RootFSolver::wrap(tmp))
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function f and
    /// the initial search interval [x lower, x upper].
    #[doc(alias = "gsl_root_fsolver_set")]
    pub fn set<F: Fn(f64) -> f64 + 'a>(
        &mut self,
        f: F,
        x_lower: f64,
        x_upper: f64,
    ) -> Result<(), Value> {
        self.inner_call = wrap_callback!(f, F + 'a);
        self.inner_closure = Some(Box::new(f));

        let ret = unsafe {
            sys::gsl_root_fsolver_set(self.unwrap_unique(), &mut self.inner_call, x_lower, x_upper)
        };
        result_handler!(ret, ())
    }

    /// The following function drives the iteration of each algorithm. Each function performs one
    /// iteration to update the state of any solver of the corresponding type. The same func-
    /// tion works for all solvers so that different methods can be substituted at runtime without
    /// modifications to the code.
    ///
    /// This function performs a single iteration of the solver s. If the iteration encounters
    /// an unexpected problem then an error code will be returned.
    ///
    /// The solver maintains a current best estimate of the root at all times. The bracketing
    /// solvers also keep track of the current best interval bounding the root.
    #[doc(alias = "gsl_root_fsolver_iterate")]
    pub fn iterate(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_root_fsolver_iterate(self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// Returns the solver type name.
    #[doc(alias = "gsl_root_fsolver_name")]
    pub fn name(&self) -> String {
        unsafe {
            let tmp = sys::gsl_root_fsolver_name(self.unwrap_shared());

            String::from_utf8_lossy(std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
        }
    }

    /// This function returns the current estimate of the root for the solver s.
    #[doc(alias = "gsl_root_fsolver_root")]
    pub fn root(&self) -> f64 {
        unsafe { sys::gsl_root_fsolver_root(self.unwrap_shared()) }
    }

    /// These functions return the current bracketing interval for the solver s.
    #[doc(alias = "gsl_root_fsolver_x_lower")]
    pub fn x_lower(&self) -> f64 {
        unsafe { sys::gsl_root_fsolver_x_lower(self.unwrap_shared()) }
    }

    /// These functions return the current bracketing interval for the solver s.
    #[doc(alias = "gsl_root_fsolver_x_upper")]
    pub fn x_upper(&self) -> f64 {
        unsafe { sys::gsl_root_fsolver_x_upper(self.unwrap_shared()) }
    }
}

ffi_wrapper!(
    RootFdfSolverType,
    *const sys::gsl_root_fdfsolver_type,
    "The root polishing algorithms described in this section require an initial guess for the
location of the root. There is no absolute guarantee of convergence—the function must be
suitable for this technique and the initial guess must be sufficiently close to the root
for it to work. When these conditions are satisfied then convergence is quadratic.
These algorithms make use of both the function and its derivative."
);

impl RootFdfSolverType {
    /// Newton’s Method is the standard root-polishing algorithm. The algorithm begins
    /// with an initial guess for the location of the root. On each iteration, a line tangent to
    /// the function f is drawn at that position. The point where this line crosses the x-axis
    /// becomes the new guess.
    #[doc(alias = "gsl_root_fdfsolver_newton")]
    pub fn newton() -> RootFdfSolverType {
        ffi_wrap!(gsl_root_fdfsolver_newton)
    }

    /// The secant method is a simplified version of Newton’s method which does not require
    /// the computation of the derivative on every step.
    #[doc(alias = "gsl_root_fdfsolver_secant")]
    pub fn secant() -> RootFdfSolverType {
        ffi_wrap!(gsl_root_fdfsolver_secant)
    }

    /// The Steffenson Method 1 provides the fastest convergence of all the routines. It com-
    /// bines the basic Newton algorithm with an Aitken “delta-squared” acceleration.
    #[doc(alias = "gsl_root_fdfsolver_steffenson")]
    pub fn steffenson() -> RootFdfSolverType {
        ffi_wrap!(gsl_root_fdfsolver_steffenson)
    }
}

ffi_wrapper!(
    RootFdfSolver<'a>,
    *mut sys::gsl_root_fdfsolver,
    gsl_root_fdfsolver_free
    ;inner_call: sys::gsl_function_fdf_struct => sys::gsl_function_fdf_struct{f: None, df: None, fdf: None, params: std::ptr::null_mut()};
    ;inner_f_closure: Option<Box<dyn Fn(f64) -> f64 + 'a>> => None;
    ;inner_df_closure: Option<Box<dyn Fn(f64) -> f64 + 'a>> => None;
    ;inner_fdf_closure: Option<Box<dyn Fn(f64, &mut f64, &mut f64) + 'a>> => None;,
    "This is a workspace for finding roots using methods which do require derivatives."
);

impl<'a> RootFdfSolver<'a> {
    /// This function returns a pointer to a newly allocated instance of a derivative-based
    /// solver of type T.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    #[doc(alias = "gsl_root_fdfsolver_alloc")]
    pub fn new(t: RootFdfSolverType) -> Option<RootFdfSolver<'a>> {
        let tmp = unsafe { sys::gsl_root_fdfsolver_alloc(t.unwrap_shared()) };

        if tmp.is_null() {
            None
        } else {
            Some(RootFdfSolver::wrap(tmp))
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function and
    /// derivative fdf and the initial guess root.
    #[doc(alias = "gsl_root_fdfsolver_set")]
    pub fn set<
        F: Fn(f64) -> f64 + 'a,
        DF: Fn(f64) -> f64 + 'a,
        FDF: Fn(f64, &mut f64, &mut f64) + 'a,
    >(
        &mut self,
        f: F,
        df: DF,
        fdf: FDF,
        root: f64,
    ) -> Result<(), Value> {
        // convert rust functions to C
        unsafe extern "C" fn inner_f<'a, F: Fn(f64) -> f64 + 'a>(
            x: c_double,
            params: *mut c_void,
        ) -> f64 {
            let f: &F = &*(params as *const F);
            f(x)
        }

        unsafe extern "C" fn inner_df<'a, DF: Fn(f64) -> f64 + 'a>(
            x: c_double,
            params: *mut c_void,
        ) -> f64 {
            let df: &DF = &*(params as *const DF);
            df(x)
        }

        unsafe extern "C" fn inner_fdf<'a, FDF: Fn(f64, &mut f64, &mut f64) + 'a>(
            x: c_double,
            params: *mut c_void,
            y: *mut c_double,
            dy: *mut c_double,
        ) {
            let fdf: &FDF = &*(params as *const FDF);
            fdf(x, &mut *y, &mut *dy);
        }

        self.inner_call = sys::gsl_function_fdf {
            f: Some(inner_f::<F>),
            df: Some(inner_df::<DF>),
            fdf: Some(inner_fdf::<FDF>),
            params: &(&f, &df, &fdf) as *const _ as *mut _,
        };
        self.inner_f_closure = Some(Box::new(f));
        self.inner_df_closure = Some(Box::new(df));
        self.inner_fdf_closure = Some(Box::new(fdf));

        let ret = unsafe {
            sys::gsl_root_fdfsolver_set(self.unwrap_unique(), &mut self.inner_call, root)
        };
        result_handler!(ret, ())
    }

    /// The following function drives the iteration of each algorithm. Each function performs one
    /// iteration to update the state of any solver of the corresponding type. The same func-
    /// tion works for all solvers so that different methods can be substituted at runtime without
    /// modifications to the code.
    ///
    /// This function performs a single iteration of the solver s. If the iteration encounters
    /// an unexpected problem then an error code will be returned.
    ///
    /// The solver maintains a current best estimate of the root at all times. The bracketing
    /// solvers also keep track of the current best interval bounding the root.
    #[doc(alias = "gsl_root_fdfsolver_iterate")]
    pub fn iterate(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_root_fdfsolver_iterate(self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// Returns the solver type name.
    #[doc(alias = "gsl_root_fdfsolver_name")]
    pub fn name(&self) -> Option<&str> {
        unsafe {
            let ptr = sys::gsl_root_fdfsolver_name(self.unwrap_shared());

            if ptr.is_null() {
                return None;
            }

            let mut len = 0;
            while *ptr.add(len) != 0 {
                len += 1;
            }

            let slice = std::slice::from_raw_parts(ptr as *const _, len);
            std::str::from_utf8(slice).ok()
        }
    }

    /// This function returns the current estimate of the root for the solver s.
    #[doc(alias = "gsl_root_fdfsolver_root")]
    pub fn root(&self) -> f64 {
        unsafe { sys::gsl_root_fdfsolver_root(self.unwrap_shared()) }
    }
}

#[cfg(any(test, doctest))]
mod test {
    /// This doc block will be used to ensure that the closure can't be set everywhere!
    ///
    /// ```compile_fail
    /// use crate::rgsl::*;
    ///
    /// fn set(root: &mut RootFSolver) {
    ///     let y = "lalal".to_owned();
    ///     root.set(
    ///         {|x| x * x - 5.;
    ///         println!("==> {:?}", y);
    ///         }, 0.0, 5.0);
    /// }
    ///
    /// let mut root = RootFSolver::new(RootFSolverType::brent()).unwrap();
    /// set(&mut root);
    /// ```
    ///
    /// Same but a working version:
    ///
    /// ```
    /// use crate::rgsl::*;
    ///
    /// fn set(root: &mut RootFSolver) {
    ///     root.set(|x| x * x - 5., 0.0, 5.0);
    /// }
    ///
    /// let mut root = RootFSolver::new(RootFSolverType::brent()).unwrap();
    /// set(&mut root);
    /// let status = root.iterate();
    /// ```
    use super::*;
    use crate::roots::{test_delta, test_interval};

    // support functions
    fn quadratic_test_fn(x: f64) -> f64 {
        x.powf(2.0) - 5.0
    }

    fn quadratic_test_fn_df(x: f64) -> f64 {
        2.0 * x
    }

    fn quadratic_test_fn_fdf(x: f64, y: &mut f64, dy: &mut f64) {
        *y = x.powf(2.0) - 5.0;
        *dy = 2.0 * x;
    }

    #[test]
    fn test_root() {
        let mut root = RootFSolver::new(RootFSolverType::brent()).unwrap();
        root.set(quadratic_test_fn, 0.0, 5.0).unwrap();

        let max_iter = 10usize;
        let epsabs = 0.0001;
        let epsrel = 0.0000001;

        let mut status = Value::Continue;
        let mut iter = 0usize;

        println!("Testing: {}", root.name());

        println!("iter, \t [x_lo, x_hi], \t min, \t error");
        while matches!(status, Value::Continue) && iter < max_iter {
            root.iterate().unwrap();

            // test for convergence
            let r = root.root();
            let x_lo = root.x_lower();
            let x_hi = root.x_upper();

            status = test_interval(x_lo, x_hi, epsabs, epsrel);

            // check if iteration succeeded
            if status == Value::Success {
                println!("Converged");
            }

            println!(
                "{} \t [{:.5}, {:.5}] \t {:.5} \t {:.5}",
                iter,
                x_lo,
                x_hi,
                r,
                x_hi - x_lo
            );
            iter += 1;
        }
        assert!(matches!(status, Value::Success))
    }

    #[test]
    fn test_root_fdf() {
        //guess value
        let r_expected = 5.0_f64.sqrt();
        let guess_value = 1.0;

        // setup solver
        let mut root = RootFdfSolver::new(RootFdfSolverType::steffenson()).unwrap();
        root.set(
            quadratic_test_fn,
            quadratic_test_fn_df,
            quadratic_test_fn_fdf,
            guess_value,
        )
        .unwrap();

        // set up iterations
        let max_iter = 20usize;
        let epsabs = 0.0001;
        let epsrel = 0.0000001;

        let mut status = Value::Continue;
        let mut iter = 0usize;

        println!("Testing: {}", root.name().unwrap());

        println!("iter, \t root, \t rel error \t abs error");

        let mut x = guess_value;
        while matches!(status, Value::Continue) && iter < max_iter {
            root.iterate().unwrap();

            // test for convergence
            let x_0 = x;
            x = root.root();
            // check if iteration succeeded
            status = test_delta(x, x_0, epsabs, epsrel);

            if matches!(status, Value::Success) {
                println!("Converged");
            }

            // print results
            println!(
                "{} \t {:.5} \t {:.5} \t {:.5}",
                iter,
                x,
                x - x_0,
                x - r_expected
            );
            iter += 1;
        }
        assert!(matches!(status, Value::Success))
    }
}
