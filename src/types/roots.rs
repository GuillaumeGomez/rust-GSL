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

use ffi;

/// The root bracketing algorithms described in this section require an initial interval which is
/// guaranteed to contain a root—if a and b are the endpoints of the interval then f (a) must
/// differ in sign from f (b). This ensures that the function crosses zero at least once in the
/// interval. If a valid initial interval is used then these algorithm cannot fail, provided the
/// function is well-behaved.
///
/// Note that a bracketing algorithm cannot find roots of even degree, since these do not
/// cross the x-axis.
pub struct RootFSolverType {
    s: *mut ffi::gsl_root_fsolver_type,
}

impl ffi::FFI<ffi::gsl_root_fsolver_type> for RootFSolverType {
    fn wrap(r: *mut ffi::gsl_root_fsolver_type) -> RootFSolverType {
        RootFSolverType { s: r }
    }

    fn soft_wrap(r: *mut ffi::gsl_root_fsolver_type) -> RootFSolverType {
        Self::wrap(r)
    }

    fn unwrap_shared(s: &RootFSolverType) -> *const ffi::gsl_root_fsolver_type {
        s.s as *const _
    }

    fn unwrap_unique(s: &mut RootFSolverType) -> *mut ffi::gsl_root_fsolver_type {
        s.s
    }
}

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
    pub fn bisection() -> RootFSolverType {
        RootFSolverType { s: unsafe { ffi::gsl_root_fsolver_bisection } }
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
    pub fn brent() -> RootFSolverType {
        RootFSolverType { s: unsafe { ffi::gsl_root_fsolver_brent } }
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
    pub fn falsepos() -> RootFSolverType {
        RootFSolverType { s: unsafe { ffi::gsl_root_fsolver_falsepos } }
    }
}

pub use ffi::gsl_function as RootFunction;

pub struct RootFSolver {
    s: *mut ffi::gsl_root_fsolver,
}


impl RootFSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    pub fn new(t: &RootFSolverType) -> Option<RootFSolver> {
        let tmp = unsafe { ffi::gsl_root_fsolver_alloc(ffi::FFI::unwrap_shared(t)) };

        if tmp.is_null() {
            None
        } else {
            Some(RootFSolver { s: tmp })
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function f and
    /// the initial search interval [x lower, x upper].
    pub fn set(&mut self, f: &mut RootFunction, x_lower: f64, x_upper: f64) -> ::Value {
        ::Value::from(unsafe {
            ffi::gsl_root_fsolver_set(self.s, f as *mut RootFunction, x_lower, x_upper)
        })
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
    pub fn iterate(&mut self) -> ::Value {
        ::Value::from(unsafe { ffi::gsl_root_fsolver_iterate(self.s) })
    }

    /// Returns the solver type name.
    pub fn name(&self) -> String {
        unsafe {
            let tmp = ffi::gsl_root_fsolver_name(self.s);

            String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
        }
    }

    /// This function returns the current estimate of the root for the solver s.
    pub fn root(&self) -> f64 {
        unsafe { ffi::gsl_root_fsolver_root(self.s) }
    }

    /// These functions return the current bracketing interval for the solver s.
    pub fn x_lower(&self) -> f64 {
        unsafe { ffi::gsl_root_fsolver_x_lower(self.s) }
    }

    /// These functions return the current bracketing interval for the solver s.
    pub fn x_upper(&self) -> f64 {
        unsafe { ffi::gsl_root_fsolver_x_upper(self.s) }
    }
}

impl Drop for RootFSolver {
    fn drop(&mut self) {
        if !self.s.is_null() {
            unsafe {
                ffi::gsl_root_fsolver_free(self.s);
            }
            self.s = ::std::ptr::null_mut();
        }
    }
}

/// The root polishing algorithms described in this section require an initial guess for the
/// location of the root. There is no absolute guarantee of convergence—the function must be
/// suitable for this technique and the initial guess must be sufficiently close to the root
/// for it to work. When these conditions are satisfied then convergence is quadratic.
///
/// These algorithms make use of both the function and its derivative.
pub struct RootFdfSolverType {
    s: *mut ffi::gsl_root_fdfsolver_type,
}

impl ffi::FFI<ffi::gsl_root_fdfsolver_type> for RootFdfSolverType {
    fn wrap(r: *mut ffi::gsl_root_fdfsolver_type) -> RootFdfSolverType {
        RootFdfSolverType { s: r }
    }

    fn soft_wrap(r: *mut ffi::gsl_root_fdfsolver_type) -> RootFdfSolverType {
        Self::wrap(r)
    }

    fn unwrap_shared(s: &RootFdfSolverType) -> *const ffi::gsl_root_fdfsolver_type {
        s.s as *const _
    }

    fn unwrap_unique(s: &mut RootFdfSolverType) -> *mut ffi::gsl_root_fdfsolver_type {
        s.s
    }
}

impl RootFdfSolverType {
    /// Newton’s Method is the standard root-polishing algorithm. The algorithm begins
    /// with an initial guess for the location of the root. On each iteration, a line tangent to
    /// the function f is drawn at that position. The point where this line crosses the x-axis
    /// becomes the new guess.
    pub fn newton() -> RootFdfSolverType {
        RootFdfSolverType { s: unsafe { ffi::gsl_root_fdfsolver_newton } }
    }

    /// The secant method is a simplified version of Newton’s method which does not require
    /// the computation of the derivative on every step.
    pub fn secant() -> RootFdfSolverType {
        RootFdfSolverType { s: unsafe { ffi::gsl_root_fdfsolver_secant } }
    }

    /// The Steffenson Method 1 provides the fastest convergence of all the routines. It com-
    /// bines the basic Newton algorithm with an Aitken “delta-squared” acceleration.
    pub fn steffenson() -> RootFdfSolverType {
        RootFdfSolverType { s: unsafe { ffi::gsl_root_fdfsolver_steffenson } }
    }
}

pub use ffi::gsl_function_fdf as RootFunctionFdf;

pub struct RootFdfSolver {
    s: *mut ffi::gsl_root_fdfsolver,
}

impl RootFdfSolver {
    /// This function returns a pointer to a newly allocated instance of a derivative-based
    /// solver of type T.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    pub fn new(t: &RootFdfSolverType) -> Option<RootFdfSolver> {
        let tmp = unsafe { ffi::gsl_root_fdfsolver_alloc(ffi::FFI::unwrap_shared(t)) };

        if tmp.is_null() {
            None
        } else {
            Some(RootFdfSolver { s: tmp })
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function and
    /// derivative fdf and the initial guess root.
    pub fn set(&mut self, f: &mut RootFunctionFdf, root: f64) -> ::Value {
        ::Value::from(unsafe {
            ffi::gsl_root_fdfsolver_set(self.s, f as *mut RootFunctionFdf, root)
        })
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
    pub fn iterate(&mut self) -> ::Value {
        ::Value::from(unsafe { ffi::gsl_root_fdfsolver_iterate(self.s) })
    }

    /// Returns the solver type name.
    pub fn name(&self) -> String {
        unsafe {
            let tmp = ffi::gsl_root_fdfsolver_name(self.s);

            String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
        }
    }

    /// This function returns the current estimate of the root for the solver s.
    pub fn root(&self) -> f64 {
        unsafe { ffi::gsl_root_fdfsolver_root(self.s) }
    }
}

impl Drop for RootFdfSolver {
    fn drop(&mut self) {
        if !self.s.is_null() {
            unsafe {
                ffi::gsl_root_fdfsolver_free(self.s);
            }
            self.s = ::std::ptr::null_mut();
        }
    }
}
