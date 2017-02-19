//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Nonlinear Least-Squares Fitting

This chapter describes functions for multidimensional nonlinear least-squares fitting. The library
provides low level components for a variety of iterative solvers and convergence tests. These can
be combined by the user to achieve the desired solution, with full access to the intermediate steps
of the iteration. Each class of methods uses the same framework, so that you can switch between
solvers at runtime without needing to recompile your program. Each instance of a solver keeps track
of its own state, allowing the solvers to be used in multi-threaded programs.

##Overview

The problem of multidimensional nonlinear least-squares fitting requires the minimization of the
squared residuals of n functions, f_i, in p parameters, x_i,

\Phi(x) = (1/2) || F(x) ||^2
        = (1/2) \sum_{i=1}^{n} f_i(x_1, ..., x_p)^2
All algorithms proceed from an initial guess using the linearization,

\psi(p) = || F(x+p) || ~=~ || F(x) + J p ||
where x is the initial point, p is the proposed step and J is the Jacobian matrix J_{ij} = d f_i /
d x_j. Additional strategies are used to enlarge the region of convergence. These include requiring
a decrease in the norm ||F|| on each step or using a trust region to avoid steps which fall outside
the linear regime.

To perform a weighted least-squares fit of a nonlinear model Y(x,t) to data (t_i, y_i) with
independent Gaussian errors \sigma_i, use function components of the following form,

f_i = (Y(x, t_i) - y_i) / \sigma_i
Note that the model parameters are denoted by x in this chapter since the non-linear least-squares
algorithms are described geometrically (i.e. finding the minimum of a surface). The independent
variable of any data to be fitted is denoted by t.

With the definition above the Jacobian is J_{ij} =(1 / \sigma_i) d Y_i / d x_j, where Y_i =
Y(x,t_i).

##High Level Driver

These routines provide a high level wrapper that combine the iteration and convergence testing for
easy use.
*/

use ffi;
use libc::c_void;
use VectorF64;

pub struct MultiFitFSolverType {
    s: *mut ffi::gsl_multifit_fsolver_type,
}

impl ffi::FFI<ffi::gsl_multifit_fsolver_type> for MultiFitFSolverType {
    fn wrap(r: *mut ffi::gsl_multifit_fsolver_type) -> MultiFitFSolverType {
        MultiFitFSolverType {
            s: r,
        }
    }

    fn soft_wrap(r: *mut ffi::gsl_multifit_fsolver_type) -> MultiFitFSolverType {
        Self::wrap(r)
    }

    fn unwrap(s: &MultiFitFSolverType) -> *mut ffi::gsl_multifit_fsolver_type {
        s.s
    }
}

pub struct MultiFitFSolver {
    s: *mut ffi::gsl_multifit_fsolver,
}

impl MultiFitFSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n
    /// observations and p parameters. The number of observations n must be greater than or equal to
    /// parameters p.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    pub fn new(t: &MultiFitFSolverType, n: usize, p: usize) -> Option<MultiFitFSolver> {
        let tmp = unsafe { ffi::gsl_multifit_fsolver_alloc(ffi::FFI::unwrap(t), n, p) };

        if tmp.is_null() {
            None
        } else {
            Some(MultiFitFSolver {
                s: tmp,
            })
        }
    }

    pub fn set(&self, f: &mut MultiFitFunction, x: &VectorF64) -> ::Value {
        unsafe { ffi::gsl_multifit_fsolver_set(self.s, f, ffi::FFI::unwrap(x)) }
    }

    pub fn iterate(&self) -> ::Value {
        unsafe { ffi::gsl_multifit_fsolver_iterate(self.s) }
    }

    pub fn name(&self) -> String {
        unsafe {
            let tmp = ffi::gsl_multifit_fsolver_name(self.s);

            String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
        }
    }

    pub fn position(&self) -> VectorF64 {
        unsafe { ffi::FFI::wrap(ffi::gsl_multifit_fsolver_position(self.s)) }
    }
}

impl Drop for MultiFitFSolver {
    fn drop(&mut self) {
        unsafe { ffi::gsl_multifit_fsolver_free(self.s) };
        self.s = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_multifit_fsolver> for MultiFitFSolver {
    fn wrap(s: *mut ffi::gsl_multifit_fsolver) -> MultiFitFSolver {
        MultiFitFSolver {
            s: s
        }
    }

    fn soft_wrap(s: *mut ffi::gsl_multifit_fsolver) -> MultiFitFSolver {
        Self::wrap(s)
    }

    fn unwrap(s: &MultiFitFSolver) -> *mut ffi::gsl_multifit_fsolver {
        s.s
    }
}

#[repr(C)]
pub struct MultiFitFunction {
    pub f: Option<extern "C" fn(x: *const ffi::gsl_vector, params: *mut c_void,
                                f: *mut ffi::gsl_vector) -> ::Value>,
    /// number of functions
    pub n: usize,
    /// number of independent variables
    pub p: usize,
    pub params: *mut c_void,
}

pub struct MultiFitFdfSolver {
    intern: *mut ffi::gsl_multifit_fdfsolver,
}

impl MultiFitFdfSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n
    /// observations and p parameters. The number of observations n must be greater than or equal
    /// to parameters p.
    pub fn new(_type: &MultiFitFdfSolverType, n: usize, p: usize) -> Option<MultiFitFdfSolver> {
        let s = unsafe {
            ffi::gsl_multifit_fdfsolver_alloc(
                _type.intern as *const ffi::gsl_multifit_fdfsolver_type, n, p)
        };
        if s.is_null() {
            None
        } else {
            Some(MultiFitFdfSolver {
                intern: s,
            })
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function f and
    /// the initial guess x.
    pub fn set(&mut self, f: &mut MultiFitFunctionFdf, x: &::VectorF64) -> ::Value {
        unsafe { ffi::gsl_multifit_fdfsolver_set(self.intern, f.to_raw(), ffi::FFI::unwrap(x)) }
    }

    pub fn x(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.intern).x) }
    }

    pub fn f(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.intern).f) }
    }

    pub fn J(&self) -> ::MatrixF64 {
        unsafe { ffi::FFI::soft_wrap((*self.intern).J) }
    }

    pub fn dx(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.intern).dx) }
    }

    pub fn name(&self) -> String {
        unsafe {
            let tmp = ffi::gsl_multifit_fdfsolver_name(self.intern);

            String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
        }
    }

    /// This function performs a single iteration of the solver s. If the iteration encounters an
    /// unexpected problem then an error code will be returned. The solver maintains a current
    /// estimate of the best-fit parameters at all times.
    pub fn iterate(&mut self) -> ::Value {
        unsafe { ffi::gsl_multifit_fdfsolver_iterate(self.intern) }
    }

    /// This function returns the current position (i.e. best-fit parameters) s->x of the solver s.
    pub fn position(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::wrap(ffi::gsl_multifit_fdfsolver_position(self.intern)) }
    }

    /// These functions iterate the solver s for a maximum of maxiter iterations. After each
    /// iteration, the system is tested for convergence using gsl_multifit_test_delta with the
    /// error tolerances epsabs and epsrel.
    #[allow(unused_assignments)]
    pub fn driver(&mut self, max_iter: usize, epsabs: f64, epsrel: f64) -> ::Value {
        let mut status = ::Value::Failure;

        if !self.intern.is_null() {
            let mut iter = 0usize;
            loop {
                status = self.iterate();

                if status != ::Value::Success {
                    break
                }

                /* test for convergence */
                status = unsafe { ffi::gsl_multifit_test_delta((*self.intern).dx, (*self.intern).x,
                                                               epsabs, epsrel) };
                iter += 1;
                if status != ::Value::Continue || iter >= max_iter {
                    break
                }
            }
        }

        status
    }
}

impl Drop for MultiFitFdfSolver {
    fn drop(&mut self) {
        if !self.intern.is_null() {
            unsafe { ffi::gsl_multifit_fdfsolver_free(self.intern); }
            self.intern = ::std::ptr::null_mut();
        }
    }
}

#[allow(dead_code)]
pub struct MultiFitFdfSolverType {
    intern: *mut ffi::gsl_multifit_fdfsolver_type,
}

impl MultiFitFdfSolverType {
    pub fn lmder() -> MultiFitFdfSolverType {
        MultiFitFdfSolverType {
            intern: unsafe { ffi::gsl_multifit_fdfsolver_lmder },
        }
    }

    pub fn lmsder() -> MultiFitFdfSolverType {
        MultiFitFdfSolverType {
            intern: unsafe { ffi::gsl_multifit_fdfsolver_lmsder },
        }
    }
}

pub struct MultiFitFunctionFdf {
    pub f: Option<Box<Fn(::VectorF64, ::VectorF64) -> ::Value>>,
    pub df: Option<Box<Fn(::VectorF64, ::MatrixF64) -> ::Value>>,
    pub fdf: Option<Box<Fn(::VectorF64, ::VectorF64, ::MatrixF64) -> ::Value>>,
    pub n: usize,
    pub p: usize,
    intern: ffi::gsl_multifit_function_fdf,
}

impl MultiFitFunctionFdf {
    pub fn new(n: usize, p: usize) -> MultiFitFunctionFdf {
        MultiFitFunctionFdf {
            f: None,
            df: None,
            fdf: None,
            n: n,
            p: p,
            intern: ffi::gsl_multifit_function_fdf {
                f: Some(f),
                df: Some(df),
                fdf: Some(fdf),
                n: n,
                p: p,
                params: ::std::ptr::null_mut(),
            },
        }
    }

    fn to_raw(&mut self) -> *mut ffi::gsl_multifit_function_fdf {
        self.intern.n = self.n;
        self.intern.p = self.p;
        self.intern.params = self as *mut MultiFitFunctionFdf as *mut c_void;
        &mut self.intern
    }
}

extern "C" fn f(x: *mut ffi::gsl_vector, params: *mut c_void,
                pf: *mut ffi::gsl_vector) -> ::Value {
    unsafe {
        let t = params as *mut MultiFitFunctionFdf;
        if let Some(ref f) = (*t).f {
            f(ffi::FFI::soft_wrap(x), ffi::FFI::soft_wrap(pf))
        } else {
            ::Value::Success
        }
    }
}

extern "C" fn df(x: *mut ffi::gsl_vector, params: *mut c_void,
                 pdf: *mut ffi::gsl_matrix) -> ::Value {
    unsafe {
        let t = params as *mut MultiFitFunctionFdf;
        if let Some(ref df) = (*t).df {
            df(ffi::FFI::soft_wrap(x), ffi::FFI::soft_wrap(pdf))
        } else {
            ::Value::Success
        }
    }
}

extern "C" fn fdf(x: *mut ffi::gsl_vector, params: *mut c_void, pf: *mut ffi::gsl_vector,
                  pdf: *mut ffi::gsl_matrix) -> ::Value {
    unsafe {
        let t = params as *mut MultiFitFunctionFdf;
        if let Some(ref fdf) = (*t).fdf {
            fdf(ffi::FFI::soft_wrap(x), ffi::FFI::soft_wrap(pf), ffi::FFI::soft_wrap(pdf))
        } else {
            ::Value::Success
        }
    }
}
