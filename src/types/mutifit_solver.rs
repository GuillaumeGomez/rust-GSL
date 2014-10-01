//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Nonlinear Least-Squares Fitting

This chapter describes functions for multidimensional nonlinear least-squares fitting. The library provides low level components for a 
variety of iterative solvers and convergence tests. These can be combined by the user to achieve the desired solution, with full access to 
the intermediate steps of the iteration. Each class of methods uses the same framework, so that you can switch between solvers at runtime 
without needing to recompile your program. Each instance of a solver keeps track of its own state, allowing the solvers to be used in 
multi-threaded programs.

##Overview

The problem of multidimensional nonlinear least-squares fitting requires the minimization of the squared residuals of n functions, f_i, in 
p parameters, x_i,

\Phi(x) = (1/2) || F(x) ||^2
        = (1/2) \sum_{i=1}^{n} f_i(x_1, ..., x_p)^2 
All algorithms proceed from an initial guess using the linearization,

\psi(p) = || F(x+p) || ~=~ || F(x) + J p ||
where x is the initial point, p is the proposed step and J is the Jacobian matrix J_{ij} = d f_i / d x_j. Additional strategies are used to 
enlarge the region of convergence. These include requiring a decrease in the norm ||F|| on each step or using a trust region to avoid steps 
which fall outside the linear regime.

To perform a weighted least-squares fit of a nonlinear model Y(x,t) to data (t_i, y_i) with independent Gaussian errors \sigma_i, use 
function components of the following form,

f_i = (Y(x, t_i) - y_i) / \sigma_i
Note that the model parameters are denoted by x in this chapter since the non-linear least-squares algorithms are described geometrically 
(i.e. finding the minimum of a surface). The independent variable of any data to be fitted is denoted by t.

With the definition above the Jacobian is J_{ij} =(1 / \sigma_i) d Y_i / d x_j, where Y_i = Y(x,t_i).

##High Level Driver

These routines provide a high level wrapper that combine the iteration and convergence testing for easy use.
!*/

/*use ffi;
use enums;

fn intern_fit_function(x: *const gsl_vector, params: *mut c_void, f: *mut gsl_vector) -> enums::Value {

}

pub struct MultiFitFunction<'r, T> {
    pub f: fn(x: &::VectorF64, params: &mut T, f: &::VectorF64) -> enums::Value,
    /// number of functions
    pub n: u64,
    /// number of independent variables
    pub p: u64,
    pub params: &'r mut T
}

pub struct MultiFitFdfSolver<T> {
    s: *mut ffi::gsl_multifit_fdfsolver,
    params: *mut c_void,
    f: Option<fn(x: &::VectorF64, params: &mut T, f: &::VectorF64) -> enums::Value>,
    c: ffi::gsl_multifit_function_fdf
}

impl<T> MultiFitFdfSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n observations and p parameters. The number
    /// of observations n must be greater than or equal to parameters p.
    /// 
    /// If there is insufficient memory to create the solver then the function returns a null pointer and the error handler is invoked with
    /// an error code of enums::NoMem.
    pub fn new(t: &MultiFitFdfSolverType, n: u64, p: u64) -> Option<MultiFitFdfSolver<T>> {
        let tmp = unsafe { ffi::gsl_multifit_fdfsolver_alloc(ffi::FFI::unwrap(t) as *const ffi::gsl_multifit_fdfsolver_type, n, p) };

        if tmp.is_null() {
            None
        } else {
            Some(MultiFitFdfSolver {
                s: tmp,
                params: None,
                f: None,
                c: ffi::gsl_multifit_fdfunction {
                    f: intern_fit_function,
                    n: 0,
                    p: 0,
                    params: *mut ::std::ptr::mut_null()
                }
            })
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function f and the initial guess x.
    pub fn set<U>(&mut self, f: &'r MultiFitFunction<U>, x: &::VectorF64) -> enums::Value {
        self.params = Some(::std::mem::transmute(f.params));
        self.f = Some(f.f);
        self.c.n = f.n;
        self.c.p = f.p;
        self.c.params = ::std::mem::transmute(self);
        unsafe { ffi::gsl_multifit_fdfsolver_set(self.s, &mut self.c, ffi::FFI::unwrap(x)) }
    }

    pub fn name(&self) -> Option<String> {
        let tmp = unsafe { ffi::gsl_multifit_fdfsolver_name(self.s as *const ffi::gsl_multifit_fdfsolver) };

        if tmp.is_null() {
            None
        } else {
            Some(::std::str::raw::from_buf(tmp))
        }
    }

    /// This function takes as input the current position x with the function values computed at the current position f, along with fdf
    /// which specifies the fit function and parameters and approximates the n-by-p Jacobian J using forward finite differences: J_ij = d
    /// f_i(x,params) / d x_j = (f_i(x^*,params) - f_i(x,params)) / d x_j. where x^* has the jth element perturbed by \Delta x_j and \Delta
    /// x_j = \epsilon |x_j|, where \epsilon is the square root of the machine precision ::DBL_EPSILON.
    //pub fn dif_df(&self, )

    /// This function performs a single iteration of the solver s. If the iteration encounters an unexpected problem then an error code
    /// will be returned. The solver maintains a current estimate of the best-fit parameters at all times.
    pub fn iterate(&self) -> enums::Value {
        unsafe { ffi::gsl_multifit_fdfsolver_iterate(self.s) }
    }

    /// This function returns the current position (i.e. best-fit parameters) s->x of the solver s.
    pub fn position(&self) -> ::VectorF64 {
        let tmp = unsafe { ffi::gsl_multifit_fdfsolver_position(self.s as *const ffi::gsl_multifit_fdfsolver) };

        ::types::vector::wrap(tmp)
    }

    /// These functions iterate the solver s for a maximum of maxiter iterations. After each iteration, the system is tested for convergence
    /// using gsl_multifit_test_delta with the error tolerances epsabs and epsrel.
    pub fn driver(&self, max_iter: u64, epsabs: f64, epsrel: f64) -> enums::Value {
        unsafe { ffi::gsl_multifit_fdfsolver_driver(self.s, max_iter, epsabs, epsrel) }
    }
}

impl<T> Drop for MultiFitFdfSolver<T> {
    fn drop(&mut self) {
        unsafe { ffi::gsl_multifit_fdfsolver_free(self.s) };
        self.s = ::std::ptr::null_mut();
    }
}*/

/*impl ffi::FFI<ffi::gsl_multifit_fdfsolver> for MultiFitFdfSolver {
    fn wrap(s: *mut ffi::gsl_multifit_fdfsolver) -> MultiFitFdfSolver {
        MultiFitFdfSolver {
            s: s
        }
    }

    fn unwrap(s: &MultiFitFdfSolver) -> *mut ffi::gsl_multifit_fdfsolver {
        s.s
    }
}*/