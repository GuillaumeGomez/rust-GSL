//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Nonlinear Least-Squares Fitting

This chapter describes functions for multidimensional nonlinear least-squares fitting. The library
provides low level components for a variety of iterative solvers and convergence tests. These can
be combined by the user to achieve the desired solution, with full access to the intermediate steps
of the iteration. Each class of methods uses the same framework, so that you can switch between
solvers at runtime without needing to recompile your program. Each instance of a solver keeps track
of its own state, allowing the solvers to be used in multi-threaded programs.

## Overview

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

## High Level Driver

These routines provide a high level wrapper that combine the iteration and convergence testing for
easy use.
*/

/*
C Equivalent code:

```
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

struct data {
    size_t n;
    double * y;
    double * sigma;
};

void print_state(size_t iter, gsl_multifit_fdfsolver * s) {
    printf("iter: %3u x = % 15.8f % 15.8f % 15.8f |f(x)| = %g\n",
           iter,
           gsl_vector_get (s->x, 0),
           gsl_vector_get (s->x, 1),
           gsl_vector_get (s->x, 2),
           gsl_blas_dnrm2 (s->f));
}

int expb_f(const gsl_vector * x, void *params,
           gsl_vector * f) {
    size_t n = ((struct data *)params)->n;
    double *y = ((struct data *)params)->y;
    double *sigma = ((struct data *) params)->sigma;

    double A = gsl_vector_get (x, 0);
    double lambda = gsl_vector_get (x, 1);
    double b = gsl_vector_get (x, 2);

    size_t i;

    for (i = 0; i < n; i++) {
        /* Model Yi = A * exp(-lambda * i) + b */
        double t = i;
        double Yi = A * exp (-lambda * t) + b;
        gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }

    return GSL_SUCCESS;
}

int expb_df(const gsl_vector * x, void *params,
            gsl_matrix * J) {
    size_t n = ((struct data *)params)->n;
    double *sigma = ((struct data *) params)->sigma;

    double A = gsl_vector_get (x, 0);
    double lambda = gsl_vector_get (x, 1);

    size_t i;

    for (i = 0; i < n; i++) {
        /* Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = (Yi - yi)/sigma[i],      */
        /*       Yi = A * exp(-lambda * i) + b  */
        /* and the xj are the parameters (A,lambda,b) */
        double t = i;
        double s = sigma[i];
        double e = exp(-lambda * t);
        gsl_matrix_set (J, i, 0, e/s);
        gsl_matrix_set (J, i, 1, -t * A * e/s);
        gsl_matrix_set (J, i, 2, 1/s);

    }
    return GSL_SUCCESS;
}

int expb_fdf(const gsl_vector * x, void *params,
             gsl_vector * f, gsl_matrix * J) {
    expb_f (x, params, f);
    expb_df (x, params, J);

    return GSL_SUCCESS;
}

int main(void) {
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = 40;
    const size_t p = 3;

    gsl_matrix *covar = gsl_matrix_alloc (p, p);

    double y[n], sigma[n];

    struct data d = { n, y, sigma};

    gsl_multifit_function_fdf f;

    double x_init[3] = { 1.0, 0.0, 0.0 };

    gsl_vector_view x = gsl_vector_view_array (x_init, p);

    const gsl_rng_type * type;
    gsl_rng * r;

    gsl_rng_env_setup();

    type = gsl_rng_default;
    r = gsl_rng_alloc (type);

    f.f = &expb_f;
    f.df = &expb_df;
    f.fdf = &expb_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    /* This is the data to be fitted */

    for (i = 0; i < n; i++) {
        double t = i;
        y[i] = 1.0 + 5 * exp (-0.1 * t) + gsl_ran_gaussian(r, 0.1);
        sigma[i] = 0.1;
        printf("data: %d %g %g\n", i, y[i], sigma[i]);
    }

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

    print_state (iter, s);

    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);

        printf ("status = %s\n", gsl_strerror (status));

        print_state (iter, s);

        if (status)
          break;

        status = gsl_multifit_test_delta (s->dx, s->x,
                          1e-4, 1e-4);
    } while (status == GSL_CONTINUE && iter < 500);

    gsl_multifit_covar (s->J, 0.0, covar);

    gsl_matrix_fprintf (stdout, covar, "%g");

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    printf("A      = %.5f +/- %.5f\n", FIT(0), ERR(0));
    printf("lambda = %.5f +/- %.5f\n", FIT(1), ERR(1));
    printf("b      = %.5f +/- %.5f\n", FIT(2), ERR(2));

    printf ("status = %s\n", gsl_strerror (status));

    gsl_multifit_fdfsolver_free (s);
    return 0;
}
```
*/

use ffi::{self, FFI};
use std::os::raw::{c_int, c_void};
use VectorF64;

ffi_wrapper!(MultiFitFSolverType, *mut sys::gsl_multifit_fsolver_type);

// pub struct MultiFitFunction<F: Fn(x: &::VectorF64, f: &mut ::VectorF64)> {
//     pub f: Box<F>,
//     /// Number of functions.
//     pub n: usize,
//     /// Number of independent variables.
//     pub p,
// }

pub struct MultiFitFunction(pub sys::gsl_multifit_function);

ffi_wrapper!(
    MultiFitFSolver,
    *mut sys::gsl_multifit_fsolver,
    gsl_multifit_fsolver_free
);

impl MultiFitFSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n
    /// observations and p parameters. The number of observations n must be greater than or equal to
    /// parameters p.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    #[doc(alias = "gsl_multifit_fsolver_alloc")]
    pub fn new(t: &MultiFitFSolverType, n: usize, p: usize) -> Option<MultiFitFSolver> {
        let tmp = unsafe { sys::gsl_multifit_fsolver_alloc(t.unwrap_shared(), n, p) };

        if tmp.is_null() {
            None
        } else {
            Some(MultiFitFSolver::wrap(tmp))
        }
    }

    #[doc(alias = "gsl_multifit_fsolver_set")]
    pub fn set(&mut self, f: &mut MultiFitFunction, x: &mut VectorF64) -> ::Value {
        // unsafe {
        //     let func = (*self.0).function;
        //     if !func.is_null() {
        //         Box::from_raw((*func).params);
        //     }
        // }
        ::Value::from(unsafe {
            sys::gsl_multifit_fsolver_set(self.unwrap_unique(), &mut f.0, x.unwrap_shared())
        })
    }

    #[doc(alias = "gsl_multifit_fsolver_iterate")]
    pub fn iterate(&mut self) -> ::Value {
        ::Value::from(unsafe { sys::gsl_multifit_fsolver_iterate(self.unwrap_unique()) })
    }

    #[doc(alias = "gsl_multifit_fsolver_name")]
    pub fn name(&self) -> String {
        unsafe {
            let tmp = sys::gsl_multifit_fsolver_name(self.unwrap_shared());

            String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
        }
    }

    #[doc(alias = "gsl_multifit_fsolver_position")]
    pub fn position(&self) -> VectorF64 {
        unsafe { ffi::FFI::wrap(sys::gsl_multifit_fsolver_position(self.unwrap_shared())) }
    }
}

ffi_wrapper!(
    MultiFitFdfSolver,
    *mut sys::gsl_multifit_fdfsolver,
    gsl_multifit_fdfsolver_free
);

impl MultiFitFdfSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n
    /// observations and p parameters. The number of observations n must be greater than or equal
    /// to parameters p.
    #[doc(alias = "gsl_multifit_fdfsolver_alloc")]
    pub fn new(_type: &MultiFitFdfSolverType, n: usize, p: usize) -> Option<MultiFitFdfSolver> {
        let s = unsafe { sys::gsl_multifit_fdfsolver_alloc(_type.unwrap_shared(), n, p) };
        if s.is_null() {
            None
        } else {
            Some(MultiFitFdfSolver::wrap(s))
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function f and
    /// the initial guess x.
    #[doc(alias = "gsl_multifit_fdfsolver_set")]
    pub fn set(&mut self, f: &mut MultiFitFunctionFdf, x: &::VectorF64) -> ::Value {
        ::Value::from(unsafe {
            sys::gsl_multifit_fdfsolver_set(self.unwrap_unique(), f.to_raw(), x.unwrap_shared())
        })
    }

    pub fn x(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.unwrap_shared()).x) }
    }

    pub fn f(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.unwrap_shared()).f) }
    }

    pub fn dx(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.unwrap_shared()).dx) }
    }

    pub fn g(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.unwrap_shared()).g) }
    }

    pub fn sqrt_wts(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::soft_wrap((*self.unwrap_shared()).sqrt_wts) }
    }

    #[doc(alias = "gsl_multifit_fdfsolver_name")]
    pub fn name(&self) -> String {
        unsafe {
            let tmp = sys::gsl_multifit_fdfsolver_name(self.unwrap_shared());

            String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
        }
    }

    /// This function performs a single iteration of the solver s. If the iteration encounters an
    /// unexpected problem then an error code will be returned. The solver maintains a current
    /// estimate of the best-fit parameters at all times.
    #[doc(alias = "gsl_multifit_fdfsolver_iterate")]
    pub fn iterate(&mut self) -> ::Value {
        ::Value::from(unsafe { sys::gsl_multifit_fdfsolver_iterate(self.unwrap_unique()) })
    }

    /// This function returns the current position (i.e. best-fit parameters) s->x of the solver s.
    #[doc(alias = "gsl_multifit_fdfsolver_position")]
    pub fn position(&self) -> ::VectorF64 {
        unsafe { ffi::FFI::wrap(sys::gsl_multifit_fdfsolver_position(self.unwrap_shared())) }
    }

    /// These functions iterate the solver s for a maximum of maxiter iterations. After each
    /// iteration, the system is tested for convergence using gsl_multifit_test_delta with the
    /// error tolerances epsabs and epsrel.
    #[allow(unused_assignments)]
    #[doc(alias = "gsl_multifit_test_delta")]
    pub fn driver(&mut self, max_iter: usize, epsabs: f64, epsrel: f64) -> ::Value {
        let mut status = ::Value::Failure;
        let ptr = self.unwrap_shared();

        if !ptr.is_null() {
            let mut iter = 0usize;
            loop {
                status = self.iterate();

                if status != ::Value::Success {
                    break;
                }

                /* test for convergence */
                status = ::Value::from(unsafe {
                    sys::gsl_multifit_test_delta((*ptr).dx, (*ptr).x, epsabs, epsrel)
                });
                iter += 1;
                if status != ::Value::Continue || iter >= max_iter {
                    break;
                }
            }
        }

        status
    }
}

ffi_wrapper!(
    MultiFitFdfSolverType,
    *const sys::gsl_multifit_fdfsolver_type
);

impl MultiFitFdfSolverType {
    pub fn lmder() -> MultiFitFdfSolverType {
        ffi_wrap!(gsl_multifit_fdfsolver_lmder)
    }

    pub fn lmsder() -> MultiFitFdfSolverType {
        ffi_wrap!(gsl_multifit_fdfsolver_lmsder)
    }
}

pub struct MultiFitFunctionFdf {
    pub f: Option<Box<dyn Fn(::VectorF64, ::VectorF64) -> ::Value>>,
    pub df: Option<Box<dyn Fn(::VectorF64, ::MatrixF64) -> ::Value>>,
    pub fdf: Option<Box<dyn Fn(::VectorF64, ::VectorF64, ::MatrixF64) -> ::Value>>,
    pub n: usize,
    pub p: usize,
    intern: sys::gsl_multifit_function_fdf,
}

impl MultiFitFunctionFdf {
    pub fn new(n: usize, p: usize, nevalf: usize, nevaldf: usize) -> MultiFitFunctionFdf {
        MultiFitFunctionFdf {
            f: None,
            df: None,
            fdf: None,
            n,
            p,
            intern: sys::gsl_multifit_function_fdf {
                f: Some(f),
                df: Some(df),
                fdf: Some(fdf),
                n,
                p,
                params: ::std::ptr::null_mut(),
                nevalf,
                nevaldf,
            },
        }
    }

    #[allow(clippy::wrong_self_convention)]
    fn to_raw(&mut self) -> *mut sys::gsl_multifit_function_fdf {
        self.intern.n = self.n;
        self.intern.p = self.p;
        self.intern.params = self as *mut MultiFitFunctionFdf as *mut c_void;
        &mut self.intern
    }
}

unsafe extern "C" fn f(
    x: *const sys::gsl_vector,
    params: *mut c_void,
    pf: *mut sys::gsl_vector,
) -> c_int {
    let t = params as *mut MultiFitFunctionFdf;
    if let Some(ref i_f) = (*t).f {
        i_f(
            ffi::FFI::soft_wrap(x as usize as *mut _),
            ffi::FFI::soft_wrap(pf),
        )
        .into()
    } else {
        ::Value::Success.into()
    }
}

unsafe extern "C" fn df(
    x: *const sys::gsl_vector,
    params: *mut c_void,
    pdf: *mut sys::gsl_matrix,
) -> c_int {
    let t = params as *mut MultiFitFunctionFdf;
    if let Some(ref i_df) = (*t).df {
        i_df(
            ffi::FFI::soft_wrap(x as usize as *mut _),
            ffi::FFI::soft_wrap(pdf),
        )
        .into()
    } else {
        ::Value::Success.into()
    }
}

unsafe extern "C" fn fdf(
    x: *const sys::gsl_vector,
    params: *mut c_void,
    pf: *mut sys::gsl_vector,
    pdf: *mut sys::gsl_matrix,
) -> c_int {
    let t = params as *mut MultiFitFunctionFdf;
    if let Some(ref i_fdf) = (*t).fdf {
        i_fdf(
            ffi::FFI::soft_wrap(x as usize as *mut _),
            ffi::FFI::soft_wrap(pf),
            ffi::FFI::soft_wrap(pdf),
        )
        .into()
    } else {
        ::Value::Success.into()
    }
}
