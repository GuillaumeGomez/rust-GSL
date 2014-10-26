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

use ffi;
use enums;

pub struct MultiFitFunction<'r, T> {
    pub f: fn(x: &::VectorF64, params: &mut T, f: &::VectorF64) -> enums::value::Value,
    /// number of functions
    pub n: u64,
    /// number of independent variables
    pub p: u64,
    pub params: &'r mut T
}

pub struct MultiFitFdfSolver<T> {
    s: *mut ffi::gsl_multifit_fdfsolver,
    params: *mut c_void,
    f: Option<fn(x: &::VectorF64, params: &mut T, f: &::VectorF64) -> enums::value::Value>,
    c: ffi::gsl_multifit_function_fdf
}

impl<T> MultiFitFdfSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n observations and p parameters. The number
    /// of observations n must be greater than or equal to parameters p.
    /// 
    /// If there is insufficient memory to create the solver then the function returns a null pointer and the error handler is invoked with
    /// an error code of ::NoMem.
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
    pub fn set<U>(&mut self, f: &'r MultiFitFunction<U>, x: &::VectorF64) -> enums::value::Value {
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
    pub fn iterate(&self) -> enums::value::Value {
        unsafe { ffi::gsl_multifit_fdfsolver_iterate(self.s) }
    }

    /// This function returns the current position (i.e. best-fit parameters) s->x of the solver s.
    pub fn position(&self) -> ::VectorF64 {
        let tmp = unsafe { ffi::gsl_multifit_fdfsolver_position(self.s as *const ffi::gsl_multifit_fdfsolver) };

        ::types::vector::wrap(tmp)
    }

    /// These functions iterate the solver s for a maximum of maxiter iterations. After each iteration, the system is tested for convergence
    /// using gsl_multifit_test_delta with the error tolerances epsabs and epsrel.
    pub fn driver(&self, max_iter: u64, epsabs: f64, epsrel: f64) -> enums::value::Value {
        unsafe { ffi::gsl_multifit_fdfsolver_driver(self.s, max_iter, epsabs, epsrel) }
    }
}

impl<T> Drop for MultiFitFdfSolver<T> {
    fn drop(&mut self) {
        unsafe { ffi::gsl_multifit_fdfsolver_free(self.s) };
        self.s = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_multifit_fdfsolver> for MultiFitFdfSolver {
    fn wrap(s: *mut ffi::gsl_multifit_fdfsolver) -> MultiFitFdfSolver {
        MultiFitFdfSolver {
            s: s
        }
    }

    fn unwrap(s: &MultiFitFdfSolver) -> *mut ffi::gsl_multifit_fdfsolver {
        s.s
    }
}

pub struct MultiFitFdfSolverType {
    name: &'static str,
    size: u64,
    alloc: fn(state: &mut LmderStateT, n: u64, p: u64) -> enums::value::Value,
    set: fn(state: &mut LmderStateT, fdf: MultiFitFunctionFdf, x: &mut VectorF64, f: &mut VectorF64, j: MatrixF64,
        dx: &mut VectorF64) -> enums::value::Value,
    iterate: fn(state: &mut LmderStateT, fdf: MultiFitFunctionFdf, x: &mut VectorF64, f: VectorF64, j: &mut MatrixF64,
        dx: &mut VectorF64) -> enums::value::Value,
    free: fn(state: &mut LmderStateT)
}

struct MultiFitFunctionFdf<T> {
    f: Option<fn<T>(x: &::VectorF64, params: &mut T, f: &mut ::VectorF64) -> enums::value::Value>,
    df: Option<fn<T>(x: &::VectorF64, params: &mut T, df: &mut ::MatrixF64) -> enums::value::Value>,
    fdf: Option<fn<T>(x: &::VectorF64, params: &mut T, f: &mut ::VectorF64, df: &mut ::MatrixF64) -> enums::value::Value>,
    n: u64,
    p: u64,
    params: &'r mut T
}

struct LmderStateT {
    iter: u64,
    xnorm: f64,
    fnorm: f64,
    delta: f64,
    par: f64,
    r: ::MatrixF64,
    tau: ::VectorF64,
    diag: ::VectorF64,
    qtf: ::VectorF64,
    newton: ::VectorF64,
    gradient: ::VectorF64,
    x_trial: ::VectorF64,
    f_trial: ::VectorF64,
    df: ::VectorF64,
    sdiag: ::VectorF64,
    rptdx: ::VectorF64,
    w: ::VectorF64,
    work1: ::VectorF64,
    perm: ::Permutation
}

fn lmder_alloc(state: &mut LmderStateT, size_t n, size_t p) -> enums::value::Value {
    state.r = match ::MatrixF64::new(n, p) {
        Some(m) => m,
        None => {
            rgsl_error!("failed to allocate space for r", enums::value::NoMem)
        }
    };

    state.tau = match ::VectorF64::new(n.min(p)) {
        Some(t) => t,
        None => {
            ::std::mem::drop(state.r);
            rgsl_error!("failed to allocate space for tau", enums::value::NoMem)
        }
    };

    state.diag = match ::VectorF64;:new(p) {
        Some(d) => d,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            rgsl_error!("failed to allocate space for diag", enums::value::NoMem)
        }
    };

    state.qtf = match ::VectorF64::new(n) {
        Some(q) => q,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            rgsl_error!("failed to allocate space for qtf", enums::value::NoMem)
        }
    };

    state.newton = match ::VectorF64::new(p) {
        Some(n) => n,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            rgsl_error!("failed to allocate space for newton", enums::value::NoMem)
        }
    };

    state.gradient = match ::VectorF64::new(p) {
        Some(g) => g,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            rgsl_error!("failed to allocate space for gradient", enums::value::NoMem)
        }
    };

    state.x_trial = match ::VectorF64::new(p) {
        Some(x) => x,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            rgsl_error!("failed to allocate space for x_trial", enums::value::NoMem)
        }
    };

    state.f_trial = ::VectorF64::new(n) {
        Some(f) => f,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            rgsl_error!("failed to allocate space for f_trial", enums::value::NoMem)
        }
    };

    state.df = ::VectorF64::new(n) {
        Some(d) => d,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            ::std::mem::drop(state.f_trial);
            rgsl_error!("failed to allocate space for df", enums::value::NoMem)
        }
    };

    state.sdiag = ::VectorF64::new(p) {
        Some(s) => s,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            ::std::mem::drop(state.f_trial);
            ::std::mem::drop(state.df);
            rgsl_error!("failed to allocate space for sdiag", enums::value::NoMem)
        }
    };

    state.rptdx = ::VectorF64::new(n) {
        Some(s) => s,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            ::std::mem::drop(state.f_trial);
            ::std::mem::drop(state.df);
            ::std::mem::drop(state.sdiag);
            rgsl_error!("failed to allocate space for rptdx", enums::value::NoMem)
        }
    };

    state.w = ::VectorF64::new(n) {
        Some(w) => w,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            ::std::mem::drop(state.f_trial);
            ::std::mem::drop(state.df);
            ::std::mem::drop(state.sdiag);
            ::std::mem::drop(state.rptdx);
            rgsl_error!("failed to allocate space for w", enums::value::NoMem)
        }
    };

    state.work1 = ::VectorF64::new(p) {
        Some(w) => w,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            ::std::mem::drop(state.f_trial);
            ::std::mem::drop(state.df);
            ::std::mem::drop(state.sdiag);
            ::std::mem::drop(state.rptdx);
            ::std::mem::drop(state.w);
            rgsl_error!("failed to allocate space for work1", enums::value::NoMem)
        }
    };

    state.perm = ::Permutation::new(p) {
        Some(p) => p,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            ::std::mem::drop(state.f_trial);
            ::std::mem::drop(state.df);
            ::std::mem::drop(state.sdiag);
            ::std::mem::drop(state.rptdx);
            ::std::mem::drop(state.w);
            ::std::mem::drop(state.work1);
            rgsl_error!("failed to allocate space for perm", enums::value::NoMem)
        }
    }

    return enums::value::Success;
}

fn lmder_set(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut VectorF64, f: &mut VectorF64, J: &mut MatrixF64,
    dx: &mut VectorF64) -> enums::value::Value {
    set(vstate, fdf, x, f, J, dx, 0)
}

fn lmsder_set(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut VectorF64, f: &mut VectorF64, J: &mut MatrixF64,
    dx: &mut VectorF64) -> enums::value::Value {
    set(vstate, fdf, x, f, J, dx, 1)
}

fn lmder_iterate(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut VectorF64, f: &mut VectorF64, J: &mut MatrixF64,
    dx: &mut VectorF64) -> enums::value::Value {
    iterate(vstate, fdf, x, f, J, dx, 0)
}

fn lmsder_iterate(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut VectorF64, f: &mut VectorF64, J: &mut MatrixF64,
    dx: &mut VectorF64) -> enums::value::Value {
    iterate(vstate, fdf, x, f, J, dx, 1)
}

fn compute_diag(j: &MatrixF64, diag: &mut VectorF64) {
    let n = j.size1();
    let p = j.size2();

    for j in range(0, p) {
        let mut sum = 0f64;

        for i in range(0, n) {
            let jij = j.get(i, j);

            sum += jij * jij;
        }
        if sum == 0f64 {
            sum = 0f64;
        }

        diag.set(j, sum.sqrt());
    }
}

fn update_diag(j: &::MatrixF64, diag: &mut ::VectorF64) {
    let n = diag.size();

    for j in range(0, n) {
        let mut sum = 0f64;

        for i in range(0, n) {
          let jij = j.get(i, j);
          
          sum += jij * jij;
        }
        if sum == 0f64 {
            sum = 1f64;
        }

        let cnorm = sum.sqrt();
        let diagj = diag.get(j);

        if cnorm > diagj {
            diag.set(j, cnorm);
        }
    }
}

fn scaled_enorm(d: &::VectorF64, f: &::VectorF64) -> f64 {
    let mut e2 = 0f64;
    let n = f.size();

    for i in range(0, n) {
        let fi = f.get(i);
        let di = d.get(i);
        let u = di * fi;

        e2 += u * u;
    }
    e2.sqrt()
}

fn compute_delta(diag: &mut ::VectorF64, x: &mut ::VectorF64) -> f64 {
    let Dx = scaled_enorm(diag, x);
    let factor = 100f64;  /* generally recommended value from MINPACK */

    if Dx > 0 {factor * Dx} else {factor}
}

fn set(state: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, t: &mut ::VectorF64, j: &mut ::MatrixF64,
    dx: &mut ::VectorF64, scale: i32) -> enums::value::Value {
    let mut signum = 0i32;

    /* Evaluate function at x */
    /* return immediately if evaluation raised error */
    {
        let status = if fdf.fdf.is_some() {
            fdf.fdf(x, fdf.params, f, j)
        } else {
            /* finite difference approximation */
            gsl_multifit_fdfsolver_dif_fdf(x, fdf, f, J)
        }

        if status != enums::value::Success
            return status;
    }

    state.par = 0;
    state.iter = 1;
    state.fnorm = enorm (f);

    dx.set_all(0f64);

    /* store column norms in diag */
    if scale != 0 {
        compute_diag(J, diag);
    } else {
        diag.set_all(1f64);
    }

    /* set delta to 100 |D x| or to 100 if |D x| is zero */
    state.xnorm = scaled_enorm(diag, x);
    state.delta = compute_delta(diag, x);

    /* Factorize J into QR decomposition */
    r.copy_from(j);
    ::linear_algebra::QRPT_decomp(r, tau, perm, &mut signum, work1);

    state.rptdx.set_zero();
    state.w.set_zero();

    /* Zero the trial vector, as in the alloc function */

    state.f_trials.set_zero();

    /*#ifdef DEBUG
    printf("r = "); gsl_matrix_fprintf(stdout, r, "%g");
    printf("perm = "); gsl_permutation_fprintf(stdout, perm, "%d");
    printf("tau = "); gsl_vector_fprintf(stdout, tau, "%g");
    #endif*/

    enums::value::Success
}

fn free() {
    ::std::mem::drop(state.r);
    ::std::mem::drop(state.tau);
    ::std::mem::drop(state.diag);
    ::std::mem::drop(state.qtf);
    ::std::mem::drop(state.newton);
    ::std::mem::drop(state.gradient);
    ::std::mem::drop(state.x_trial);
    ::std::mem::drop(state.f_trial);
    ::std::mem::drop(state.df);
    ::std::mem::drop(state.sdiag);
    ::std::mem::drop(state.rptdx);
    ::std::mem::drop(state.w);
    ::std::mem::drop(state.work1);
}

fn enorm(f: &::VectorF64) -> f64 {
    ::blas::level1::dnrm2(f)
}

fn compute_gradient_direction(r: &::MatrixF64, p: &::Permutation, qtf: &::VectorF64, diag: &::VectorF64, g: &mut ::VectorF64) {
    let n = r.size2();

    for j in range(0, n) {
        let mut sum = 0f64;

        for i in range(0, j + 1) {
            sum += r.get(i, j) * qtf.get(i);
        }

        {
            let pj = p.get(j);
            let dpj = diag.get(pj);

            g.set(j, sum / dpj);
        }
    }
}

fn compute_trial_step(x: &mut ::VectorF64, dx: &mut ::VectorF64, x_trial: &mut ::VectorF64) {
    let n = x.size();

    for i in range(0, n) {
        let pi = dx.get(i);
        let xi = x.get(i);

        x_trial.set(i, xi + pi);
    }
}

fn compute_actual_reduction(fnorm: f64, fnorm1: f64) -> f64 {
    if 0.1f64 * fnorm1 < fnorm {
        let u = fnorm1 / fnorm;

        1f64 - u * u
    } else {
        -1f64
    }
}

fn compute_rptdx(r: &::MatrixF64, p: &::Permutation, dx: &::VectorF64, rptdx: &mut ::VectorF64) {
    let n = dx.size;

    for i in range(0, n) {
        let mut sum = 0f64;

        for j in range(i, n) {
            let pj = p.get(j);

            sum += r.get(i, j) * dx.get(pj);
        }

        rptdx.set(i, sum);
    }
}

fn iterate(state: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut VectorF64, f: &mut VectorF64, J: &mut MatrixF64,
    dx: &mut VectorF64, scale: i32) -> enums::value::Value {
    let mut prered = 0f64;
    let mut actred = 0f64;
    let mut pnorm = 0f64;
    let mut fnorm1 = 0f64;
    let mut fnorm1p = 0f64;
    let mut gnorm = 0f64;
    let mut dirder = 0f64;

    let mut iter = 0i32;

    let mut p1 = 0.1f64;
    let mut p25 = 0.25f64;
    let mut p5 = 0.5f64;
    let mut p75 = 0.75f64;
    let mut p0001 = 0.0001f64;

    if state.fnorm == 0f64 {
        return enums::value::Success;
    }

    /* Compute qtf = Q^T f */

    qtf.copy_from(f);
    ::linear_algebra::QR_QTvec(r, tau, qtf);

    /* Compute norm of scaled gradient */

    compute_gradient_direction(r, perm, qtf, diag, gradient);

    { 
        size_t iamax = ::blas::level1::idamax(gradient);

        gnorm = unsafe { fabsf64(gradient.get(iamax) / state.fnorm) };
    }

    /* Determine the Levenberg-Marquardt parameter */

    loop {
        iter += 1;

        {
            let status = lmpar(r, perm, qtf, diag, state.delta, &mut (state.par), newton, gradient, sdiag, dx, w);

            if status != enums::value::Success {
                return status;
            }
        }

        /* Take a trial step */

        dx.scale(-1f64); /* reverse the step to go downhill */

        compute_trial_step(x, dx, state.x_trial);

        pnorm = scaled_enorm(diag, dx);

        if state.iter == 1 {
            if pnorm < state.delta {
        /*#ifdef DEBUG
              printf("set delta = pnorm = %g\n" , pnorm);
        #endif*/
                state.delta = pnorm;
            }
        }

        /* Evaluate function at x + p */
        /* return immediately if evaluation raised error */
        {
            let status = fdf.f(x_trial, f.params, f_trial);
            
            if status != enums::value::Success {
                return status;
            }
        }

        fnorm1 = enorm(f_trial);

        /* Compute the scaled actual reduction */

        actred = compute_actual_reduction(state.fnorm, fnorm1);

        /*#ifdef DEBUG
            printf("lmiterate: fnorm = %g fnorm1 = %g  actred = %g\n", state->fnorm, fnorm1, actred);
            printf("r = "); gsl_matrix_fprintf(stdout, r, "%g");
            printf("perm = "); gsl_permutation_fprintf(stdout, perm, "%d");
            printf("dx = "); gsl_vector_fprintf(stdout, dx, "%g");
        #endif*/

        /* Compute rptdx = R P^T dx, noting that |J dx| = |R P^T dx| */

        compute_rptdx(r, perm, dx, rptdx);

        /*#ifdef DEBUG
        printf("rptdx = "); gsl_vector_fprintf(stdout, rptdx, "%g");
        #endif*/

        fnorm1p = enorm(rptdx);

        /* Compute the scaled predicted reduction = |J dx|^2 + 2 par |D dx|^2 */

        { 
            let t1 = fnorm1p / state.fnorm;
            let t2 = (state.par.sqrt() * pnorm) / state.fnorm;

            prered = t1 * t1 + t2 * t2 / p5;
            dirder = -(t1 * t1 + t2 * t2);
        }

        /* compute the ratio of the actual to predicted reduction */

        let ratio = if prered > 0f64 {
            actred / prered
        } else {
            0f64
        };

        /*#ifdef DEBUG
        printf("lmiterate: prered = %g dirder = %g ratio = %g\n", prered, dirder,ratio);
        #endif*/

        /* update the step bound */
        if ratio > p25 {
    /*#ifdef DEBUG
          printf("ratio > p25\n");
    #endif*/
            if state.par == 0f64 || ratio >= p75 {
                state->delta = pnorm / p5;
                state->par *= p5;
    /*#ifdef DEBUG
              printf("updated step bounds: delta = %g, par = %g\n", state->delta, state->par);
    #endif*/
            }
        } else {
            let temp = if actred >= 0f64 {p5} else {p5 * dirder / (dirder + p5 * actred)};

    /*#ifdef DEBUG
          printf("ratio < p25\n");
    #endif*/

            if p1 * fnorm1 >= state->fnorm || temp < p1 {
                temp = p1;
            }

            state.delta = temp * state.delta.min(pnorm / p1);

            state.par /= temp;
    /*#ifdef DEBUG
          printf("updated step bounds: delta = %g, par = %g\n", state->delta, state->par);
    #endif*/
        }


        /* test for successful iteration, termination and stringent tolerances */

        if ratio >= p0001 {
            x.copy_from(x_trial);
            f.copy_from(f_trial);

            /* return immediately if evaluation raised error */
            {
                let status = if fdf.df {
                    fdf.df(x_trial, fdf.params, j)
                } else {
                    gsl_multifit_fdfsolver_dif_df(x_trial, fdf, f_trial, j)
                };

                if status != enums::value::Success {
                    return status;
                }
            }

            /* wa2_j  = diag_j * x_j */
            state.xnorm = scaled_enorm(diag, x);
            state.fnorm = fnorm1;
            state.iter += 1;

            /* Rescale if necessary */
            if scale != 0 {
                update_diag(j, diag);
            }

            {
                let mut signum = 0;

                r.copy_from(j);
                ::linear_algebra::QRPT_decomp(r, tau, perm, &mut signum, work1);
            }

            return enums::value::Success;
        } else if actred.fabs() <= ::DBL_EPSILON  && prered <= ::DBL_EPSILON  && p5 * ratio <= 1f64 {
            return enums::value::TolF;
        } else if state.delta <= ::DBL_EPSILON * state.xnorm {
            return enums::value::TolX;
        } else if gnorm <= ::DBL_EPSILON {
            return enums::value::TolG;
        } else if iter >= 10 {
            /* Repeat inner loop if unsuccessful */
            break;
        }
    }

    return enums::value::NoProg;
}