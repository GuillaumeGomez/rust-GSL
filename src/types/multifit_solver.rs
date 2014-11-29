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
use libc::c_void;

pub struct MultiFitFunction<'r, T> {
    pub f: fn(x: &::VectorF64, params: &mut T, f: &::VectorF64) -> ::Value,
    /// number of functions
    pub n: u64,
    /// number of independent variables
    pub p: u64,
    pub params: &'r mut T
}

pub struct MultiFitFdfSolver<T> {
    _type: MultiFitFdfSolverType,
    fdf: *mut c_void,//MultiFitFunctionFdf<T>,
    x: ::VectorF64,
    f: ::VectorF64,
    j: ::MatrixF64,
    dx: ::VectorF64,
    state: *mut c_void
}

impl<T> MultiFitFdfSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n observations and p parameters. The number
    /// of observations n must be greater than or equal to parameters p.
    pub fn new(_type: &MultiFitFdfSolverType, n: u64, p: u64) -> Option<MultiFitFdfSolver<T>> {
        let r = MultiFitFdfSolver {
            _type: *_type,
            fdf: ::std::ptr::null_mut(),
            x: ::VectorF64::new(p).unwrap(),
            f: ::VectorF64::new(n).unwrap(),
            j: ::MatrixF64::new(n, p).unwrap(),
            dx: ::VectorF64::new(p).unwrap(),
            state: ::std::ptr::null_mut()
        };
        if _type.alloc(r.state, n, p) == ::Value::Success {
            Some(r)
        } else {
            None
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function f and the initial guess x.
    pub fn set<U>(&mut self, f: &MultiFitFunction<U>, x: &::VectorF64) -> ::Value {
        if self.f.len() != f.n {
            rgsl_error!("function size does not match solver", ::Value::BadLen);
            return ::Value::BadLen;
        }

        if self.x.len() != x.len() {
            rgsl_error!("vector length does not match solver", ::Value::BadLen);
            return ::Value::BadLen;
        }  

        self.fdf = *f;
        self.x.copy_from(x);

        self._type.set(self.state, self.fdf, self.x, self.f, self.j, self.dx)
    }

    pub fn name(&self) -> String {
        self._type.name
    }

    /// This function performs a single iteration of the solver s. If the iteration encounters an unexpected problem then an error code
    /// will be returned. The solver maintains a current estimate of the best-fit parameters at all times.
    pub fn iterate(&self) -> ::Value {
        self._type.iterate(self.state, self.fdf, self.x, self.f, self.j, self.dx)
    }

    /// This function returns the current position (i.e. best-fit parameters) s->x of the solver s.
    pub fn position<'r>(&'r self) -> &'r ::VectorF64 {
        self.x
    }

    /// These functions iterate the solver s for a maximum of maxiter iterations. After each iteration, the system is tested for convergence
    /// using gsl_multifit_test_delta with the error tolerances epsabs and epsrel.
    pub fn driver(&self, max_iter: u64, epsabs: f64, epsrel: f64) -> ::Value {
        let mut status = ::Value::Success;
        let mut iter = 0u64;

        loop {
            status = self.iterate();

            if status != ::Value::Success {
                break;
            }

            /* test for convergence */
            status = gsl_multifit_test_delta(self.dx, self.x, epsabs, epsrel);
            iter += 1;
            if status != ::Value::Continue || iter >= max_iter {
                break;
            }
        }

        status
    }
}

impl<T> Drop for MultiFitFdfSolver<T> {
    fn drop(&mut self) {
        //self.s = ::std::ptr::null_mut();
    }
}

pub struct MultiFitFdfSolverType {
    name: &'static str,
    //size: u64,
    alloc: fn(state: &mut LmderStateT, n: u64, p: u64) -> ::Value,
    set: fn(state: &mut LmderStateT, fdf: MultiFitFunctionFdf, x: &mut ::VectorF64, f: &mut ::VectorF64, j: ::MatrixF64,
        dx: &mut ::VectorF64) -> ::Value,
    iterate: fn(state: &mut LmderStateT, fdf: MultiFitFunctionFdf, x: &mut ::VectorF64, f: ::VectorF64, j: &mut ::MatrixF64,
        dx: &mut ::VectorF64) -> ::Value,
    free: fn(state: &mut LmderStateT)
}

impl MultiFitFdfSolverType {
    pub fn lmder_type() -> MultiFitFdfSolverType {
        MultiFitFdfSolverType {
            name: "lmder",
            alloc: lmder_alloc,
            set: lmder_set,
            iterate: lmder_iterate,
            free: lmder_free
        }
    }

    pub fn lmsder_type() -> MultiFitFdfSolverType {
        MultiFitFdfSolverType {
            name: "lmsder",
            alloc: lmder_alloc,
            set: lmsder_set,
            iterate: lmder_iterate,
            free: lmder_free
        }
    }
}

struct MultiFitFunctionFdf<T> {
    f: Option<fn(x: &::VectorF64, params: &mut T, f: &mut ::VectorF64) -> ::Value>,
    df: Option<fn(x: &::VectorF64, params: &mut T, df: &mut ::MatrixF64) -> ::Value>,
    fdf: Option<fn(x: &::VectorF64, params: &mut T, f: &mut ::VectorF64, df: &mut ::MatrixF64) -> ::Value>,
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

fn lmder_alloc(state: &mut LmderStateT, n: u64, p: u64) -> ::Value {
    state.r = match ::MatrixF64::new(n, p) {
        Some(m) => m,
        None => {
            rgsl_error!("failed to allocate space for r", ::Value::NoMem)
        }
    };

    state.tau = match ::VectorF64::new(n.min(p)) {
        Some(t) => t,
        None => {
            ::std::mem::drop(state.r);
            rgsl_error!("failed to allocate space for tau", ::Value::NoMem)
        }
    };

    state.diag = match ::VectorF64::new(p) {
        Some(d) => d,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            rgsl_error!("failed to allocate space for diag", ::Value::NoMem)
        }
    };

    state.qtf = match ::VectorF64::new(n) {
        Some(q) => q,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            rgsl_error!("failed to allocate space for qtf", ::Value::NoMem)
        }
    };

    state.newton = match ::VectorF64::new(p) {
        Some(n) => n,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            rgsl_error!("failed to allocate space for newton", ::Value::NoMem)
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
            rgsl_error!("failed to allocate space for gradient", ::Value::NoMem)
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
            rgsl_error!("failed to allocate space for x_trial", ::Value::NoMem)
        }
    };

    state.f_trial = match ::VectorF64::new(n) {
        Some(f) => f,
        None => {
            ::std::mem::drop(state.r);
            ::std::mem::drop(state.tau);
            ::std::mem::drop(state.diag);
            ::std::mem::drop(state.qtf);
            ::std::mem::drop(state.newton);
            ::std::mem::drop(state.gradient);
            ::std::mem::drop(state.x_trial);
            rgsl_error!("failed to allocate space for f_trial", ::Value::NoMem)
        }
    };

    state.df = match ::VectorF64::new(n) {
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
            rgsl_error!("failed to allocate space for df", ::Value::NoMem)
        }
    };

    state.sdiag = match ::VectorF64::new(p) {
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
            rgsl_error!("failed to allocate space for sdiag", ::Value::NoMem)
        }
    };

    state.rptdx = match ::VectorF64::new(n) {
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
            rgsl_error!("failed to allocate space for rptdx", ::Value::NoMem)
        }
    };

    state.w = match ::VectorF64::new(n) {
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
            rgsl_error!("failed to allocate space for w", ::Value::NoMem)
        }
    };

    state.work1 = match ::VectorF64::new(p) {
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
            rgsl_error!("failed to allocate space for work1", ::Value::NoMem)
        }
    };

    state.perm = match ::Permutation::new(p) {
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
            rgsl_error!("failed to allocate space for perm", ::Value::NoMem)
        }
    };

    ::Value::Success
}

fn lmder_free(state: &mut LmderStateT) {
    ::std::mem::drop(state.perm);
    ::std::mem::drop(state.work1);
    ::std::mem::drop(state.w);
    ::std::mem::drop(state.rptdx);
    ::std::mem::drop(state.sdiag);
    ::std::mem::drop(state.df);
    ::std::mem::drop(state.f_trial);
    ::std::mem::drop(state.x_trial);
    ::std::mem::drop(state.gradient);
    ::std::mem::drop(state.newton);
    ::std::mem::drop(state.qtf);
    ::std::mem::drop(state.diag);
    ::std::mem::drop(state.tau);
    ::std::mem::drop(state.r);
}

fn lmder_set<T>(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    set(vstate, fdf, x, f, J, dx, 0)
}

fn lmsder_set<T>(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    set(vstate, fdf, x, f, J, dx, 1)
}

fn lmder_iterate<T>(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    iterate(vstate, fdf, x, f, J, dx, 0)
}

fn lmsder_iterate<T>(vstate: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    iterate(vstate, fdf, x, f, J, dx, 1)
}

fn compute_diag(j: &::MatrixF64, diag: &mut ::VectorF64) {
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

fn set<T>(state: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, j: &mut ::MatrixF64,
    dx: &mut ::VectorF64, scale: i32) -> ::Value {
    let mut signum = 0i32;

    /* Evaluate function at x */
    /* return immediately if evaluation raised error */
    {
        let status = if fdf.fdf.is_some() {
            fdf.fdf(x, fdf.params, f, j)
        } else {
            /* finite difference approximation */
            gsl_multifit_fdfsolver_dif_fdf(x, fdf, f, j)
        };

        if status != ::Value::Success {
            return status;
        }
    }

    state.par = 0;
    state.iter = 1;
    state.fnorm = enorm(f);

    dx.set_all(0f64);

    /* store column norms in diag */
    if scale != 0 {
        compute_diag(j, state.diag);
    } else {
        state.diag.set_all(1f64);
    }

    /* set delta to 100 |D x| or to 100 if |D x| is zero */
    state.xnorm = scaled_enorm(state.diag, x);
    state.delta = compute_delta(state.diag, x);

    /* Factorize J into QR decomposition */
    state.r.copy_from(j);
    ::linear_algebra::QRPT_decomp(state.r, state.tau, state.perm, &mut signum, state.work1);

    state.rptdx.set_zero();
    state.w.set_zero();

    /* Zero the trial vector, as in the alloc function */

    state.f_trials.set_zero();

    /*#ifdef DEBUG
    printf("r = "); gsl_matrix_fprintf(stdout, r, "%g");
    printf("perm = "); gsl_permutation_fprintf(stdout, perm, "%d");
    printf("tau = "); gsl_vector_fprintf(stdout, tau, "%g");
    #endif*/

    ::Value::Success
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

fn iterate<T>(state: &LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64, scale: i32) -> ::Value {
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
        return ::Value::Success;
    }

    /* Compute qtf = Q^T f */

    state.qtf.copy_from(f);
    ::linear_algebra::QR_QTvec(state.r, state.tau, state.qtf);

    /* Compute norm of scaled gradient */

    compute_gradient_direction(state.r, state.perm, state.qtf, state.diag, state.gradient);

    { 
        let iamax = ::blas::level1::idamax(state.gradient);

        gnorm = unsafe { (state.gradient.get(iamax) / state.fnorm).fabs() };
    }

    /* Determine the Levenberg-Marquardt parameter */

    loop {
        iter += 1;

        {
            let status = lmpar(state.r, state.perm, state.qtf, state.diag, state.delta, &mut (state.par), state.newton,
                state.gradient, state.sdiag, dx, state.w);

            if status != ::Value::Success {
                return status;
            }
        }

        /* Take a trial step */

        dx.scale(-1f64); /* reverse the step to go downhill */

        compute_trial_step(x, dx, state.x_trial);

        pnorm = scaled_enorm(state.diag, dx);

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
            let status = fdf.f(state.x_trial, f.params, state.f_trial);
            
            if status != ::Value::Success {
                return status;
            }
        }

        fnorm1 = enorm(state.f_trial);

        /* Compute the scaled actual reduction */

        actred = compute_actual_reduction(state.fnorm, fnorm1);

        /*#ifdef DEBUG
            printf("lmiterate: fnorm = %g fnorm1 = %g  actred = %g\n", state->fnorm, fnorm1, actred);
            printf("r = "); gsl_matrix_fprintf(stdout, r, "%g");
            printf("perm = "); gsl_permutation_fprintf(stdout, perm, "%d");
            printf("dx = "); gsl_vector_fprintf(stdout, dx, "%g");
        #endif*/

        /* Compute rptdx = R P^T dx, noting that |J dx| = |R P^T dx| */

        compute_rptdx(state.r, state.perm, dx, state.rptdx);

        /*#ifdef DEBUG
        printf("rptdx = "); gsl_vector_fprintf(stdout, rptdx, "%g");
        #endif*/

        fnorm1p = enorm(state.rptdx);

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
                state.delta = pnorm / p5;
                state.par *= p5;
    /*#ifdef DEBUG
              printf("updated step bounds: delta = %g, par = %g\n", state->delta, state->par);
    #endif*/
            }
        } else {
            let temp = if actred >= 0f64 {p5} else {p5 * dirder / (dirder + p5 * actred)};

    /*#ifdef DEBUG
          printf("ratio < p25\n");
    #endif*/

            if p1 * fnorm1 >= state.fnorm || temp < p1 {
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
            x.copy_from(state.x_trial);
            f.copy_from(state.f_trial);

            /* return immediately if evaluation raised error */
            {
                let status = if fdf.df {
                    fdf.df(state.x_trial, fdf.params, J)
                } else {
                    gsl_multifit_fdfsolver_dif_df(state.x_trial, fdf, state.f_trial, J)
                };

                if status != ::Value::Success {
                    return status;
                }
            }

            /* wa2_j  = diag_j * x_j */
            state.xnorm = scaled_enorm(state.diag, x);
            state.fnorm = fnorm1;
            state.iter += 1;

            /* Rescale if necessary */
            if scale != 0 {
                update_diag(J, state.diag);
            }

            {
                let mut signum = 0;

                state.r.copy_from(J);
                ::linear_algebra::QRPT_decomp(state.r, state.tau, state.perm, &mut signum, state.work1);
            }

            return ::Value::Success;
        } else if actred.fabs() <= ::DBL_EPSILON  && prered <= ::DBL_EPSILON  && p5 * ratio <= 1f64 {
            return ::Value::TolF;
        } else if state.delta <= ::DBL_EPSILON * state.xnorm {
            return ::Value::TolX;
        } else if gnorm <= ::DBL_EPSILON {
            return ::Value::TolG;
        } else if iter >= 10 {
            /* stop inner loop if successful */
            break;
        }
    }

    return ::Value::NoProg;
}

fn gsl_multifit_test_delta(dx: &::VectorF64, x: &::VectorF64, epsabs: f64, epsrel: f64) -> ::Value {
    let mut ok = false;

    if epsrel < 0f64 {
        rgsl_error!("relative tolerance is negative", ::Value::BadTol);
        return ::Value::BadTol;
    }

    for i in range(0, x.size) {
        let xi = x.get(i);
        let dxi = dx.get(i);
        let tolerance = epsabs + epsrel * xi.fabs();

        if dxi.fabs() < tolerance {
            ok = 1;
        } else {
            ok = false;
            break;
        }
    }

    if ok {
        ::Value::Success
    } else {
        ::Value::Continue
    }
}

/*
gsl_multifit_fdfsolver_dif_fdf()
  Compute function values (analytic) and approximate Jacobian using finite
differences
Inputs: x      - parameter vector
        fdf    - fdf
        f      - (output) function values f_i(x)
        J      - (output) approximate Jacobian matrix
Return: success or error
*/

fn gsl_multifit_fdfsolver_dif_fdf<T>(x: &::VectorF64, fdf: &MultiFitFunctionFdf<T>, f: &mut ::VectorF64,
    j: &mut ::MatrixF64) -> ::Value {
    let mut status = ::Value::Success;

    status = fdf.f(x, fdf.params, f); // GSL_MULTIFIT_FN_EVAL_F(fdf, x, f);
    if status != ::Value::Success {
        return status;
    }

    status = fdjac(x, fdf, f, j);
    if status != ::Value::Success {
        return status;
    }

    status
}

/*
fdjac()
  Compute approximate Jacobian using forward differences
Inputs: x   - parameter vector
        fdf - fdf struct
        f   - (input) vector of function values f_i(x)
        J   - (output) Jacobian matrix
Return: success or error
*/

fn fdjac<T>(x: &mut ::VectorF64, fdf: &MultiFitFunctionFdf<T>, f: &::VectorF64, jm: &mut ::MatrixF64) -> ::Value {
    let mut status = ::Value::Success;
    let epsfcn = 0f64;
    let eps = (epsfcn.max(::DBL_EPSILON)).sqrtf();

    for j in range (0, fdf.p) {
        let xj = x.get(j);

        /* use column j of J as temporary storage for f(x + dx) */
        let v = jm.column(j);

        let mut h = eps * xj.fabs();
        if h == 0f64 {
            h = eps;
        }

        /* perturb x_j to compute forward difference */
        x.set(j, xj + h);

        status += fdf.f(x, fdf.params, &v.vector); //GSL_MULTIFIT_FN_EVAL_F (fdf, x, &v.vector);
        if status != ::Value::Success {
            return status;
        }

        /* restore x_j */
        x.set(j, xj);

        h = 1f64 / h;
        for i in range(0, fdf.n) {
            let fnext = v.vector.get(i);
            let fi = f.get(i);

            jm.set(i, j, (fnext - fi) * h);
        }
    }

    status
}

/*
gsl_multifit_fdfsolver_dif_df()
  Compute approximate Jacobian using finite differences
Inputs: x   - parameter vector
        fdf - fdf
        f   - (input) function values f_i(x)
        J   - (output) approximate Jacobian matrix
Return: success or error
*/

fn gsl_multifit_fdfsolver_dif_df<T>(x: &::VectorF64, fdf: &MultiFitFunctionFdf<T>, f: &::VectorF64, J: &::MatrixF64) -> ::Value {
    fdjac(x, fdf, f, J)
}

fn lmpar(r: &mut ::MatrixF64, perm: &::Permutation, qtf: &::VectorF64, diag: &::VectorF64, delta: f64, par_inout: &mut f64,
    newton: &mut ::VectorF64, gradient: &mut ::VectorF64, sdiag: &mut ::VectorF64, x: &mut ::VectorF64,
    w: &mut ::VectorF64) -> ::Value {
    double dxnorm, gnorm, fp, fp_old, par_lower, par_upper, par_c;

    let mut par = *par_inout;
    size_t iter = 0;

    /*#ifdef DEBUG
    printf("ENTERING lmpar\n");
    #endif*/


    compute_newton_direction (r, perm, qtf, newton);

    /*#ifdef DEBUG
    printf ("newton = ");
    gsl_vector_fprintf (stdout, newton, "%g");
    printf ("\n");

    printf ("diag = ");
    gsl_vector_fprintf (stdout, diag, "%g");
    printf ("\n");
    #endif*/

    /* Evaluate the function at the origin and test for acceptance of
     the Gauss-Newton direction. */

    dxnorm = scaled_enorm(diag, newton);

    fp = dxnorm - delta;

    /*#ifdef DEBUG
    printf ("dxnorm = %g, delta = %g, fp = %g\n", dxnorm, delta, fp);
    #endif*/

    if fp <= 0.1 * delta {
        gsl_vector_memcpy (x, newton);
        /*#ifdef DEBUG
          printf ("took newton (fp = %g, delta = %g)\n", fp, delta);
        #endif*/

        *par_inout = 0;

        return ::Value::Success;
    }

    /*#ifdef DEBUG
    printf ("r = ");
    gsl_matrix_fprintf (stdout, r, "%g");
    printf ("\n");

    printf ("newton = ");
    gsl_vector_fprintf (stdout, newton, "%g");
    printf ("\n");

    printf ("dxnorm = %g\n", dxnorm);
    #endif*/

    compute_newton_bound (r, newton, dxnorm, perm, diag, w);

    /*#ifdef DEBUG
    printf("perm = "); gsl_permutation_fprintf(stdout, perm, "%d");

    printf ("diag = ");
    gsl_vector_fprintf (stdout, diag, "%g");
    printf ("\n");

    printf ("w = ");
    gsl_vector_fprintf (stdout, w, "%g");
    printf ("\n");
    #endif*/


    {
        double wnorm = enorm (w);
        double phider = wnorm * wnorm;

        /* w == zero if r rank-deficient, 
           then set lower bound to zero form MINPACK, lmder.f 
           Hans E. Plesser 2002-02-25 (hans.plesser@itf.nlh.no) */
        if wnorm > 0 {
            par_lower = fp / (delta * phider);
        } else {
            par_lower = 0.0;
        }
    }

    /*#ifdef DEBUG
    printf("par       = %g\n", par      );
    printf("par_lower = %g\n", par_lower);
    #endif*/

    compute_gradient_direction(r, perm, qtf, diag, gradient);

    gnorm = enorm(gradient);

    /*#ifdef DEBUG
    printf("gradient = "); gsl_vector_fprintf(stdout, gradient, "%g"); printf("\n");
    printf("gnorm = %g\n", gnorm);
    #endif*/

    par_upper =  gnorm / delta;

    if par_upper == 0 {
        par_upper = ::DBL_MIN / GSL_MIN_DBL(delta, 0.1);
    }

    /*#ifdef DEBUG
    printf("par_upper = %g\n", par_upper);
    #endif*/

    if par > par_upper {
        /*#ifdef DEBUG
        printf("set par to par_upper\n");
        #endif*/

        par = par_upper;
    } else if par < par_lower {
        /*#ifdef DEBUG
        printf("set par to par_lower\n");
        #endif*/

        par = par_lower;
    }

    if par == 0 {
        par = gnorm / dxnorm;
        /*#ifdef DEBUG
        printf("set par to gnorm/dxnorm = %g\n", par);
        #endif*/
    }

    /* Beginning of iteration */

iteration:

    iter++;

    /*#ifdef DEBUG
    printf("lmpar iteration = %d\n", iter);
    #endif*/

    //#ifdef BRIANSFIX
    /* Seems like this is described in the paper but not in the MINPACK code */

    // I keep the BRIANSFIX for the moment
    if par < par_lower || par > par_upper {
        par = (par_lower * par_upper).sqrtf().max(0.001f64 * par_upper);
    }
    //#endif

    /* Evaluate the function at the current value of par */

    if par == 0 {
        par = (0.001f64 * par_upper).max(::DBL_MIN);
        /*#ifdef DEBUG
          printf("par = 0, set par to  = %g\n", par);
        #endif*/
    }

    /* Compute the least squares solution of [ R P x - Q^T f, sqrt(par) D x]
     for A = Q R P^T */

    /*#ifdef DEBUG
    printf ("calling qrsolv with par = %g\n", par);
    #endif*/

    {
        let sqrt_par = sqrt(par);

        qrsolv(r, perm, sqrt_par, diag, qtf, x, sdiag, w);
    }

    dxnorm = scaled_enorm (diag, x);

    fp_old = fp;

    fp = dxnorm - delta;

    /*#ifdef DEBUG
    printf ("After qrsolv dxnorm = %g, delta = %g, fp = %g\n", dxnorm, delta, fp);
    printf ("sdiag = ") ; gsl_vector_fprintf(stdout, sdiag, "%g"); printf("\n");
    printf ("x = ") ; gsl_vector_fprintf(stdout, x, "%g"); printf("\n");
    printf ("r = ") ; gsl_matrix_fprintf(stdout, r, "%g"); printf("\nXXX\n");
    #endif*/

    /* If the function is small enough, accept the current value of par */

    if fp.fabs() <= 0.1f64 * delta {
        goto line220;
    }

    if par_lower == 0 && fp <= fp_old && fp_old < 0 {
        goto line220;
    }

    /* Check for maximum number of iterations */

    if iter == 10 {
        goto line220;
    }

    /* Compute the Newton correction */

    compute_newton_correction(r, sdiag, perm, x, dxnorm, diag, w);

    /*#ifdef DEBUG
    printf ("newton_correction = ");
    gsl_vector_fprintf(stdout, w, "%g"); printf("\n");
    #endif*/

    {
        let wnorm = enorm(w);
        par_c = fp / (delta * wnorm * wnorm);
    }

    /*#ifdef DEBUG
    printf("fp = %g\n", fp);
    printf("par_lower = %g\n", par_lower);
    printf("par_upper = %g\n", par_upper);
    printf("par_c = %g\n", par_c);
    #endif*/


    /* Depending on the sign of the function, update par_lower or par_upper */

    if fp > 0 {
        if par > par_lower {
            par_lower = par;
            /*#ifdef DEBUG
              printf("fp > 0: set par_lower = par = %g\n", par);
            #endif*/
        }
    } else if fp < 0 {
        if par < par_upper {
            /*#ifdef DEBUG
              printf("fp < 0: set par_upper = par = %g\n", par);
            #endif*/
            par_upper = par;
        }
    }

    /* Compute an improved estimate for par */

    /*#ifdef DEBUG
      printf("improved estimate par = MAX(%g, %g) \n", par_lower, par+par_c);
    #endif*/

    par = par_lower.max(par + par_c);

    /*#ifdef DEBUG
      printf("improved estimate par = %g \n", par);
    #endif*/

    goto iteration;

    line220:

    /*#ifdef DEBUG
    printf("LEAVING lmpar, par = %g\n", par);
    #endif*/

    *par_inout = par;

    ::Value::Success
}