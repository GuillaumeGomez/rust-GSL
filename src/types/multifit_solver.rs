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
use libc::c_void;

/*pub struct MultiFitFSolver {
    s: *mut ffi::gsl_multifit_fsolver,
}

impl MultiFitFSolver {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n
    /// observations and p parameters. The number of observations n must be greater than or equal to
    /// parameters p.
    ///
    /// If there is insufficient memory to create the solver then the function returns a null
    /// pointer and the error handler is invoked with an error code of `Value::NoMemory`.
    pub fn new(t: *mut gsl_multifit_fsolver_type, n: usize, p: usize) -> Option<MultiFitFSolver> {
        let tmp = unsafe { ffi::gsl_multifit_fsolver_alloc(ffi::unwrap(t), n, p) };

        if tmp.is_null() {
            None
        } else {
            Some(MultiFitFSolver {
                s: tmp
            })
        }
    }

    pub fn set(&self, gsl_multifit_function * f, const gsl_vector * x) -> ::Value {
        unsafe { gsl_multifit_fsolver_set(self.s, gsl_multifit_function * f, const gsl_vector * x)
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

    fn unwrap(s: &MultiFitFSolver) -> *mut ffi::gsl_multifit_fsolver {
        s.s
    }
}*/

pub struct MultiFitFunction<'r, T:'r> {
    pub f: fn(x: &::VectorF64, params: &mut T, f: &::VectorF64) -> ::Value,
    /// number of functions
    pub n: usize,
    /// number of independent variables
    pub p: usize,
    pub params: &'r mut T
}

pub struct MultiFitFdfSolver<'r, T:'r> {
    _type: &'r MultiFitFdfSolverType<T>,
    fdf: *mut c_void, //MultiFitFunctionFdf<'r, T>,
    pub x: ::VectorF64,
    pub f: ::VectorF64,
    pub j: ::MatrixF64,
    pub dx: ::VectorF64,
    state: LmderStateT
}

impl<'r, T> MultiFitFdfSolver<'r, T> {
    /// This function returns a pointer to a newly allocated instance of a solver of type T for n observations and p parameters. The number
    /// of observations n must be greater than or equal to parameters p.
    pub fn new(_type: &'r MultiFitFdfSolverType<T>, n: usize, p: usize) -> Option<MultiFitFdfSolver<'r, T>> {
        let mut r = MultiFitFdfSolver {
            _type: _type,
            fdf: ::std::ptr::null_mut(),
            x: ::VectorF64::new(p).unwrap(),
            f: ::VectorF64::new(n).unwrap(),
            j: ::MatrixF64::new(n, p).unwrap(),
            dx: ::VectorF64::new(p).unwrap(),
            state: LmderStateT::new(n, p)
        };
        if ((*_type).alloc)(&mut r.state, n, p) == ::Value::Success {
            Some(r)
        } else {
            None
        }
    }

    /// This function initializes, or reinitializes, an existing solver s to use the function f and the initial guess x.
    pub fn set(&mut self, f: &mut MultiFitFunctionFdf<'r, T>, x: &::VectorF64) -> ::Value {
        if self.f.len() != f.n {
            rgsl_error!("function size does not match solver", ::Value::BadLength);
            return ::Value::BadLength;
        }

        if self.x.len() != x.len() {
            rgsl_error!("vector length does not match solver", ::Value::BadLength);
            return ::Value::BadLength;
        }

        self.fdf = unsafe { ::std::mem::transmute(f) };
        self.x.copy_from(x);

        unsafe {
            (self._type.set)(&mut self.state, ::std::mem::transmute(self.fdf), &mut self.x, &mut self.f, &mut self.j, &mut self.dx)
        }
    }

    pub fn name(&self) -> &'static str {
        self._type.name
    }

    /// This function performs a single iteration of the solver s. If the iteration encounters an unexpected problem then an error code
    /// will be returned. The solver maintains a current estimate of the best-fit parameters at all times.
    pub fn iterate(&mut self) -> ::Value {
        unsafe {
            (self._type.iterate)(&mut self.state, ::std::mem::transmute(self.fdf), &mut self.x, &mut self.f, &mut self.j, &mut self.dx)
        }
    }

    /// This function returns the current position (i.e. best-fit parameters) s->x of the solver s.
    pub fn position(&'r self) -> &'r ::VectorF64 {
        &self.x
    }

    /// These functions iterate the solver s for a maximum of maxiter iterations. After each iteration, the system is tested for convergence
    /// using gsl_multifit_test_delta with the error tolerances epsabs and epsrel.
    #[allow(unused_assignments)]
    pub fn driver(&mut self, max_iter: usize, epsabs: f64, epsrel: f64) -> ::Value {
        let mut status = ::Value::Success;
        let mut iter = 0usize;

        loop {
            status = self.iterate();

            if status != ::Value::Success {
                break;
            }

            /* test for convergence */
            status = gsl_multifit_test_delta(&self.dx, &self.x, epsabs, epsrel);
            iter += 1;
            if status != ::Value::Continue || iter >= max_iter {
                break;
            }
        }

        status
    }
}

impl<'r, T> Drop for MultiFitFdfSolver<'r, T> {
    fn drop(&mut self) {
        //self.s = ::std::ptr::null_mut();
    }
}

#[allow(dead_code)]
pub struct MultiFitFdfSolverType<T> {
    name: &'static str,
    //size: usize,
    alloc: fn(state: &mut LmderStateT, n: usize, p: usize) -> ::Value,
    set: fn(state: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, j: &mut ::MatrixF64,
        dx: &mut ::VectorF64) -> ::Value,
    iterate: fn(state: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64,
        j: &mut ::MatrixF64, dx: &mut ::VectorF64) -> ::Value,
    free: fn(state: &mut LmderStateT)
}

impl<T> MultiFitFdfSolverType<T> {
    pub fn lmder() -> MultiFitFdfSolverType<T> {
        MultiFitFdfSolverType {
            name: "lmder",
            alloc: lmder_alloc,
            set: lmder_set,
            iterate: lmder_iterate,
            free: lmder_free
        }
    }

    pub fn lmsder() -> MultiFitFdfSolverType<T> {
        MultiFitFdfSolverType {
            name: "lmsder",
            alloc: lmder_alloc,
            set: lmsder_set,
            iterate: lmder_iterate,
            free: lmder_free
        }
    }
}

pub struct MultiFitFunctionFdf<'r, T:'r> {
    pub f: fn(x: &::VectorF64, params: &mut T, f: &mut ::VectorF64) -> ::Value,
    pub df: Option<fn(x: &::VectorF64, params: &mut T, df: &mut ::MatrixF64) -> ::Value>,
    pub fdf: Option<fn(x: &::VectorF64, params: &mut T, f: &mut ::VectorF64, df: &mut ::MatrixF64) -> ::Value>,
    pub n: usize,
    pub p: usize,
    pub params: &'r mut T
}

#[allow(dead_code)]
struct LmderStateT {
    iter: usize,
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

impl LmderStateT {
    fn new(n: usize, p: usize) -> LmderStateT {
        LmderStateT {
            iter: 0usize,
            xnorm: 0f64,
            fnorm: 0f64,
            delta: 0f64,
            par: 0f64,
            r: ::MatrixF64::new(n, p).unwrap(),
            tau: ::VectorF64::new(if p < n {p} else {n}).unwrap(),
            diag: ::VectorF64::new(p).unwrap(),
            qtf: ::VectorF64::new(n).unwrap(),
            newton: ::VectorF64::new(p).unwrap(),
            gradient: ::VectorF64::new(p).unwrap(),
            x_trial: ::VectorF64::new(p).unwrap(),
            f_trial: ::VectorF64::new(n).unwrap(),
            df: ::VectorF64::new(n).unwrap(),
            sdiag: ::VectorF64::new(p).unwrap(),
            rptdx: ::VectorF64::new(n).unwrap(),
            w: ::VectorF64::new(n).unwrap(),
            work1: ::VectorF64::new(p).unwrap(),
            perm: ::Permutation::new(p).unwrap()
        }
    }
}

#[allow(unused_variables)]
// I'm not sure if I'll keep this function or not...
fn lmder_alloc(state: &mut LmderStateT, n: usize, p: usize) -> ::Value {
    /*state.r = match ::MatrixF64::new(n, p) {
        Some(m) => m,
        None => {
            rgsl_error!("failed to allocate space for r", ::Value::NoMem);
            ::MatrixF64::new(0, 0).unwrap() // to allow compilation
        }
    };

    state.tau = match ::VectorF64::new(if n < p {n} else {p}) {
        Some(t) => t,
        None => {
            rgsl_error!("failed to allocate space for tau", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.diag = match ::VectorF64::new(p) {
        Some(d) => d,
        None => {
            rgsl_error!("failed to allocate space for diag", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.qtf = match ::VectorF64::new(n) {
        Some(q) => q,
        None => {
            rgsl_error!("failed to allocate space for qtf", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.newton = match ::VectorF64::new(p) {
        Some(n) => n,
        None => {
            rgsl_error!("failed to allocate space for newton", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.gradient = match ::VectorF64::new(p) {
        Some(g) => g,
        None => {
            rgsl_error!("failed to allocate space for gradient", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.x_trial = match ::VectorF64::new(p) {
        Some(x) => x,
        None => {
            rgsl_error!("failed to allocate space for x_trial", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.f_trial = match ::VectorF64::new(n) {
        Some(f) => f,
        None => {
            rgsl_error!("failed to allocate space for f_trial", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.df = match ::VectorF64::new(n) {
        Some(d) => d,
        None => {
            rgsl_error!("failed to allocate space for df", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.sdiag = match ::VectorF64::new(p) {
        Some(s) => s,
        None => {
            rgsl_error!("failed to allocate space for sdiag", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.rptdx = match ::VectorF64::new(n) {
        Some(s) => s,
        None => {
            rgsl_error!("failed to allocate space for rptdx", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.w = match ::VectorF64::new(n) {
        Some(w) => w,
        None => {
            rgsl_error!("failed to allocate space for w", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.work1 = match ::VectorF64::new(p) {
        Some(w) => w,
        None => {
            rgsl_error!("failed to allocate space for work1", ::Value::NoMem);
            ::VectorF64::new(0).unwrap() // to allow compilation
        }
    };

    state.perm = match ::Permutation::new(p) {
        Some(p) => p,
        None => {
            rgsl_error!("failed to allocate space for perm", ::Value::NoMem);
            ::Permutation::new(0).unwrap() // to allow compilation
        }
    };*/

    ::Value::Success
}

#[allow(unused_variables)]
fn lmder_free(state: &mut LmderStateT) {
    /*::std::mem::drop(state.perm);
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
    ::std::mem::drop(state.r);*/
}

fn lmder_set<T>(vstate: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    set(vstate, fdf, x, f, J, dx, 0)
}

fn lmsder_set<T>(vstate: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    set(vstate, fdf, x, f, J, dx, 1)
}

fn lmder_iterate<T>(vstate: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    iterate(vstate, fdf, x, f, J, dx, 0)
}

#[allow(dead_code)]
fn lmsder_iterate<T>(vstate: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64) -> ::Value {
    iterate(vstate, fdf, x, f, J, dx, 1)
}

fn compute_diag(J: &::MatrixF64, diag: &mut ::VectorF64) {
    let n = J.size1();
    let p = J.size2();

    for j in 0..p {
        let mut sum = 0f64;

        for i in 0..n {
            let jij = J.get(i, j);

            sum += jij * jij;
        }
        if sum == 0f64 {
            sum = 0f64;
        }

        diag.set(j, sum.sqrt());
    }
}

fn update_diag(J: &::MatrixF64, diag: &mut ::VectorF64) {
    let n = diag.len();

    for j in 0..n {
        let mut sum = 0f64;

        for i in 0..n {
          let jij = J.get(i, j);

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
    let n = f.len();

    for i in 0..n {
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

    if Dx > 0f64 {
        factor * Dx
    } else {
        factor
    }
}

fn set<T>(state: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, j: &mut ::MatrixF64,
    dx: &mut ::VectorF64, scale: i32) -> ::Value {
    let mut signum = 0i32;

    /* Evaluate function at x */
    /* return immediately if evaluation raised error */
    {
        let status = match fdf.fdf {
            Some(func) => {
                func(x, fdf.params, f, j)
            }
            None => {
                /* finite difference approximation */
                gsl_multifit_fdfsolver_dif_fdf(x, fdf, f, j)
            }
        };

        if status != ::Value::Success {
            return status;
        }
    }

    state.par = 0f64;
    state.iter = 1;
    state.fnorm = enorm(f);

    dx.set_all(0f64);

    /* store column norms in diag */
    if scale != 0 {
        compute_diag(j, &mut state.diag);
    } else {
        state.diag.set_all(1f64);
    }

    /* set delta to 100 |D x| or to 100 if |D x| is zero */
    state.xnorm = scaled_enorm(&state.diag, x);
    state.delta = compute_delta(&mut state.diag, x);

    /* Factorize J into QR decomposition */
    state.r.copy_from(j);
    ::linear_algebra::QRPT_decomp(&state.r, &state.tau, &state.perm, &mut signum, &state.work1);

    state.rptdx.set_zero();
    state.w.set_zero();

    /* Zero the trial vector, as in the alloc function */

    state.f_trial.set_zero();

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

    for j in 0..n {
        let mut sum = 0f64;

        for i in 0..(j + 1) {
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
    let n = x.len();

    for i in 0..n {
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
    let n = dx.len();

    for i in 0..n {
        let mut sum = 0f64;

        for j in i..n {
            let pj = p.get(j);

            sum += r.get(i, j) * dx.get(pj);
        }

        rptdx.set(i, sum);
    }
}

#[allow(unused_assignments)]
fn iterate<T>(state: &mut LmderStateT, fdf: &mut MultiFitFunctionFdf<T>, x: &mut ::VectorF64, f: &mut ::VectorF64, J: &mut ::MatrixF64,
    dx: &mut ::VectorF64, scale: i32) -> ::Value {
    let mut prered = 0f64;
    let mut actred = 0f64;
    let mut pnorm = 0f64;
    let mut fnorm1 = 0f64;
    let mut fnorm1p = 0f64;
    let mut gnorm = 0f64;
    let mut dirder = 0f64;

    let mut iter = 0i32;

    let p1 = 0.1f64;
    let p25 = 0.25f64;
    let p5 = 0.5f64;
    let p75 = 0.75f64;
    let p0001 = 0.0001f64;

    if state.fnorm == 0f64 {
        return ::Value::Success;
    }

    /* Compute qtf = Q^T f */

    state.qtf.copy_from(f);
    ::linear_algebra::QR_QTvec(&state.r, &state.tau, &state.qtf);

    /* Compute norm of scaled gradient */

    compute_gradient_direction(&state.r, &state.perm, &state.qtf, &state.diag, &mut state.gradient);

    {
        let iamax = ::blas::level1::idamax(&state.gradient);

        gnorm = unsafe { (state.gradient.get(iamax as usize) / state.fnorm).abs() };
    }

    /* Determine the Levenberg-Marquardt parameter */

    loop {
        iter += 1;

        {
            let status = lmpar(&mut state.r, &state.perm, &state.qtf, &state.diag, state.delta, &mut (state.par), &mut state.newton,
                &mut state.gradient, &mut state.sdiag, dx, &mut state.w);

            if status != ::Value::Success {
                return status;
            }
        }

        /* Take a trial step */

        dx.scale(-1f64); /* reverse the step to go downhill */

        compute_trial_step(x, dx, &mut state.x_trial);

        pnorm = scaled_enorm(&state.diag, dx);

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
            let status = ((*fdf).f)(&state.x_trial, fdf.params, &mut state.f_trial); //GSL_MULTIFIT_FN_EVAL_F (fdf, x_trial, f_trial);

            if status != ::Value::Success {
                return status;
            }
        }

        fnorm1 = enorm(&state.f_trial);

        /* Compute the scaled actual reduction */

        actred = compute_actual_reduction(state.fnorm, fnorm1);

        /*#ifdef DEBUG
            printf("lmiterate: fnorm = %g fnorm1 = %g  actred = %g\n", state->fnorm, fnorm1, actred);
            printf("r = "); gsl_matrix_fprintf(stdout, r, "%g");
            printf("perm = "); gsl_permutation_fprintf(stdout, perm, "%d");
            printf("dx = "); gsl_vector_fprintf(stdout, dx, "%g");
        #endif*/

        /* Compute rptdx = R P^T dx, noting that |J dx| = |R P^T dx| */

        compute_rptdx(&state.r, &state.perm, dx, &mut state.rptdx);

        /*#ifdef DEBUG
        printf("rptdx = "); gsl_vector_fprintf(stdout, rptdx, "%g");
        #endif*/

        fnorm1p = enorm(&state.rptdx);

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
            let mut temp = if actred >= 0f64 {p5} else {p5 * dirder / (dirder + p5 * actred)};

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
            x.copy_from(&state.x_trial);
            f.copy_from(&state.f_trial);

            /* return immediately if evaluation raised error */
            {
                let status = match fdf.df {
                    Some(func) => {
                        func(&state.x_trial, fdf.params, J)
                    }
                    None => {
                        gsl_multifit_fdfsolver_dif_df(&mut state.x_trial, fdf, &mut state.f_trial, J)
                    }
                };

                if status != ::Value::Success {
                    return status;
                }
            }

            /* wa2_j  = diag_j * x_j */
            state.xnorm = scaled_enorm(&state.diag, x);
            state.fnorm = fnorm1;
            state.iter += 1;

            /* Rescale if necessary */
            if scale != 0 {
                update_diag(J, &mut state.diag);
            }

            {
                let mut signum = 0;

                state.r.copy_from(J);
                ::linear_algebra::QRPT_decomp(&state.r, &state.tau, &state.perm, &mut signum, &state.work1);
            }

            return ::Value::Success;
        } else if actred.abs() <= ::DBL_EPSILON  && prered <= ::DBL_EPSILON  && p5 * ratio <= 1f64 {
            return ::Value::ToleranceF;
        } else if state.delta <= ::DBL_EPSILON * state.xnorm {
            return ::Value::ToleranceX;
        } else if gnorm <= ::DBL_EPSILON {
            return ::Value::ToleranceG;
        } else if iter >= 10 {
            /* stop inner loop if successful */
            break;
        }
    }

    return ::Value::NoProgress;
}

fn gsl_multifit_test_delta(dx: &::VectorF64, x: &::VectorF64, epsabs: f64, epsrel: f64) -> ::Value {
    let mut ok = false;

    if epsrel < 0f64 {
        rgsl_error!("relative tolerance is negative", ::Value::BadTolerance);
        return ::Value::BadTolerance;
    }

    for i in 0..x.len() {
        let xi = x.get(i);
        let dxi = dx.get(i);
        let tolerance = epsabs + epsrel * xi.abs();

        if dxi.abs() < tolerance {
            ok = true;
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

fn gsl_multifit_fdfsolver_dif_fdf<T>(x: &mut ::VectorF64, fdf: &mut MultiFitFunctionFdf<T>, f: &mut ::VectorF64,
    j: &mut ::MatrixF64) -> ::Value {
    let mut status = ((*fdf).f)(x, fdf.params, f); // GSL_MULTIFIT_FN_EVAL_F(fdf, x, f);
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

fn fdjac<T>(x: &mut ::VectorF64, fdf: &mut MultiFitFunctionFdf<T>, f: &::VectorF64, jm: &mut ::MatrixF64) -> ::Value {
    let mut status = ::Value::Success;
    let epsfcn = 0f64;
    let eps = (epsfcn.max(::DBL_EPSILON)).sqrt();

    for j in 0..fdf.p {
        let xj = x.get(j);

        /* use column j of J as temporary storage for f(x + dx) */
        let (mut v, _) = jm.get_col(j).unwrap();

        let mut h = eps * xj.abs();
        if h == 0f64 {
            h = eps;
        }

        /* perturb x_j to compute forward difference */
        x.set(j, xj + h);

        status = ((*fdf).f)(x, fdf.params, &mut v); //GSL_MULTIFIT_FN_EVAL_F (fdf, x, &v.vector);
        if status != ::Value::Success {
            return status;
        }

        /* restore x_j */
        x.set(j, xj);

        h = 1f64 / h;
        for i in 0..fdf.n {
            let fnext = v.get(i);
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

fn gsl_multifit_fdfsolver_dif_df<T>(x: &mut ::VectorF64, fdf: &mut MultiFitFunctionFdf<T>, f: &::VectorF64,
    J: &mut ::MatrixF64) -> ::Value {
    fdjac(x, fdf, f, J)
}

#[allow(unused_assignments)]
fn lmpar(r: &mut ::MatrixF64, perm: &::Permutation, qtf: &::VectorF64, diag: &::VectorF64, delta: f64, par_inout: &mut f64,
    newton: &mut ::VectorF64, gradient: &mut ::VectorF64, sdiag: &mut ::VectorF64, x: &mut ::VectorF64,
    w: &mut ::VectorF64) -> ::Value {
    let mut dxnorm = 0f64;
    let mut gnorm = 0f64;
    let mut fp = 0f64;
    let mut fp_old = 0f64;
    let mut par_lower = 0f64;
    let mut par_upper = 0f64;
    let mut par_c = 0f64;

    let mut par = *par_inout;
    let mut iter = 0usize;

    /*#ifdef DEBUG
    printf("ENTERING lmpar\n");
    #endif*/


    compute_newton_direction(r, perm, qtf, newton);

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
        x.copy_from(newton);
        /*#ifdef DEBUG
          printf ("took newton (fp = %g, delta = %g)\n", fp, delta);
        #endif*/

        *par_inout = 0f64;

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

    compute_newton_bound(r, newton, dxnorm, perm, diag, w);

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
        let wnorm = enorm(w);
        let phider = wnorm * wnorm;

        /* w == zero if r rank-deficient,
           then set lower bound to zero form MINPACK, lmder.f
           Hans E. Plesser 2002-02-25 (hans.plesser@itf.nlh.no) */
        par_lower = if wnorm > 0f64 {
            fp / (delta * phider)
        } else {
            0f64
        };
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

    if par_upper == 0f64 {
        par_upper = ::DBL_MIN / delta.min(0.1f64);
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

    if par == 0f64 {
        par = gnorm / dxnorm;
        /*#ifdef DEBUG
        printf("set par to gnorm/dxnorm = %g\n", par);
        #endif*/
    }

    /* Beginning of iteration */
    loop {
        iter += 1;

        /*#ifdef DEBUG
        printf("lmpar iteration = %d\n", iter);
        #endif*/

        //#ifdef BRIANSFIX
        /* Seems like this is described in the paper but not in the MINPACK code */

        // I keep the BRIANSFIX for the moment
        if par < par_lower || par > par_upper {
            par = (par_lower * par_upper).sqrt().max(0.001f64 * par_upper);
        }
        //#endif

        /* Evaluate the function at the current value of par */

        if par == 0f64 {
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
            let sqrt_par = par.sqrt();

            qrsolv(r, perm, sqrt_par, diag, qtf, x, sdiag, w);
        }

        dxnorm = scaled_enorm(diag, x);

        fp_old = fp;

        fp = dxnorm - delta;

        /*#ifdef DEBUG
        printf ("After qrsolv dxnorm = %g, delta = %g, fp = %g\n", dxnorm, delta, fp);
        printf ("sdiag = ") ; gsl_vector_fprintf(stdout, sdiag, "%g"); printf("\n");
        printf ("x = ") ; gsl_vector_fprintf(stdout, x, "%g"); printf("\n");
        printf ("r = ") ; gsl_matrix_fprintf(stdout, r, "%g"); printf("\nXXX\n");
        #endif*/

        /* If the function is small enough, accept the current value of par */

        if fp.abs() <= 0.1f64 * delta {
            break;
        }

        if par_lower == 0f64 && fp <= fp_old && fp_old < 0f64 {
            break;
        }

        /* Check for maximum number of iterations */

        if iter == 10 {
            break;
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

        if fp > 0f64 {
            if par > par_lower {
                par_lower = par;
                /*#ifdef DEBUG
                  printf("fp > 0: set par_lower = par = %g\n", par);
                #endif*/
            }
        } else if fp < 0f64 {
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

    }

    /*#ifdef DEBUG
    printf("LEAVING lmpar, par = %g\n", par);
    #endif*/

    *par_inout = par;

    ::Value::Success
}

fn compute_newton_direction(r: &::MatrixF64, perm: &::Permutation, qtf: &::VectorF64, x: &mut ::VectorF64) {
    /* Compute and store in x the Gauss-Newton direction. If the
     Jacobian is rank-deficient then obtain a least squares
     solution. */

    let n = r.size2();

    for i in 0..n {
        let qtfi = qtf.get(i);

        x.set(i, qtfi);
    }

    let nsing = count_nsing(r);

    /*#ifdef DEBUG
    printf("nsing = %d\n", nsing);
    printf("r = "); gsl_matrix_fprintf(stdout, r, "%g"); printf("\n");
    printf("qtf = "); gsl_vector_fprintf(stdout, x, "%g"); printf("\n");
    #endif*/

    if nsing < n {
        for i in nsing..n {
            x.set(i, 0f64);
        }
    }

    if nsing > 0usize {
        for j in nsing..0 {
            let rjj = r.get(j, j);
            let temp = x.get(j) / rjj;

            x.set(j, temp);

            for i in 0..j {
                let rij = r.get(i, j);
                let xi = x.get(i);

                x.set(i, xi - rij * temp);
            }
        }
    }

    perm.permute_vector_inverse(x);
}

fn compute_newton_bound(r: &::MatrixF64, x: &::VectorF64, dxnorm: f64, perm: &::Permutation, diag: &::VectorF64,
    w: &mut ::VectorF64) {
    /* If the jacobian is not rank-deficient then the Newton step
     provides a lower bound for the zero of the function. Otherwise
     set this bound to zero. */

    let n = r.size2();

    let nsing = count_nsing(r);

    if nsing < n {
        w.set_zero();
        return;
    }

    for i in 0..n {
        let pi = perm.get(i);

        let dpi = diag.get(pi);
        let xpi = x.get(pi);

        w.set(i, dpi * (dpi * xpi / dxnorm));
    }

    for j in 0..n {
        let mut sum = 0f64;

        for i in 0..j {
            sum += r.get(i, j) * w.get(i);
        }

        {
            let rjj = r.get(j, j);
            let wj = w.get(j);

            w.set(j, (wj - sum) / rjj);
        }
    }
}

/* This function computes the solution to the least squares system
   phi = [ A x =  b , lambda D x = 0 ]^2

   where A is an M by N matrix, D is an N by N diagonal matrix, lambda
   is a scalar parameter and b is a vector of length M.
   The function requires the factorization of A into A = Q R P^T,
   where Q is an orthogonal matrix, R is an upper triangular matrix
   with diagonal elements of non-increasing magnitude and P is a
   permuation matrix. The system above is then equivalent to
   [ R z = Q^T b, P^T (lambda D) P z = 0 ]
   where x = P z. If this system does not have full rank then a least
   squares solution is obtained.  On output the function also provides
   an upper triangular matrix S such that
   P^T (A^T A + lambda^2 D^T D) P = S^T S
   Parameters,

   r: On input, contains the full upper triangle of R. On output the
   strict lower triangle contains the transpose of the strict upper
   triangle of S, and the diagonal of S is stored in sdiag.  The full
   upper triangle of R is not modified.
   p: the encoded form of the permutation matrix P. column j of P is
   column p[j] of the identity matrix.
   lambda, diag: contains the scalar lambda and the diagonal elements
   of the matrix D
   qtb: contains the product Q^T b
   x: on output contains the least squares solution of the system
   wa: is a workspace of length N
   */
fn qrsolv(r: &mut ::MatrixF64, p: &::Permutation, lambda: f64, diag: &::VectorF64, qtb: &::VectorF64,
    x: &mut ::VectorF64, sdiag: &mut ::VectorF64, wa: &mut ::VectorF64) -> ::Value {
    let n = r.size2();

    /* Copy r and qtb to preserve input and initialise s. In particular,
     save the diagonal elements of r in x */

    let mut j = 0;
    while j < n {
        let rjj = r.get(j, j);
        let qtbj = qtb.get(j);

        let mut i = j + 1;
        while i < n {
            let rji = r.get(j, i);

            r.set(i, j, rji);
            i += 1;
        }

        x.set(j, rjj);
        wa.set(j, qtbj);
        j += 1;
    }

    /* Eliminate the diagonal matrix d using a Givens rotation */
    j = 0;
    while j < n {
        let pj = p.get(j);

        let diagpj = lambda * diag.get(pj);

        if diagpj == 0f64 {
            continue;
        }

        sdiag.set(j, diagpj);

        let mut k = j + 1;
        while k < n {
            sdiag.set(k, 0f64);
            k += 1;
        }

        /* The transformations to eliminate the row of d modify only a
           single element of qtb beyond the first n, which is initially
           zero */

        let mut qtbpj = 0f64;
        k = j;
        while k < n {
            /* Determine a Givens rotation which eliminates the
               appropriate element in the current row of d */

            let wak = wa.get(k);
            let rkk = r.get(k, k);
            let sdiagk = sdiag.get(k);

            if sdiagk == 0f64 {
                continue;
            }

            let (sine, cosine) = if rkk.abs() < sdiagk.abs() {
                let cotangent = rkk / sdiagk;
                let t_sine = 0.5f64 / (0.25f64 + 0.25f64 * cotangent * cotangent).sqrt();

                (t_sine, t_sine * cotangent)
            } else {
                let tangent = sdiagk / rkk;
                let t_cos = 0.5f64 / (0.25f64 + 0.25f64 * tangent * tangent).sqrt();

                (t_cos * tangent, t_cos)
            };

            /* Compute the modified diagonal element of r and the
               modified element of [qtb,0] */

            {
                let new_rkk = cosine * rkk + sine * sdiagk;
                let new_wak = cosine * wak + sine * qtbpj;

                qtbpj = -sine * wak + cosine * qtbpj;

                r.set(k, k, new_rkk);
                wa.set(k, new_wak);
            }

            /* Accumulate the transformation in the row of s */
            let mut i = k + 1;
            while i < n {
                let rik = r.get(i, k);
                let sdiagi = sdiag.get(i);

                let new_rik = cosine * rik + sine * sdiagi;
                let new_sdiagi = -sine * rik + cosine * sdiagi;

                r.set(i, k, new_rik);
                sdiag.set(i, new_sdiagi);
                i += 1;
            }
            k += 1;
        }

        /* Store the corresponding diagonal element of s and restore the
           corresponding diagonal element of r */

        {
            let rjj = r.get(j, j);
            let xj = x.get(j);

            sdiag.set(j, rjj);
            r.set(j, j, xj);
        }
        j += 1;

    }

    /* Solve the triangular system for z. If the system is singular then
       obtain a least squares solution */

    let mut nsing = n;
    j = 0;
    while j < n {
        let sdiagj = sdiag.get(j);

        if sdiagj == 0f64 {
            nsing = j;
            break;
        }
        j += 1;
    }

    j = nsing;
    while j < n {
        wa.set(j, 0f64);
        j += 1;
    }

    let mut k = 0;
    while k < nsing {
        let mut sum = 0f64;

        j = (nsing - 1) - k;

        let mut i = j + 1;
        while i < nsing {
            sum += r.get(i, j) * wa.get(i);
            i += 1;
        }

        {
            let waj = wa.get(j);
            let sdiagj = sdiag.get(j);

            wa.set(j, (waj - sum) / sdiagj);
        }
        k += 1;
    }

    /* Permute the components of z back to the components of x */
    j = 0;
    while j < n {
        let pj = p.get(j);
        let waj = wa.get(j);

        x.set(pj, waj);
        j += 1;
    }

    ::Value::Success
}

fn compute_newton_correction(r: &::MatrixF64, sdiag: &::VectorF64, p: &::Permutation, x: &mut ::VectorF64, dxnorm: f64,
    diag: &::VectorF64, w: &mut ::VectorF64) {
    let n = r.size2();

    let mut i = 0;
    while i < n {
      let pi = p.get(i);

      let dpi = diag.get(pi);
      let xpi = x.get(pi);

      w.set(i, dpi * (dpi * xpi) / dxnorm);
      i += 1;
    }

    let mut j = 0;
    while j < n {
        let sj = sdiag.get(j);
        let wj = w.get(j);

        let tj = wj / sj;

        w.set(j, tj);

        let mut i = j + 1;
        while i < n {
            let rij = r.get (i, j);
            let wi = w.get (i);

            w.set (i, wi - rij * tj);
            i += 1;
        }
        j += 1;
    }
}

fn count_nsing(r: &::MatrixF64) -> usize {
    /* Count the number of nonsingular entries. Returns the index of the
       first entry which is singular. */

    let n = r.size2();
    let mut i = 0usize;

    while i < n {
        let rii = r.get(i, i);

        if rii == 0f64 {
            break;
        }
        i += 1;
    }

    i
}
