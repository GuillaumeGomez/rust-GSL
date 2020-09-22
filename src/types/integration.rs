//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use enums;
use ffi::FFI;

/// The QAG algorithm is a simple adaptive integration procedure. The integration region is divided
/// into subintervals, and on each iteration the subinterval with the largest estimated error is
/// bisected. This reduces the overall error rapidly, as the subintervals become concentrated
/// around local difficulties in the integrand. These subintervals are managed by a
/// gsl_integration_workspace struct, which handles the memory for the subinterval ranges, results
/// and error estimates.
pub struct IntegrationWorkspace {
    w: *mut sys::gsl_integration_workspace,
}

impl IntegrationWorkspace {
    /// This function allocates a workspace sufficient to hold n double precision intervals, their
    /// integration results and error estimates. One workspace may be used multiple times as all
    /// necessary reinitialization is performed automatically by the integration routines.
    pub fn new(n: u64) -> Option<IntegrationWorkspace> {
        let tmp = unsafe { sys::gsl_integration_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(IntegrationWorkspace { w: tmp })
        }
    }

    /// This function applies an integration rule adaptively until an estimate of the integral of f
    /// over (a,b) is achieved within the desired absolute and relative error limits, epsabs and
    /// epsrel. The function returns the final approximation, result, and an estimate of the
    /// absolute error, abserr. The integration rule is determined by the value of key, which should
    /// be chosen from the following symbolic names,
    ///
    /// GSL_INTEG_GAUSS15  (key = 1)
    ///
    /// GSL_INTEG_GAUSS21  (key = 2)
    ///
    /// GSL_INTEG_GAUSS31  (key = 3)
    ///
    /// GSL_INTEG_GAUSS41  (key = 4)
    ///
    /// GSL_INTEG_GAUSS51  (key = 5)
    ///
    /// GSL_INTEG_GAUSS61  (key = 6)
    ///
    /// corresponding to the 15f64, 21f64, 31f64, 41f64, 51 and 61 point Gauss-Kronrod rules. The
    /// higher-order rules give better accuracy for smooth functions, while lower-order rules save
    /// time when the function contains local difficulties, such as discontinuities.
    ///
    /// On each iteration the adaptive integration strategy bisects the interval with the largest
    /// error estimate. The subintervals and their results are stored in the memory provided by
    /// workspace. The maximum number of subintervals is given by limit, which may not exceed the
    /// allocated size of the workspace.
    ///
    /// Returns `(result, abs_err)` if everything went fine.
    pub fn qag<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
        key: enums::GaussKonrodRule,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qag(
                &function,
                a,
                b,
                epsabs,
                epsrel,
                limit,
                key.into(),
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }

    /// This function applies the Gauss-Kronrod 21-point integration rule adaptively until an
    /// estimate of the integral of f over (a,b) is achieved within the desired absolute and
    /// relative error limits, epsabs and epsrel. The results are extrapolated using the
    /// epsilon-algorithm, which accelerates the convergence of the integral in the presence of
    /// discontinuities and integrable singularities. The function returns the final approximation
    /// from the extrapolation, result, and an estimate of the absolute error, abserr. The
    /// subintervals and their results are stored in the memory provided by workspace. The maximum
    /// number of subintervals is given by limit, which may not exceed the allocated size of the
    /// workspace.
    ///
    /// Returns `(result, abs_err)` if everything went fine.
    pub fn qags<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qags(
                &function,
                a,
                b,
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }

    /// This function applies the adaptive integration algorithm QAGS taking account of the
    /// user-supplied locations of singular points. The array pts of length npts should contain the
    /// endpoints of the integration ranges defined by the integration region and locations of the
    /// singularities.
    ///
    /// For example, to integrate over the region (a,b) with break-points at x_1, x_2, x_3
    /// (where a < x_1 < x_2 < x_3 < b) the following pts array should be used
    ///
    /// pts[0] = a
    /// pts[1] = x_1
    /// pts[2] = x_2
    /// pts[3] = x_3
    /// pts[4] = b
    /// with npts = 5.
    ///
    /// If you know the locations of the singular points in the integration region then this routine
    /// will be faster than QAGS.
    ///
    /// Returns `(result, abs_err)` if everything went fine.
    pub fn qagp<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        pts: &mut [f64],
        epsabs: f64,
        epsrel: f64,
        limit: u64,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qagp(
                &function,
                pts.as_mut_ptr(),
                pts.len() as _,
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }

    /// This function computes the integral of the function f over the infinite interval
    /// `(-\infty,+\infty)`. The integral is mapped onto the semi-open interval `(0,1]` using the
    /// transformation:
    ///
    /// ```text
    /// x = (1-t)/t,
    ///
    /// \int_{-\infty}^{+\infty} dx f(x) =
    ///      \int_0^1 dt (f((1-t)/t) + f((-1+t)/t))/t^2.
    /// ```
    ///
    /// It is then integrated using the QAGS algorithm. The normal 21-point Gauss-Kronrod rule of
    /// QAGS is replaced by a 15-point rule, because the transformation can generate an integrable
    /// singularity at the origin. In this case a lower-order rule is more efficient.
    ///
    /// Returns `(result, abs_err)` if everything went fine.
    pub fn qagi<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qagi(
                &mut function,
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }

    /// This function computes the integral of the function f over the semi-infinite interval
    /// `(a,+\infty)`. The integral is mapped onto the semi-open interval `(0,1]` using the
    /// transformation:
    ///
    /// ```text
    /// x = a + (1-t)/t,
    ///
    /// \int_{a}^{+\infty} dx f(x) =
    ///      \int_0^1 dt f(a + (1-t)/t)/t^2
    /// ```
    ///
    /// and then integrated using the QAGS algorithm.
    ///
    /// Returns `(result, abs_err)` if everything went fine.
    pub fn qagiu<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        a: f64,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qagiu(
                &mut function,
                a,
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }

    /// This function computes the integral of the function f over the semi-infinite interval
    /// `(-\infty,b)`. The integral is mapped onto the semi-open interval `(0,1]` using the
    /// transformation:
    ///
    /// ```text
    ///  x = b - (1-t)/t,
    ///
    /// \int_{-\infty}^{b} dx f(x) =
    ///      \int_0^1 dt f(b - (1-t)/t)/t^2
    /// ```
    ///
    /// and then integrated using the QAGS algorithm.
    ///
    /// Returns `(result, abs_err)` if everything went fine.
    pub fn qagil<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        b: f64,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qagil(
                &mut function,
                b,
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }

    /// This function computes the Cauchy principal value of the integral of f over `(a,b)`, with a
    /// singularity at c,
    ///
    /// ```text
    /// I = \int_a^b dx f(x) / (x - c)
    /// ```
    ///
    /// The adaptive bisection algorithm of QAG is used, with modifications to ensure that
    /// subdivisions do not occur at the singular point x = c.
    ///
    /// When a subinterval contains the point x = c or is close to it then a special 25-point
    /// modified Clenshaw-Curtis rule is used to control the singularity. Further away from the
    /// singularity the algorithm uses an ordinary 15-point Gauss-Kronrod integration rule.
    ///
    /// Returns `(result, abs_err)` if everything went fine.
    pub fn qawc<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        a: f64,
        b: f64,
        c: f64,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qawc(
                &mut function,
                a,
                b,
                c,
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }
}

impl Drop for IntegrationWorkspace {
    fn drop(&mut self) {
        unsafe { sys::gsl_integration_workspace_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl FFI<sys::gsl_integration_workspace> for IntegrationWorkspace {
    fn wrap(w: *mut sys::gsl_integration_workspace) -> IntegrationWorkspace {
        IntegrationWorkspace { w: w }
    }

    fn soft_wrap(w: *mut sys::gsl_integration_workspace) -> IntegrationWorkspace {
        Self::wrap(w)
    }

    fn unwrap_shared(w: &IntegrationWorkspace) -> *const sys::gsl_integration_workspace {
        w.w as *const _
    }

    fn unwrap_unique(w: &mut IntegrationWorkspace) -> *mut sys::gsl_integration_workspace {
        w.w
    }
}

/// The QAWS algorithm is designed for integrands with algebraic-logarithmic singularities at the
/// end-points of an integration region. In order to work efficiently the algorithm requires a
/// precomputed table of Chebyshev moments.
pub struct IntegrationQawsTable {
    w: *mut sys::gsl_integration_qaws_table,
}

impl IntegrationQawsTable {
    /// This function allocates space for a gsl_integration_qaws_table struct describing a singular
    /// weight function W(x) with the parameters `alpha`, `beta`, `mu` and `nu`,
    ///
    /// ```text
    /// W(x) = (x-a)^alpha (b-x)^beta log^mu (x-a) log^nu (b-x)
    /// ```
    ///
    /// where `alpha > -1f64`, `beta > -1f64`, and `mu = 0, 1`, `nu = 0, 1`. The weight function can
    /// take four different forms depending on the values of `mu` and `nu`,
    ///
    /// ```text
    /// W(x) = (x-a)^alpha (b-x)^beta                   (mu = 0, nu = 0)
    /// W(x) = (x-a)^alpha (b-x)^beta log(x-a)          (mu = 1, nu = 0)
    /// W(x) = (x-a)^alpha (b-x)^beta log(b-x)          (mu = 0, nu = 1)
    /// W(x) = (x-a)^alpha (b-x)^beta log(x-a) log(b-x) (mu = 1, nu = 1)
    /// ```
    ///
    /// The singular points (a,b) do not have to be specified until the integral is computed, where
    /// they are the endpoints of the integration range.
    ///
    /// The function returns a pointer to the newly allocated table gsl_integration_qaws_table if no
    /// errors were detected, and 0 in the case of error.
    pub fn new(alpha: f64, beta: f64, mu: i32, nu: i32) -> Option<IntegrationQawsTable> {
        let tmp = unsafe { sys::gsl_integration_qaws_table_alloc(alpha, beta, mu, nu) };

        if tmp.is_null() {
            None
        } else {
            Some(IntegrationQawsTable { w: tmp })
        }
    }

    /// This function modifies the parameters (\alpha, \beta, \mu, \nu)
    pub fn set(&mut self, alpha: f64, beta: f64, mu: i32, nu: i32) -> ::Value {
        ::Value::from(unsafe { sys::gsl_integration_qaws_table_set(self.w, alpha, beta, mu, nu) })
    }

    /// This function computes the integral of the function f(x) over the interval (a,b) with the
    /// singular weight function `(x-a)^\alpha (b-x)^\beta \log^\mu (x-a) \log^\nu (b-x)`. The
    /// parameters of the weight function (\alpha, \beta, \mu, \nu) are taken from the table self.
    /// The integral is,
    ///
    /// ```text
    /// I = \int_a^b dx f(x) (x-a)^alpha (b-x)^beta log^mu (x-a) log^nu (b-x).
    /// ```
    ///
    /// The adaptive bisection algorithm of QAG is used. When a subinterval contains one of the
    /// endpoints then a special 25-point modified Clenshaw-Curtis rule is used to control the
    /// singularities. For subintervals which do not include the endpoints an ordinary 15-point
    /// Gauss-Kronrod integration rule is used.
    ///
    /// Returns `(result, abs_err)`
    pub fn qaws<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
        workspace: &mut IntegrationWorkspace,
    ) -> Result<(f64, f64), ::Value> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut function = wrap_callback!(f, F);

        let ret = unsafe {
            sys::gsl_integration_qaws(
                &mut function,
                a,
                b,
                FFI::unwrap_unique(self),
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(workspace),
                &mut result,
                &mut abs_err,
            )
        };
        result!(ret, (result, abs_err))
    }
}

impl Drop for IntegrationQawsTable {
    fn drop(&mut self) {
        unsafe { sys::gsl_integration_qaws_table_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl FFI<sys::gsl_integration_qaws_table> for IntegrationQawsTable {
    fn wrap(w: *mut sys::gsl_integration_qaws_table) -> IntegrationQawsTable {
        IntegrationQawsTable { w: w }
    }

    fn soft_wrap(w: *mut sys::gsl_integration_qaws_table) -> IntegrationQawsTable {
        Self::wrap(w)
    }

    fn unwrap_shared(w: &IntegrationQawsTable) -> *const sys::gsl_integration_qaws_table {
        w.w as *const _
    }

    fn unwrap_unique(w: &mut IntegrationQawsTable) -> *mut sys::gsl_integration_qaws_table {
        w.w
    }
}

/// The QAWO algorithm is designed for integrands with an oscillatory factor, \sin(\omega x) or
/// \cos(\omega x). In order to work efficiently the algorithm requires a table of Chebyshev moments
/// which must be pre-computed with calls to the functions below.
pub struct IntegrationQawoTable {
    w: *mut sys::gsl_integration_qawo_table,
}

impl IntegrationQawoTable {
    /// This function allocates space for a gsl_integration_qawo_table struct and its associated
    /// workspace describing a sine or cosine weight function W(x) with the parameters (\omega, L),
    ///
    /// ```text
    /// W(x) = sin(omega x)
    /// W(x) = cos(omega x)
    /// ```
    ///
    /// The parameter L must be the length of the interval over which the function will be
    /// integrated L = b - a. The choice of sine or cosine is made with the parameter sine which
    /// should be chosen from one of the two following symbolic values:
    ///
    /// ```text
    /// ::Cosine
    /// ::IntegrationQawo::Sine
    /// ```
    ///
    /// The gsl_integration_qawo_table is a table of the trigonometric coefficients required in the
    /// integration process. The parameter n determines the number of levels of coefficients that
    /// are computed. Each level corresponds to one bisection of the interval L, so that n levels
    /// are sufficient for subintervals down to the length L/2^n. The integration routine
    /// gsl_integration_qawo returns the error ::Table if the number of levels is insufficient for
    /// the requested accuracy.
    pub fn new(
        omega: f64,
        l: f64,
        sine: ::IntegrationQawo,
        n: u64,
    ) -> Option<IntegrationQawoTable> {
        let tmp = unsafe { sys::gsl_integration_qawo_table_alloc(omega, l, sine.into(), n) };

        if tmp.is_null() {
            None
        } else {
            Some(IntegrationQawoTable { w: tmp })
        }
    }

    /// This function changes the parameters omega, L and sine of the existing self workspace.
    pub fn set(&mut self, omega: f64, l: f64, sine: ::IntegrationQawo) -> ::Value {
        ::Value::from(unsafe { sys::gsl_integration_qawo_table_set(self.w, omega, l, sine.into()) })
    }

    /// This function allows the length parameter l of the self workspace to be changed.
    pub fn set_length(&mut self, l: f64) -> ::Value {
        ::Value::from(unsafe { sys::gsl_integration_qawo_table_set_length(self.w, l) })
    }

    /// This function uses an adaptive algorithm to compute the integral of f over (a,b) with the
    /// weight function \sin(\omega x) or \cos(\omega x) defined by the table `wf`,
    ///
    /// I = \int_a^b dx f(x) sin(omega x)
    /// I = \int_a^b dx f(x) cos(omega x)
    ///
    /// The results are extrapolated using the epsilon-algorithm to accelerate the convergence of
    /// the integral. The function returns the final approximation from the extrapolation, result,
    /// and an estimate of the absolute error, abserr. The subintervals and their results are
    /// stored in the memory provided by workspace. The maximum number of subintervals is given by
    /// limit, which may not exceed the allocated size of the workspace.
    ///
    /// Those subintervals with “large” widths d where d\omega > 4 are computed using a 25-point
    /// Clenshaw-Curtis integration rule, which handles the oscillatory behavior. Subintervals with
    /// a "small" widths where d\omega < 4 are computed using a 15-point Gauss-Kronrod integration.
    ///
    /// Returns `(result, abserr)` if everything went fine.
    pub fn qawo<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        a: f64,
        epsabs: f64,
        epsrel: f64,
        limit: u64,
        workspace: &mut IntegrationWorkspace,
    ) -> Result<(f64, f64), ::Value> {
        let mut function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abserr = 0.;

        let ret = unsafe {
            sys::gsl_integration_qawo(
                &mut function,
                a,
                epsabs,
                epsrel,
                limit,
                FFI::unwrap_unique(workspace),
                FFI::unwrap_unique(self),
                &mut result,
                &mut abserr,
            )
        };
        result!(ret, (result, abserr))
    }
}

impl Drop for IntegrationQawoTable {
    fn drop(&mut self) {
        unsafe { sys::gsl_integration_qawo_table_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl FFI<sys::gsl_integration_qawo_table> for IntegrationQawoTable {
    fn wrap(w: *mut sys::gsl_integration_qawo_table) -> IntegrationQawoTable {
        IntegrationQawoTable { w: w }
    }

    fn soft_wrap(w: *mut sys::gsl_integration_qawo_table) -> IntegrationQawoTable {
        Self::wrap(w)
    }

    fn unwrap_shared(w: &IntegrationQawoTable) -> *const sys::gsl_integration_qawo_table {
        w.w as *const _
    }

    fn unwrap_unique(w: &mut IntegrationQawoTable) -> *mut sys::gsl_integration_qawo_table {
        w.w
    }
}

/// CQUAD is a new doubly-adaptive general-purpose quadrature routine which can handle most types of
/// singularities, non-numerical function values such as Inf or NaN, as well as some divergent
/// integrals. It generally requires more function evaluations than the integration routines in
/// QUADPACK, yet fails less often for difficult integrands.
///
/// The underlying algorithm uses a doubly-adaptive scheme in which Clenshaw-Curtis quadrature rules
/// of increasing degree are used to compute the integral in each interval. The L_2-norm of the
/// difference between the underlying interpolatory polynomials of two successive rules is used as
/// an error estimate. The interval is subdivided if the difference between two successive rules is
/// too large or a rule of maximum degree has been reached.
pub struct CquadWorkspace {
    w: *mut sys::gsl_integration_cquad_workspace,
}

impl CquadWorkspace {
    /// This function allocates a workspace sufficient to hold the data for n intervals. The number
    /// n is not the maximum number of intervals that will be evaluated. If the workspace is full,
    /// intervals with smaller error estimates will be discarded. A minimum of 3 intervals
    /// is required and for most functions, a workspace of size 100 is sufficient.
    pub fn new(n: u64) -> Option<CquadWorkspace> {
        let tmp = unsafe { sys::gsl_integration_cquad_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(CquadWorkspace { w: tmp })
        }
    }

    /// This function computes the integral of f over (a,b) within the desired absolute and relative
    /// error limits, epsabs and epsrel using the CQUAD algorithm. The function returns the final
    /// approximation, result, an estimate of the absolute error, abserr, and the number of function
    /// evaluations required, nevals.
    ///
    /// The CQUAD algorithm divides the integration region into subintervals, and in each iteration,
    /// the subinterval with the largest estimated error is processed. The algorithm uses
    /// Clenshaw-Curits quadrature rules of degree 4, 8, 16 and 32 over 5, 9, 17 and 33 nodes
    /// respectively. Each interval is initialized with the lowest-degree rule. When an interval is
    /// processed, the next-higher degree rule is evaluated and an error estimate is computed based
    /// on the L_2-norm of the difference between the underlying interpolating polynomials of both
    /// rules. If the highest-degree rule has already been used, or the interpolatory polynomials
    /// differ significantly, the interval is bisected.
    ///
    /// The subintervals and their results are stored in the memory provided by workspace. If the
    /// error estimate or the number of function evaluations is not needed, the pointers abserr and
    /// nevals can be set to NULL (not in rgsl).
    ///
    /// Returns `(result, abs_err, n_evals)` if everything went fine.
    pub fn cquad<F: Fn(f64) -> f64>(
        &mut self,
        f: F,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
    ) -> Result<(f64, f64, u64), ::Value> {
        let mut function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut n_evals = 0;

        let ret = unsafe {
            sys::gsl_integration_cquad(
                &mut function,
                a,
                b,
                epsabs,
                epsrel,
                FFI::unwrap_unique(self),
                &mut result,
                &mut abs_err,
                &mut n_evals,
            )
        };
        result!(ret, (result, abs_err, n_evals))
    }
}

impl Drop for CquadWorkspace {
    fn drop(&mut self) {
        unsafe { sys::gsl_integration_cquad_workspace_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl FFI<sys::gsl_integration_cquad_workspace> for CquadWorkspace {
    fn wrap(w: *mut sys::gsl_integration_cquad_workspace) -> CquadWorkspace {
        CquadWorkspace { w: w }
    }

    fn soft_wrap(w: *mut sys::gsl_integration_cquad_workspace) -> CquadWorkspace {
        Self::wrap(w)
    }

    fn unwrap_shared(w: &CquadWorkspace) -> *const sys::gsl_integration_cquad_workspace {
        w.w as *const _
    }

    fn unwrap_unique(w: &mut CquadWorkspace) -> *mut sys::gsl_integration_cquad_workspace {
        w.w
    }
}

/// The fixed-order Gauss-Legendre integration routines are provided for fast integration of smooth
/// functions with known polynomial order. The n-point Gauss-Legendre rule is exact for polynomials
/// of order 2*n-1 or less. For example, these rules are useful when integrating basis functions to
/// form mass matrices for the Galerkin method. Unlike other numerical integration routines within
/// the library, these routines do not accept absolute or relative error bounds.
pub struct GLFixedTable {
    w: *mut sys::gsl_integration_glfixed_table,
}

impl GLFixedTable {
    /// This function determines the Gauss-Legendre abscissae and weights necessary for an n-point
    /// fixed order integration scheme. If possible, high precision precomputed coefficients are
    /// used. If precomputed weights are not available, lower precision coefficients are computed
    /// on the fly.
    pub fn new(n: u64) -> Option<GLFixedTable> {
        let tmp = unsafe { sys::gsl_integration_glfixed_table_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(GLFixedTable { w: tmp })
        }
    }

    /// For i in [0, …, t->n - 1], this function obtains the i-th Gauss-Legendre point xi and weight
    /// wi on the interval [a,b]. The points and weights are ordered by increasing point value. A
    /// function f may be integrated on [a,b] by summing wi * f(xi) over i.
    pub fn point(&self, a: f64, b: f64, i: u64, xi: &mut f64, wi: &mut f64) -> ::Value {
        ::Value::from(unsafe { sys::gsl_integration_glfixed_point(a, b, i, xi, wi, self.w) })
    }

    /// This function applies the Gauss-Legendre integration rule contained in table self and
    /// returns the result.
    pub fn glfixed<F: Fn(f64) -> f64>(&self, f: F, a: f64, b: f64) -> f64 {
        let function = wrap_callback!(f, F);
        unsafe { sys::gsl_integration_glfixed(&function, a, b, self.w) }
    }

    pub fn glfixed_point(&self, a: f64, b: f64, xi: &mut [f64], wi: &mut [f64]) -> enums::Value {
        assert!(xi.len() == wi.len());

        enums::Value::from(unsafe {
            sys::gsl_integration_glfixed_point(
                a,
                b,
                xi.len() as _,
                xi.as_mut_ptr(),
                wi.as_mut_ptr(),
                self.w,
            )
        })
    }
}

impl Drop for GLFixedTable {
    fn drop(&mut self) {
        unsafe { sys::gsl_integration_glfixed_table_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl FFI<sys::gsl_integration_glfixed_table> for GLFixedTable {
    fn wrap(w: *mut sys::gsl_integration_glfixed_table) -> GLFixedTable {
        GLFixedTable { w: w }
    }

    fn soft_wrap(w: *mut sys::gsl_integration_glfixed_table) -> GLFixedTable {
        Self::wrap(w)
    }

    fn unwrap_shared(w: &GLFixedTable) -> *const sys::gsl_integration_glfixed_table {
        w.w as *const _
    }

    fn unwrap_unique(w: &mut GLFixedTable) -> *mut sys::gsl_integration_glfixed_table {
        w.w
    }
}
