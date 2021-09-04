//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
## Introduction

Each algorithm computes an approximation to a definite integral of the form,

I = \int_a^b f(x) w(x) dx
where w(x) is a weight function (for general integrands w(x)=1). The user provides absolute and relative error bounds (epsabs, epsrel) which
specify the following accuracy requirement,

|RESULT - I|  <= max(epsabs, epsrel |I|)

where RESULT is the numerical approximation obtained by the algorithm. The algorithms attempt to estimate the absolute error ABSERR = |RESULT
- I| in such a way that the following inequality holds,

|RESULT - I| <= ABSERR <= max(epsabs, epsrel |I|)

In short, the routines return the first approximation which has an absolute error smaller than epsabs or a relative error smaller than epsrel.

Note that this is an either-or constraint, not simultaneous. To compute to a specified absolute error, set epsrel to zero. To compute to a
specified relative error, set epsabs to zero. The routines will fail to converge if the error bounds are too stringent, but always return the
best approximation obtained up to that stage.

The algorithms in QUADPACK use a naming convention based on the following letters,

Q - quadrature routine

N - non-adaptive integrator
A - adaptive integrator

G - general integrand (user-defined)
W - weight function with integrand

S - singularities can be more readily integrated
P - points of special difficulty can be supplied
I - infinite range of integration
O - oscillatory weight function, cos or sin
F - Fourier integral
C - Cauchy principal value
The algorithms are built on pairs of quadrature rules, a higher order rule and a lower order rule. The higher order rule is used to compute the
best approximation to an integral over a small range. The difference between the results of the higher order rule and the lower order rule gives
an estimate of the error in the approximation.

 * [Integrands without weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-without-weight-functions.html#Integrands-without-weight-functions)
 * [Integrands with weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-with-weight-functions.html#Integrands-with-weight-functions)
 * [Integrands with singular weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-with-singular-weight-functions.html#Integrands-with-singular-weight-functions)

## QNG non-adaptive Gauss-Kronrod integration

The QNG algorithm is a non-adaptive procedure which uses fixed Gauss-Kronrod-Patterson abscissae to sample the integrand at a maximum of 87
points. It is provided for fast integration of smooth functions.

## QAG adaptive integration

The QAG algorithm is a simple adaptive integration procedure. The integration region is divided into subintervals, and on each iteration the
subinterval with the largest estimated error is bisected. This reduces the overall error rapidly, as the subintervals become concentrated
around local difficulties in the integrand. These subintervals are managed by a gsl_integration_workspace struct, which handles the memory
for the subinterval ranges, results and error estimates.

## QAGS adaptive integration with singularities

The presence of an integrable singularity in the integration region causes an adaptive routine to concentrate new subintervals around the
singularity. As the subintervals decrease in size the successive approximations to the integral converge in a limiting fashion. This
approach to the limit can be accelerated using an extrapolation procedure. The QAGS algorithm combines adaptive bisection with the Wynn
epsilon-algorithm to speed up the integration of many types of integrable singularities.

## References and Further Reading

The following book is the definitive reference for QUADPACK, and was written by the original authors. It provides descriptions of the
algorithms, program listings, test programs and examples. It also includes useful advice on numerical integration and many references
to the numerical integration literature used in developing QUADPACK.

R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K. Kahaner. QUADPACK A subroutine package for automatic integration Springer Verlag, 1983.
The CQUAD integration algorithm is described in the following paper:

P. Gonnet, “Increasing the Reliability of Adaptive Quadrature Using Explicit Interpolants”, ACM Transactions on Mathematical Software, Volume 37
(2010), Issue 3, Article 26.
!*/

use crate::Value;
use ffi::FFI;

/// This function applies the Gauss-Kronrod 10-point, 21-point, 43-point and 87-point integration
/// rules in succession until an estimate of the integral of f over (a,b) is achieved within the
/// desired absolute and relative error limits, eps_abs and eps_rel. The function returns the final
/// approximation, result, an estimate of the absolute error, abserr and the number of function
/// evaluations used, neval. The Gauss-Kronrod rules are designed in such a way that each rule uses
/// all the results of its predecessors, in order to minimize the total number of function
/// evaluations.
///
/// Returns `(result, abs_err, n_eval)`.
#[doc(alias = "gsl_integration_qng")]
pub fn qng<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    eps_abs: f64,
    eps_rel: f64,
) -> (::Value, f64, f64, usize) {
    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut n_eval = 0;

    let ret = unsafe {
        sys::gsl_integration_qng(
            &function,
            a,
            b,
            eps_abs,
            eps_rel,
            &mut result,
            &mut abs_err,
            &mut n_eval,
        )
    };
    (::Value::from(ret), result, abs_err, n_eval)
}

/// Gauss quadrature weights and kronrod quadrature abscissae and weights as evaluated with 80
/// decimal digit arithmetic by L. W.
///
/// Fullerton, Bell Labs, Nov. 1981.
///
/// Returns `(result, abs_err, resabs, resasc)`.
#[doc(alias = "gsl_integration_qk15")]
pub fn qk15<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut resabs = 0.;
    let mut resasc = 0.;

    unsafe {
        sys::gsl_integration_qk15(
            &function,
            a,
            b,
            &mut result,
            &mut abs_err,
            &mut resabs,
            &mut resasc,
        )
    };
    (result, abs_err, resabs, resasc)
}

/// Returns `(result, abs_err, resabs, resasc)`.
#[doc(alias = "gsl_integration_qk21")]
pub fn qk21<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut resabs = 0.;
    let mut resasc = 0.;

    unsafe {
        sys::gsl_integration_qk21(
            &function,
            a,
            b,
            &mut result,
            &mut abs_err,
            &mut resabs,
            &mut resasc,
        )
    };
    (result, abs_err, resabs, resasc)
}

/// Returns `(result, abs_err, resabs, resasc)`.
#[doc(alias = "gsl_integration_qk31")]
pub fn qk31<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut resabs = 0.;
    let mut resasc = 0.;

    unsafe {
        sys::gsl_integration_qk31(
            &function,
            a,
            b,
            &mut result,
            &mut abs_err,
            &mut resabs,
            &mut resasc,
        )
    };
    (result, abs_err, resabs, resasc)
}

/// Returns `(result, abs_err, resabs, resasc)`.
#[doc(alias = "gsl_integration_qk41")]
pub fn qk41<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut resabs = 0.;
    let mut resasc = 0.;

    unsafe {
        sys::gsl_integration_qk41(
            &function,
            a,
            b,
            &mut result,
            &mut abs_err,
            &mut resabs,
            &mut resasc,
        )
    };
    (result, abs_err, resabs, resasc)
}

/// Returns `(result, abs_err, resabs, resasc)`.
#[doc(alias = "gsl_integration_qk51")]
pub fn qk51<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut resabs = 0.;
    let mut resasc = 0.;

    unsafe {
        sys::gsl_integration_qk51(
            &function,
            a,
            b,
            &mut result,
            &mut abs_err,
            &mut resabs,
            &mut resasc,
        )
    };
    (result, abs_err, resabs, resasc)
}

/// Returns `(result, abs_err, resabs, resasc)`.
#[doc(alias = "gsl_integration_qk61")]
pub fn qk61<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut resabs = 0.;
    let mut resasc = 0.;

    unsafe {
        sys::gsl_integration_qk61(
            &function,
            a,
            b,
            &mut result,
            &mut abs_err,
            &mut resabs,
            &mut resasc,
        )
    };
    (result, abs_err, resabs, resasc)
}

/// Returns `(result, abs_err, resabs, resasc)`.
#[doc(alias = "gsl_integration_qk")]
pub fn qk<F: Fn(f64) -> f64>(
    xgk: &[f64],
    wg: &[f64],
    wgk: &[f64],
    fv1: &mut [f64],
    fv2: &mut [f64],
    f: F,
    a: f64,
    b: f64,
) -> (f64, f64, f64, f64) {
    assert!(xgk.len() == wg.len());
    assert!(xgk.len() == wgk.len());
    assert!(xgk.len() == fv1.len());
    assert!(xgk.len() == fv2.len());

    let function = wrap_callback!(f, F);
    let mut result = 0.;
    let mut abs_err = 0.;
    let mut resabs = 0.;
    let mut resasc = 0.;

    unsafe {
        sys::gsl_integration_qk(
            xgk.len() as _,
            xgk.as_ptr(),
            wg.as_ptr(),
            wgk.as_ptr(),
            fv1.as_mut_ptr(),
            fv2.as_mut_ptr(),
            &function,
            a,
            b,
            &mut result,
            &mut abs_err,
            &mut resabs,
            &mut resasc,
        );
    }
    (result, abs_err, resabs, resasc)
}

/// This function attempts to compute a Fourier integral of the function f over the semi-infinite
/// interval `[a,+\infty)`.
///
/// ```text
/// I = \int_a^{+\infty} dx f(x) sin(omega x)
/// I = \int_a^{+\infty} dx f(x) cos(omega x)
/// ```
///
/// The parameter \omega and choice of \sin or \cos is taken from the table wf (the length L can
/// take any value, since it is overridden by this function to a value appropriate for the Fourier
/// integration). The integral is computed using the QAWO algorithm over each of the subintervals,
///
/// ```text
/// C_1 = [a, a + c]
/// C_2 = [a + c, a + 2 c]
/// ... = ...
/// C_k = [a + (k-1) c, a + k c]
/// ```
///
/// where c = (2 floor(|\omega|) + 1) \pi/|\omega|. The width c is chosen to cover an odd number of
/// periods so that the contributions from the intervals alternate in sign and are monotonically
/// decreasing when f is positive and monotonically decreasing. The sum of this sequence of
/// contributions is accelerated using the epsilon-algorithm.
///
/// This function works to an overall absolute tolerance of abserr. The following strategy is used:
/// on each interval C_k the algorithm tries to achieve the tolerance
///
/// ```text
/// TOL_k = u_k abserr
/// ```
///
/// where u_k = (1 - p)p^{k-1} and p = 9/10. The sum of the geometric series of contributions from
/// each interval gives an overall tolerance of abserr.
///
/// If the integration of a subinterval leads to difficulties then the accuracy requirement for
/// subsequent intervals is relaxed,
///
/// ```text
/// TOL_k = u_k max(abserr, max_{i<k}{E_i})
/// ```
///
/// where E_k is the estimated error on the interval C_k.
///
/// The subintervals and their results are stored in the memory provided by workspace. The maximum
/// number of subintervals is given by limit, which may not exceed the allocated size of the
/// workspace. The integration over each subinterval uses the memory provided by cycle_workspace
/// as workspace for the QAWO algorithm.
///
/// Returns `(result, abs_err)`.
#[doc(alias = "gsl_integration_qawf")]
pub fn qawf<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    epsabs: f64,
    limit: usize,
    workspace: &mut ::IntegrationWorkspace,
    cycle_workspace: &mut ::IntegrationWorkspace,
    wf: &mut ::IntegrationQawoTable,
) -> (Value, f64, f64) {
    let mut result = 0.;
    let mut abs_err = 0.;

    let mut function = wrap_callback!(f, F);
    let ret = unsafe {
        sys::gsl_integration_qawf(
            &mut function,
            a,
            epsabs,
            limit,
            FFI::unwrap_unique(workspace),
            FFI::unwrap_unique(cycle_workspace),
            FFI::unwrap_unique(wf),
            &mut result,
            &mut abs_err,
        )
    };
    (::Value::from(ret), result, abs_err)
}
