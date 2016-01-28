//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The error function is described in Abramowitz & Stegun, Chapter 7.

use std::mem::zeroed;
use enums;
use ffi;

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
pub fn erf(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_erf(x) }
}

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
pub fn erf_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_erf_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
pub fn erfc(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_erfc(x) }
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
pub fn erfc_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_erfc_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
pub fn log_erfc(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_log_erfc(x) }
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
pub fn log_erfc_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_log_erfc_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
pub fn erf_Z(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_erf_Z(x) }
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
pub fn erf_Z_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_erf_Z_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the upper tail of the Gaussian probability function Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2).
/// 
/// The hazard function for the normal distribution, also known as the inverse Mills’ ratio, is defined as,
/// 
/// h(x) = Z(x)/Q(x) = \sqrt{2/\pi} \exp(-x^2 / 2) / \erfc(x/\sqrt 2)
/// 
/// It decreases rapidly as x approaches -\infty and asymptotes to h(x) \sim x as x approaches +\infty.
pub fn erf_Q(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_erf_Q(x) }
}

/// This routine computes the upper tail of the Gaussian probability function Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2).
/// 
/// The hazard function for the normal distribution, also known as the inverse Mills’ ratio, is defined as,
/// 
/// h(x) = Z(x)/Q(x) = \sqrt{2/\pi} \exp(-x^2 / 2) / \erfc(x/\sqrt 2)
/// 
/// It decreases rapidly as x approaches -\infty and asymptotes to h(x) \sim x as x approaches +\infty.
pub fn erf_Q_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_erf_Q_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the hazard function for the normal distribution.
pub fn hazard(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_hazard(x) }
}

/// This routine computes the hazard function for the normal distribution.
pub fn hazard_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_hazard_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

pub fn str_error(error: ::Value) -> &'static str {
    match error {
        ::Value::Success => "Success",
        ::Value::Failure => "Failure",
        ::Value::Continue => "The iteration has not converged yet",
        ::Value::Dom => "Input domain error",
        ::Value::Range => "Output range error",
        ::Value::Fault => "Invalid pointer",
        ::Value::Inval => "Invalid argument supplied by user",
        ::Value::Failed => "generic failure",
        ::Value::Factor => "Factorization failed",
        ::Value::Sanity => "Sanity check failed - shouldn't happen",
        ::Value::NoMem => "Malloc failed",
        ::Value::BadFunc => "Problem with user-supplied function",
        ::Value::RunAway => "Iterative process is out of control",
        ::Value::MaxIter => "Exceeded max number of iterations",
        ::Value::ZeroDiv => "Tried to divide by zero",
        ::Value::BadTol => "Specified tolerance is invalid or theoretically unattainable",
        ::Value::Tol => "Failed to reach the specified tolerance",
        ::Value::UndrFlw => "Underflow",
        ::Value::OvrFlw => "Overflow",
        ::Value::Loss => "Loss of accuracy",
        ::Value::Round => "Roundoff error",
        ::Value::BadLen => "Matrix/vector sizes are not conformant",
        ::Value::NotSqr => "Matrix not square",
        ::Value::Sing => "Singularity or extremely bad function behavior detected",
        ::Value::Diverge => "Integral or series is divergent",
        ::Value::Unsup => "The required feature is not supported by this hardware platform",
        ::Value::Unimpl => "The requested feature is not (yet) implemented",
        ::Value::Cache => "Cache limit exceeded",
        ::Value::Table => "Table limit exceeded",
        ::Value::NoProg => "Iteration is not making progress towards solution",
        ::Value::NoProgJ => "Jacobian evaluations are not improving the solution",
        ::Value::TolF => "Cannot reach the specified tolerance in F",
        ::Value::TolX => "Cannot reach the specified tolerance in X",
        ::Value::TolG => "Cannot reach the specified tolerance in gradient",
        ::Value::EOF => "End of file"
    }
}