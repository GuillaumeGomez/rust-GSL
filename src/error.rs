//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The error function is described in Abramowitz & Stegun, Chapter 7.

use crate::{types, Value};
use std::ffi::CStr;
use std::mem::MaybeUninit;
use std::os::raw::{c_char, c_int};

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
#[doc(alias = "gsl_sf_erf")]
pub fn erf(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erf(x) }
}

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
#[doc(alias = "gsl_sf_erf_e")]
pub fn erf_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erf_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
#[doc(alias = "gsl_sf_erfc")]
pub fn erfc(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erfc(x) }
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
#[doc(alias = "gsl_sf_erfc_e")]
pub fn erfc_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erfc_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
#[doc(alias = "gsl_sf_log_erfc")]
pub fn log_erfc(x: f64) -> f64 {
    unsafe { sys::gsl_sf_log_erfc(x) }
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
#[doc(alias = "gsl_sf_log_erfc_e")]
pub fn log_erfc_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_log_erfc_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
#[doc(alias = "gsl_sf_erf_Z")]
pub fn erf_Z(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erf_Z(x) }
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
#[doc(alias = "gsl_sf_erf_Z_e")]
pub fn erf_Z_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erf_Z_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the upper tail of the Gaussian probability function Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2).
///
/// The hazard function for the normal distribution, also known as the inverse Mills’ ratio, is defined as,
///
/// h(x) = Z(x)/Q(x) = \sqrt{2/\pi} \exp(-x^2 / 2) / \erfc(x/\sqrt 2)
///
/// It decreases rapidly as x approaches -\infty and asymptotes to h(x) \sim x as x approaches +\infty.
#[doc(alias = "gsl_sf_erf_Q")]
pub fn erf_Q(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erf_Q(x) }
}

/// This routine computes the upper tail of the Gaussian probability function Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2).
///
/// The hazard function for the normal distribution, also known as the inverse Mills’ ratio, is defined as,
///
/// h(x) = Z(x)/Q(x) = \sqrt{2/\pi} \exp(-x^2 / 2) / \erfc(x/\sqrt 2)
///
/// It decreases rapidly as x approaches -\infty and asymptotes to h(x) \sim x as x approaches +\infty.
#[doc(alias = "gsl_sf_erf_Q_e")]
pub fn erf_Q_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erf_Q_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the hazard function for the normal distribution.
#[doc(alias = "gsl_sf_hazard")]
pub fn hazard(x: f64) -> f64 {
    unsafe { sys::gsl_sf_hazard(x) }
}

/// This routine computes the hazard function for the normal distribution.
#[doc(alias = "gsl_sf_hazard_e")]
pub fn hazard_e(x: f64) -> Result<types::Result, Value> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hazard_e(x, result.as_mut_ptr()) };

    result_handler!(ret, unsafe { result.assume_init() }.into())
}

pub fn str_error(error: crate::Value) -> &'static str {
    match error {
        crate::Value::Success => "Success",
        crate::Value::Failure => "Failure",
        crate::Value::Continue => "The iteration has not converged yet",
        crate::Value::Domain => "Input domain error",
        crate::Value::Range => "Output range error",
        crate::Value::Fault => "Invalid pointer",
        crate::Value::Invalid => "Invalid argument supplied by user",
        crate::Value::Failed => "generic failure",
        crate::Value::Factorization => "Factorization failed",
        crate::Value::Sanity => "Sanity check failed - shouldn't happen",
        crate::Value::NoMemory => "Malloc failed",
        crate::Value::BadFunction => "Problem with user-supplied function",
        crate::Value::RunAway => "Iterative process is out of control",
        crate::Value::MaxIteration => "Exceeded max number of iterations",
        crate::Value::ZeroDiv => "Tried to divide by zero",
        crate::Value::BadTolerance => {
            "Specified tolerance is invalid or theoretically unattainable"
        }
        crate::Value::Tolerance => "Failed to reach the specified tolerance",
        crate::Value::UnderFlow => "Underflow",
        crate::Value::OverFlow => "Overflow",
        crate::Value::Loss => "Loss of accuracy",
        crate::Value::Round => "Roundoff error",
        crate::Value::BadLength => "Matrix/vector sizes are not conformant",
        crate::Value::NotSquare => "Matrix not square",
        crate::Value::Singularity => "Singularity or extremely bad function behavior detected",
        crate::Value::Diverge => "Integral or series is divergent",
        crate::Value::Unsupported => {
            "The required feature is not supported by this hardware platform"
        }
        crate::Value::Unimplemented => "The requested feature is not (yet) implemented",
        crate::Value::Cache => "Cache limit exceeded",
        crate::Value::Table => "Table limit exceeded",
        crate::Value::NoProgress => "Iteration is not making progress towards solution",
        crate::Value::NoProgressJacobian => "Jacobian evaluations are not improving the solution",
        crate::Value::ToleranceF => "Cannot reach the specified tolerance in F",
        crate::Value::ToleranceX => "Cannot reach the specified tolerance in X",
        crate::Value::ToleranceG => "Cannot reach the specified tolerance in gradient",
        crate::Value::EOF => "End of file",
        crate::Value::Unknown(_) => "Unknown error",
    }
}

static mut CALLBACK: Option<fn(&str, &str, u32, crate::Value)> = None;

/// `f` is the type of GSL error handler functions. An error handler will be passed four arguments
/// which specify the reason for the error (a string), the name of the source file in which it
/// occurred (also a string), the line number in that file (an integer) and the error number (an
/// integer). The source file and line number are set at compile time using the __FILE__ and
/// __LINE__ directives in the preprocessor. An error handler function returns type void. Error
/// handler functions should be defined like this,
///
/// This function sets a new error handler, new_handler, for the GSL library routines. The previous
/// handler is returned (so that you can restore it later). Note that the pointer to a user defined
/// error handler function is stored in a static variable, so there can be only one error handler
/// per program. This function should be not be used in multi-threaded programs except to set up a
/// program-wide error handler from a master thread. The following example shows how to set and
/// restore a new error handler,
///
/// ```
/// use crate::rgsl::error::set_error_handler;
/// use crate::rgsl::Value;
///
/// fn error_handling(error_str: &str, file: &str, line: u32, error_value: Value) {
///     println!("[{:?}] '{}:{}': {}", error_value, file, line, error_str);
/// }
///
/// /* save original handler, install new handler */
/// let old_handler = set_error_handler(Some(error_handling));
///
/// /* code uses new handler */
/// // ...
///
/// /* restore original handler */
/// set_error_handler(old_handler);
/// ```
///
/// To use the default behavior (abort on error) set the error handler to NULL,
///
/// ```
/// # use crate::rgsl::error::set_error_handler;
/// let old_handler = set_error_handler(None);
/// ```
#[doc(alias = "gsl_set_error_handler")]
pub fn set_error_handler(
    f: Option<fn(&str, &str, u32, crate::Value)>,
) -> Option<fn(&str, &str, u32, crate::Value)> {
    unsafe {
        let out = CALLBACK.take();
        match f {
            Some(f) => {
                CALLBACK = Some(f);
                sys::gsl_set_error_handler(Some(inner_error_handler));
            }
            None => {
                sys::gsl_set_error_handler(None);
            }
        }
        out
    }
}

/// This function turns off the error handler by defining an error handler which does nothing. This
/// will cause the program to continue after any error, so the return values from any library
/// routines must be checked. This is the recommended behavior for production programs. The previous
/// handler is returned (so that you can restore it later).
#[doc(alias = "gsl_set_error_handler_off")]
pub fn set_error_handler_off() -> Option<fn(&str, &str, u32, crate::Value)> {
    unsafe {
        sys::gsl_set_error_handler_off();
        CALLBACK.take()
    }
}

extern "C" fn inner_error_handler(
    reason: *const c_char,
    file: *const c_char,
    line: c_int,
    gsl_errno: c_int,
) {
    unsafe {
        if let Some(ref call) = CALLBACK {
            let s = CStr::from_ptr(reason);
            let f = CStr::from_ptr(file);
            call(
                s.to_str().unwrap_or("Unknown"),
                f.to_str().unwrap_or("Unknown"),
                line as _,
                crate::Value::from(gsl_errno),
            );
        }
    }
}

#[test]
fn test_error_handler() {
    use crate::{bessel, Value};

    set_error_handler_off();
    match bessel::K0_e(1e3) {
        Err(Value::UnderFlow) => println!("K0(1e3) underflowed"),
        _ => panic!("unexpected"),
    }
}
