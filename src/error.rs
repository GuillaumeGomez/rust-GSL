//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The error function is described in Abramowitz & Stegun, Chapter 7.

use std::mem::zeroed;
use enums;
use ffi;
use libc;
use std::ffi::CStr;

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
pub fn erf(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_erf(x) }
}

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
pub fn erf_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_erf_e(x, &mut result) };

    (enums::Value::from(ret), ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
pub fn erfc(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_erfc(x) }
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
pub fn erfc_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_erfc_e(x, &mut result) };

    (enums::Value::from(ret), ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
pub fn log_erfc(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_log_erfc(x) }
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
pub fn log_erfc_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_log_erfc_e(x, &mut result) };

    (enums::Value::from(ret), ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
pub fn erf_Z(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_erf_Z(x) }
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
pub fn erf_Z_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_erf_Z_e(x, &mut result) };

    (enums::Value::from(ret), ::types::Result{val: result.val, err: result.err})
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

    (enums::Value::from(ret), ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the hazard function for the normal distribution.
pub fn hazard(x: f64) -> f64 {
    unsafe { ::ffi::gsl_sf_hazard(x) }
}

/// This routine computes the hazard function for the normal distribution.
pub fn hazard_e(x: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_hazard_e(x, &mut result) };

    (enums::Value::from(ret), ::types::Result{val: result.val, err: result.err})
}

pub fn str_error(error: ::Value) -> &'static str {
    match error {
        ::Value::Success => "Success",
        ::Value::Failure => "Failure",
        ::Value::Continue => "The iteration has not converged yet",
        ::Value::Domain => "Input domain error",
        ::Value::Range => "Output range error",
        ::Value::Fault => "Invalid pointer",
        ::Value::Invalid => "Invalid argument supplied by user",
        ::Value::Failed => "generic failure",
        ::Value::Factorization => "Factorization failed",
        ::Value::Sanity => "Sanity check failed - shouldn't happen",
        ::Value::NoMemory => "Malloc failed",
        ::Value::BadFunction => "Problem with user-supplied function",
        ::Value::RunAway => "Iterative process is out of control",
        ::Value::MaxIteration => "Exceeded max number of iterations",
        ::Value::ZeroDiv => "Tried to divide by zero",
        ::Value::BadTolerance => "Specified tolerance is invalid or theoretically unattainable",
        ::Value::Tolerance => "Failed to reach the specified tolerance",
        ::Value::UnderFlow => "Underflow",
        ::Value::OverFlow => "Overflow",
        ::Value::Loss => "Loss of accuracy",
        ::Value::Round => "Roundoff error",
        ::Value::BadLength => "Matrix/vector sizes are not conformant",
        ::Value::NotSquare => "Matrix not square",
        ::Value::Singularity => "Singularity or extremely bad function behavior detected",
        ::Value::Diverge => "Integral or series is divergent",
        ::Value::Unsupported => "The required feature is not supported by this hardware platform",
        ::Value::Unimplemented => "The requested feature is not (yet) implemented",
        ::Value::Cache => "Cache limit exceeded",
        ::Value::Table => "Table limit exceeded",
        ::Value::NoProgress => "Iteration is not making progress towards solution",
        ::Value::NoProgressJacobian => "Jacobian evaluations are not improving the solution",
        ::Value::ToleranceF => "Cannot reach the specified tolerance in F",
        ::Value::ToleranceX => "Cannot reach the specified tolerance in X",
        ::Value::ToleranceG => "Cannot reach the specified tolerance in gradient",
        ::Value::EOF => "End of file",
        ::Value::Unknown(_) => "Unknown error",
    }
}

static mut CALLBACK: Option<fn(&str, &str, u32, ::Value)> = None;

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
/// use rgsl::error::set_error_handler;
/// use rgsl::Value;
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
/// # use rgsl::error::set_error_handler;
/// let old_handler = set_error_handler(None);
/// ```
pub fn set_error_handler(
    f: Option<fn(&str, &str, u32, ::Value)>,
) -> Option<fn(&str, &str, u32, ::Value)> {
    let f = f.into();
    unsafe {
        let out = CALLBACK.take();
        match f {
            Some(f) => {
                CALLBACK = Some(f);
                ffi::gsl_set_error_handler(Some(inner_error_handler));
            }
            None => {
                ffi::gsl_set_error_handler(None);
            }
        }
        out
    }
}

/// This function turns off the error handler by defining an error handler which does nothing. This
/// will cause the program to continue after any error, so the return values from any library
/// routines must be checked. This is the recommended behavior for production programs. The previous
/// handler is returned (so that you can restore it later).
pub fn set_error_handler_off() -> Option<fn(&str, &str, u32, ::Value)> {
    unsafe {
        ffi::gsl_set_error_handler_off();
        CALLBACK.take()
    }
}

extern "C" fn inner_error_handler(
    reason: *const libc::c_char,
    file: *const libc::c_char,
    line: libc::c_int,
    gsl_errno: libc::c_int,
) {
    unsafe { 
        if let Some(ref call) = CALLBACK {
            let s = CStr::from_ptr(reason);
            let f = CStr::from_ptr(file);
            call(s.to_str().unwrap_or_else(|_| "Unknown"),
                 f.to_str().unwrap_or_else(|_| "Unknown"),
                 line as _, ::Value::from(gsl_errno));
        }
    }
}

#[test]
fn test_error_handler() {
    use ::{bessel, Value};

    set_error_handler_off();
    match bessel::K0_e(1e3) {
        (Value::UnderFlow, r) => println!("K0(1e3) underflowed: {:.3e}", r.val),
        _ => panic!("unexpected"),
    }
}
