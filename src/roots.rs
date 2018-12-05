//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;

pub fn test_interval(x_lower: f64, x_upper: f64, epsabs: f64, epsrel: f64) -> ::Value {
    ::Value::from(unsafe {
        ffi::gsl_root_test_interval(x_lower, x_upper, epsabs, epsrel)
    })
}

pub fn test_residual(f: f64, epsabs: f64) -> ::Value {
    ::Value::from(unsafe {
        ffi::gsl_root_test_residual(f, epsabs)
    })
}

pub fn test_delta(x1: f64, x0: f64, epsabs: f64, epsrel: f64) -> ::Value {
    ::Value::from(unsafe {
        ffi::gsl_root_test_delta(x1, x0, epsabs, epsrel)
    })
}
