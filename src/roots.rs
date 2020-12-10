//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[doc(alias = "gsl_root_test_interval")]
pub fn test_interval(x_lower: f64, x_upper: f64, epsabs: f64, epsrel: f64) -> ::Value {
    ::Value::from(unsafe { sys::gsl_root_test_interval(x_lower, x_upper, epsabs, epsrel) })
}

#[doc(alias = "gsl_root_test_residual")]
pub fn test_residual(f: f64, epsabs: f64) -> ::Value {
    ::Value::from(unsafe { sys::gsl_root_test_residual(f, epsabs) })
}

#[doc(alias = "gsl_root_test_delta")]
pub fn test_delta(x1: f64, x0: f64, epsabs: f64, epsrel: f64) -> ::Value {
    ::Value::from(unsafe { sys::gsl_root_test_delta(x1, x0, epsabs, epsrel) })
}
