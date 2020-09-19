//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::default::Default;

/// The error handling form of the special functions always calculate an error estimate along with the value of the result.
/// Therefore, structures are provided for amalgamating a value and error estimate.
#[derive(Clone, Copy)]
pub struct Result {
    /// Contains the value.
    pub val: f64,
    /// Contains an estimate of the absolute error in the value.
    pub err: f64,
}

impl Default for Result {
    fn default() -> Result {
        Result::new()
    }
}

impl Result {
    pub fn new() -> Result {
        Result {
            val: 0f64,
            err: 0f64,
        }
    }
}

impl From<::sys::gsl_sf_result> for Result {
    fn from(v: ::sys::gsl_sf_result) -> Self {
        Self {
            val: v.val,
            err: v.err,
        }
    }
}

/// In some cases, an overflow or underflow can be detected and handled by a function.
/// In this case, it may be possible to return a scaling exponent as well as an error/value pair in order to save the result from exceeding the dynamic range of the built-in types.
#[derive(Clone, Copy)]
pub struct ResultE10 {
    /// Contains the value.
    pub val: f64,
    /// Contains an estimate of the absolute error in the value.
    pub err: f64,
    /// Exponent field such that the actual result is obtained as result * 10^(e10).
    pub e10: i32,
}

impl Default for ResultE10 {
    fn default() -> ResultE10 {
        ResultE10::new()
    }
}

impl ResultE10 {
    pub fn new() -> ResultE10 {
        ResultE10 {
            val: 0f64,
            err: 0f64,
            e10: 0i32,
        }
    }
}

impl From<::sys::gsl_sf_result_e10> for ResultE10 {
    fn from(v: ::sys::gsl_sf_result_e10) -> Self {
        Self {
            val: v.val,
            err: v.err,
            e10: v.e10,
        }
    }
}
