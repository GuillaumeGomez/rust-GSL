/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

pub type gsl_mode_t = u32;

pub mod Gsl {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum Mode {
        PrecDouble,
        PrecSingle,
        PrecApprox
    }
}

/// The error handling form of the special functions always calculate an error estimate along with the value of the result. Therefore, structures are provided for amalgamating a value and error estimate.
pub struct GslResult {
    /// Contains the value.
    pub val: f64,
    /// Contains an estimate of the absolute error in the value.
    pub err: f64
}

impl GslResult {
    pub fn new() -> GslResult {
        GslResult {
            val: 0f64,
            err: 0f64
        }
    }
}