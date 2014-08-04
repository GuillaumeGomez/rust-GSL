//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::default::Default;

pub type gsl_mode_t = u32;
pub struct CblasIndex(pub u32);

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

impl Default for GslResult {
    fn default() -> GslResult {
        GslResult::new()
    }
}

pub struct GslComplex {
    pub data: [f64, ..2]
}

impl Default for GslComplex {
    fn default() -> GslComplex {
        GslComplex {
            data: [0f64, 0f64]
        }
    }
}

pub struct GslComplexFloat {
    pub data: [f32, ..2]
}

impl Default for GslComplexFloat {
    fn default() -> GslComplexFloat {
        GslComplexFloat {
            data: [0f32, 0f32]
        }
    }
}