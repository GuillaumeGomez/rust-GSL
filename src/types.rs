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

pub struct GslResult {
    pub val: f64,
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