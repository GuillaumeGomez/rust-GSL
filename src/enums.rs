/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

pub mod Gsl {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum Mode {
        PrecDouble,
        PrecSingle,
        PrecApprox
    }
}