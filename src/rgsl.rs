/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

#![crate_name = "rgsl"]
#![desc = "Rust binding for GSL"]
#![crate_type = "rlib"]
#![crate_type = "dylib"]

#![feature(globs)]

#![allow(non_camel_case_types)]
#![allow(non_snake_case_functions)]

extern crate libc;

pub use airy::Airy;
pub use bessel::Bessel;
pub use canonical::Canonical;

mod ffi;
pub mod types;
pub mod airy;
pub mod bessel;
pub mod canonical;

#[cfg(target_os = "linux")]
mod platform {
    #[link(name = "gsl")] extern {}
    #[link(name = "gslcblas")] extern {}
}