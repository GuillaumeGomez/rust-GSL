//
// FFI binding for the GSL library
//

#![allow(improper_ctypes)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

pub extern crate libc;

mod auto;

pub use auto::*;
