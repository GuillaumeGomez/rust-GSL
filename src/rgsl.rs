/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

/*!
rust-gsl
========

A __Rust__ binding for the [GSL library](http://www.gnu.org/software/gsl/).

##Installation

This binding requires the __GSL__ library to be installed.

To build it, please use :

```Shell
> make
```

This command build __rgsl__, the examples and the documentation.

You can build them separatly too.

```Shell
> make rgsl
> make examples
> make doc
```

##Documentation

You can access the __rgsl__ documentation locally, just build it :

```Shell
> make doc
```

Then open this file with an internet browser :
file:///{rgsl_location}/doc/rgsl/index.html

## License
__rust-GSL__ is a wrapper for __GSL__, therefore inherits the [GPL licence](http://www.gnu.org/copyleft/gpl.html).

Here is the list of all modules :
!*/

#![crate_name = "rgsl"]
#![desc = "Rust binding for GSL"]
#![crate_type = "rlib"]
#![crate_type = "dylib"]

#![feature(globs)]

#![allow(non_camel_case_types)]
#![allow(non_snake_case_functions)]
#![allow(uppercase_variables)]

extern crate libc;

pub use airy::Airy;
pub use bessel::Bessel;
pub use canonical::Canonical;
pub use cblas::Cblas;

mod ffi;
pub mod enums;
pub mod types;
pub mod airy;
pub mod bessel;
pub mod canonical;
pub mod cblas;

#[cfg(target_os = "linux")]
mod platform {
    #[link(name = "gsl")] extern {}
    #[link(name = "gslcblas")] extern {}
}