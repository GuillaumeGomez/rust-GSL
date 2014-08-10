//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

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

Since this project supports cargo, you can also build it like this :

```Shell
> cargo build
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
#![feature(macro_rules)]

#![allow(non_camel_case_types)]
#![allow(non_snake_case_functions)]
#![allow(uppercase_variables)]

extern crate libc;

pub use airy::Airy;
pub use bessel::Bessel;
pub use blas::Blas;
pub use canonical::Canonical;
pub use cblas::Cblas;
pub use clausen::Clausen;
pub use coulomb::Coulomb;
pub use elementary::Elementary;
pub use exponential_integrals::ExponentialIntegrals;
pub use fit::Fit;
pub use pow::Pow;
pub use trigonometric::Trigonometric;
pub use types::Gsl;

pub use self::enums::{
    mode,
    cblas_order,
    cblas_side,
    cblas_transpose,
    cblas_uplo,
    cblas_diag,
    gsl_value
};

pub use self::enums::{
    Mode,
    CblasOrder,
    CblasSide,
    CblasTranspose,
    CblasUplo,
    CblasDiag,
    GslValue
};

mod ffi;
pub mod enums;
pub mod types;
pub mod airy;
pub mod bessel;
pub mod blas;
pub mod canonical;
pub mod cblas;
pub mod clausen;
pub mod coulomb;
pub mod elementary;
pub mod exponential_integrals;
pub mod fit;
pub mod pow;
pub mod trigonometric;

#[cfg(target_os = "linux")]
mod platform {
    #[link(name = "gsl")] extern {}
    #[link(name = "gslcblas")] extern {}
}