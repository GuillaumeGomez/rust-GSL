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

pub use types::{
    CblasIndex,
    SF_GAMMA_XMAX,
    SF_FACT_NMAX,
    SF_DOUBLEFACT_NMAX,
    Result,
    ResultE10,
    MatrixF64,
    MatrixF32,
    MatrixComplexF64,
    MatrixComplexF32,
    VectorF64,
    VectorF32,
    VectorComplexF64,
    VectorComplexF32,
    ComplexF64,
    ComplexF32
};

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
pub mod coupling_coefficients;
pub mod dawson;
pub mod debye;
pub mod dilogarithm;
pub mod elementary;
pub mod elementary_operations;
pub mod error;
pub mod exponential;
pub mod exponential_integrals;
pub mod fit;
pub mod gamma_beta;
pub mod jacobian_elliptic;
pub mod pow;
pub mod trigonometric;

#[cfg(target_os = "linux")]
mod platform {
    #[link(name = "gsl")] extern {}
    #[link(name = "gslcblas")] extern {}
}