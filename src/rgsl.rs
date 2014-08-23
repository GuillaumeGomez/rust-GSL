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
#![allow(ctypes)]

extern crate libc;

pub use types::{
    ComplexF32,
    ComplexF64,
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
    Rng,
    RngType,
    Permutation,
    ChebSeries,
    Combination,
    PolyComplex,
    DiscreteHankel,
    EigenSymmetricWorkspace,
    EigenSymmetricVWorkspace,
    EigenHermitianWorkspace,
    EigenHermitianVWorkspace,
    EigenNonSymmWorkspace,
    EigenNonSymmVWorkspace,
    EigenGenSymmWorkspace,
    EigenGenSymmVWorkspace,
    EigenGenHermWorkspace,
    EigenGenHermVWorkspace,
    EigenGenWorkspace,
    EigenGenVWorkspace
};

pub use enums::{
    Mode,
    Value
};

pub use elementary::Elementary;
pub use pow::Pow;
pub use trigonometric::Trigonometric;
pub use types::rng;

mod ffi;
pub mod randist;
pub mod types;
pub mod enums;

pub mod airy;
pub mod bessel;
pub mod blas;
pub mod cblas;
pub mod clausen;
pub mod coulomb;
pub mod coupling_coefficients;
pub mod dawson;
pub mod debye;
pub mod dilogarithm;
pub mod elementary;
pub mod elementary_operations;
pub mod elliptic;
pub mod error;
pub mod exponential;
pub mod exponential_integrals;
pub mod fermi_dirac;
pub mod fit;
pub mod gamma_beta;
pub mod gegenbauer;
pub mod hypergeometric;
pub mod jacobian_elliptic;
pub mod laguerre;
pub mod lambert_w;
pub mod legendre;
pub mod logarithm;
pub mod numerical_differentiation;
pub mod polynomials;
pub mod pow;
pub mod power;
pub mod psi;
pub mod sort;
pub mod synchrotron;
pub mod transport;
pub mod trigonometric;
pub mod zeta;

pub type comparison_fn<T> = fn(a: &T, b: &T) -> i32;
pub type function<T> = fn(x: f64, p: &mut T) -> f64;
//pub type ComplexPackedPtr = &mut [f64];

/// The maximum x such that gamma(x) is not considered an overflow.
pub static SF_GAMMA_XMAX : f64 = 171.0;
/// The maximum n such that gsl_sf_fact(n) does not give an overflow.
pub static SF_FACT_NMAX : f64 = 170.0;
/// The maximum n such that gsl_sf_doublefact(n) does not give an overflow.
pub static SF_DOUBLEFACT_NMAX : f64 = 297.0;

pub static SF_MATHIEU_COEFF : u32 = 100;

pub static DBL_EPSILON       : f64 = 2.2204460492503131e-16;
pub static SQRT_DBL_EPSILON  : f64 = 1.4901161193847656e-08;
pub static ROOT3_DBL_EPSILON : f64 = 6.0554544523933429e-06;
pub static ROOT4_DBL_EPSILON : f64 = 1.2207031250000000e-04;
pub static ROOT5_DBL_EPSILON : f64 = 7.4009597974140505e-04;
pub static ROOT6_DBL_EPSILON : f64 = 2.4607833005759251e-03;

pub static DBL_MIN           : f64 = 2.2250738585072014e-308;
pub static SQRT_DBL_MIN      : f64 = 1.4916681462400413e-154;
pub static ROOT3_DBL_MIN     : f64 = 2.8126442852362996e-103;
pub static ROOT4_DBL_MIN     : f64 = 1.2213386697554620e-77;
pub static ROOT5_DBL_MIN     : f64 = 2.9476022969691763e-62;
pub static ROOT6_DBL_MIN     : f64 = 5.3034368905798218e-52;

//pub static DBL_MAX           : f64 = 1.7976931348623157e+308;
pub static SQRT_DBL_MAX      : f64 = 1.3407807929942596e+154;
pub static ROOT3_DBL_MAX     : f64 = 5.6438030941222897e+102;
pub static ROOT4_DBL_MAX     : f64 = 1.1579208923731620e+77;
pub static ROOT5_DBL_MAX     : f64 = 4.4765466227572707e+61;
pub static ROOT6_DBL_MAX     : f64 = 2.3756689782295612e+51;
pub static LOG_DBL_MAX       : f64 = 7.0978271289338397e+02;

#[cfg(target_os = "linux")]
mod platform {
    #[link(name = "gsl")] extern {}
    #[link(name = "gslcblas")] extern {}
}