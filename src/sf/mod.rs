//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Special functions.
//!
//! # Usage
//!
//! The special functions are available in two calling conventions, a
//! _natural_ form which returns the numerical value of the function and
//! an _error-handling_ form which returns an error code. The two types
//! of function provide alternative ways of accessing the same
//! underlying code.

pub mod airy;
pub mod bessel;
pub mod clausen;
pub mod coulomb;
pub mod coupling_coefficients;
pub mod dawson;
pub mod debye;
pub mod dilogarithm;
pub mod elementary_operations;
pub mod elliptic;
pub mod error;
pub mod exponential;
pub mod exponential_integrals;
pub mod fermi_dirac;
pub mod gamma_beta;
pub mod gegenbauer;
pub mod hypergeometric;
pub mod jacobian_elliptic;
pub mod laguerre;
pub mod lambert_w;
pub mod legendre;
pub mod logarithm;
pub mod power;
pub mod psi;
pub mod synchrotron;
pub mod transport;
pub mod trigonometric;
pub mod zeta;
