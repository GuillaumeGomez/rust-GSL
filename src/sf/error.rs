//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The error function is described in Abramowitz & Stegun, Chapter 7.

use crate::{types, Error};
use std::mem::MaybeUninit;

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
#[doc(alias = "gsl_sf_erf")]
pub fn erf(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erf(x) }
}

/// This routine computes the error function erf(x), where erf(x) = (2/\sqrt(\pi)) \int_0^x dt \exp(-t^2).
#[doc(alias = "gsl_sf_erf_e")]
pub fn erf_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erf_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
#[doc(alias = "gsl_sf_erfc")]
pub fn erfc(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erfc(x) }
}

/// This routine computes the complementary error function erfc(x) = 1 - erf(x) = (2/\sqrt(\pi)) \int_x^\infty \exp(-t^2).
#[doc(alias = "gsl_sf_erfc_e")]
pub fn erfc_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erfc_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
#[doc(alias = "gsl_sf_log_erfc")]
pub fn log_erfc(x: f64) -> f64 {
    unsafe { sys::gsl_sf_log_erfc(x) }
}

/// This routine computes the logarithm of the complementary error function \log(\erfc(x)).
#[doc(alias = "gsl_sf_log_erfc_e")]
pub fn log_erfc_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_log_erfc_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
#[doc(alias = "gsl_sf_erf_Z")]
pub fn erf_Z(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erf_Z(x) }
}

/// This routine computes the Gaussian probability density function Z(x) = (1/\sqrt{2\pi}) \exp(-x^2/2).
#[doc(alias = "gsl_sf_erf_Z_e")]
pub fn erf_Z_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erf_Z_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the upper tail of the Gaussian probability function Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2).
///
/// The hazard function for the normal distribution, also known as the inverse Mills’ ratio, is defined as,
///
/// h(x) = Z(x)/Q(x) = \sqrt{2/\pi} \exp(-x^2 / 2) / \erfc(x/\sqrt 2)
///
/// It decreases rapidly as x approaches -\infty and asymptotes to h(x) \sim x as x approaches +\infty.
#[doc(alias = "gsl_sf_erf_Q")]
pub fn erf_Q(x: f64) -> f64 {
    unsafe { sys::gsl_sf_erf_Q(x) }
}

/// This routine computes the upper tail of the Gaussian probability function Q(x) = (1/\sqrt{2\pi}) \int_x^\infty dt \exp(-t^2/2).
///
/// The hazard function for the normal distribution, also known as the inverse Mills’ ratio, is defined as,
///
/// h(x) = Z(x)/Q(x) = \sqrt{2/\pi} \exp(-x^2 / 2) / \erfc(x/\sqrt 2)
///
/// It decreases rapidly as x approaches -\infty and asymptotes to h(x) \sim x as x approaches +\infty.
#[doc(alias = "gsl_sf_erf_Q_e")]
pub fn erf_Q_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_erf_Q_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}

/// This routine computes the hazard function for the normal distribution.
#[doc(alias = "gsl_sf_hazard")]
pub fn hazard(x: f64) -> f64 {
    unsafe { sys::gsl_sf_hazard(x) }
}

/// This routine computes the hazard function for the normal distribution.
#[doc(alias = "gsl_sf_hazard_e")]
pub fn hazard_e(x: f64) -> Result<types::Result, Error> {
    let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
    let ret = unsafe { sys::gsl_sf_hazard_e(x, result.as_mut_ptr()) };

    Error::handle(ret, unsafe { result.assume_init() }.into())
}
