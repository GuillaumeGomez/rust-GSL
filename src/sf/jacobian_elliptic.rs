//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Jacobian Elliptic functions are defined in Abramowitz & Stegun, Chapter 16.
!*/

use crate::Error;

/// This function computes the Jacobian elliptic functions sn(u|m), cn(u|m), dn(u|m) by descending
/// Landen transformations.
///
/// Returns `(sn, cn, dn)`.
#[doc(alias = "gsl_sf_elljac_e")]
pub fn elljac_e(u: f64, m: f64) -> Result<(f64, f64, f64), Error> {
    let mut sn = 0.;
    let mut cn = 0.;
    let mut dn = 0.;
    let ret = unsafe { sys::gsl_sf_elljac_e(u, m, &mut sn, &mut cn, &mut dn) };
    Error::handle(ret, (sn, cn, dn))
}
