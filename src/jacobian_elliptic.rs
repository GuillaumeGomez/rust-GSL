//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Jacobian Elliptic functions are defined in Abramowitz & Stegun, Chapter 16.
!*/

use enums;

/// This function computes the Jacobian elliptic functions sn(u|m), cn(u|m), dn(u|m) by descending Landen transformations.
pub fn elljac_e(u: f64, m: f64, sn: &mut f64, cn: &mut f64, dn: &mut f64) -> enums::GslValue {
    unsafe { ::ffi::gsl_sf_elljac_e(u, m, sn, cn, dn) }
}