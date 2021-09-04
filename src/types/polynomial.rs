//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# General Polynomial Equations

The roots of polynomial equations cannot be found analytically beyond the special cases of the quadratic, cubic and quartic equation. The algorithm
described in this section uses an iterative method to find the approximate locations of roots of higher order polynomials.
!*/

use crate::Value;
use ffi::FFI;

ffi_wrapper!(
    PolyComplexWorkspace,
    *mut sys::gsl_poly_complex_workspace,
    gsl_poly_complex_workspace_free
);

impl PolyComplexWorkspace {
    /// This function allocates space for a gsl_poly_complex_workspace struct and a workspace suitable for solving a polynomial with n coefficients
    /// using the routine gsl_poly_complex_solve.
    ///
    /// The function returns a pointer to the newly allocated gsl_poly_complex_workspace if no errors were detected, and a null pointer in the case
    /// of error.
    #[doc(alias = "gsl_poly_complex_workspace_alloc")]
    pub fn new(n: usize) -> Option<Self> {
        let tmp = unsafe { sys::gsl_poly_complex_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function computes the roots of the general polynomial P(x) = a_0 + a_1 x + a_2 x^2 + ... + a_{n-1} x^{n-1} using balanced-QR reduction
    /// of the companion matrix. The parameter n specifies the length of the coefficient array. The coefficient of the highest order term must be
    /// non-zero. The function requires a workspace w of the appropriate size. The n-1 roots are returned in the packed complex array z of length
    /// 2(n-1), alternating real and imaginary parts.
    ///
    /// The function returns Success if all the roots are found. If the QR reduction does not converge, the error handler is invoked with an error
    /// code of Failed. Note that due to finite precision, roots of higher multiplicity are returned as a cluster of simple roots with reduced
    /// accuracy. The solution of polynomials with higher-order roots requires specialized algorithms that take the multiplicity structure into
    /// account (see e.g. Z. Zeng, Algorithm 835, ACM Transactions on Mathematical Software, Volume 30, Issue 2 (2004), pp 218–236).
    #[doc(alias = "gsl_poly_complex_solve")]
    pub fn solve(&mut self, a: &[f64], z: &mut [f64]) -> Value {
        Value::from(unsafe {
            sys::gsl_poly_complex_solve(
                a.as_ptr(),
                a.len() as _,
                self.unwrap_unique(),
                z.as_mut_ptr(),
            )
        })
    }
}
