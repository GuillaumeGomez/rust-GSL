//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#General Polynomial Equations

The roots of polynomial equations cannot be found analytically beyond the special cases of the quadratic, cubic and quartic equation. The algorithm 
described in this section uses an iterative method to find the approximate locations of roots of higher order polynomials.
!*/

use ffi;
use enums;

pub struct PolyComplex {
    w: *mut ffi::gsl_poly_complex_workspace
}

impl PolyComplex {
    /// This function allocates space for a gsl_poly_complex_workspace struct and a workspace suitable for solving a polynomial with n coefficients
    /// using the routine gsl_poly_complex_solve.
    /// 
    /// The function returns a pointer to the newly allocated gsl_poly_complex_workspace if no errors were detected, and a null pointer in the case
    /// of error.
    pub fn new(n: u64) -> Option<PolyComplex> {
        let tmp = unsafe { ffi::gsl_poly_complex_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(PolyComplex {
                w: tmp
            })
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
    pub fn solve(&self, a: &[f64], z: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_poly_complex_solve(a.as_ptr(), a.len() as u64, self.w, z.as_mut_ptr()) }
    }
}

impl Drop for PolyComplex {
    fn drop(&mut self) {
        unsafe { ffi::gsl_poly_complex_workspace_free(self.w) };
        self.w = ::std::ptr::mut_null();
    }
}

impl ffi::FFI<ffi::gsl_poly_complex_workspace> for PolyComplex {
    fn wrap(w: *mut ffi::gsl_poly_complex_workspace) -> PolyComplex {
        PolyComplex {
            w: w
        }
    }

    fn unwrap(w: &PolyComplex) -> *mut ffi::gsl_poly_complex_workspace {
        w.w
    }
}