//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
B-splines are commonly used as basis functions to fit smoothing curves to large data sets.
To do this, the abscissa axis is broken up into some number of intervals, where the endpoints of
each interval are called breakpoints.

These breakpoints are then converted to knots by imposing various continuity and smoothness
conditions at each interface. Given a nondecreasing knot vector t = {t_0, t_1, …, t_{n+k-1}}, the n
basis splines of order k are defined by

```latex
B_(i,1)(x) = (1, t_i <= x < t_(i+1)

             (0, else

B_(i,k)(x) = [(x - t_i)/(t_(i+k-1) - t_i)] B_(i,k-1)(x)

              + [(t_(i+k) - x)/(t_(i+k) - t_(i+1))] B_(i+1,k-1)(x)
```

for i = 0, …, n-1. The common case of cubic B-splines is given by k = 4. The above recurrence
relation can be evaluated in a numerically stable way by the de Boor algorithm.

If we define appropriate knots on an interval [a,b] then the B-spline basis functions form a
complete set on that interval. Therefore we can expand a smoothing function as

f(x) = \sum_i c_i B_(i,k)(x)

given enough (x_j, f(x_j)) data pairs. The coefficients c_i can be readily obtained from a
least-squares fit.

###References and Further Reading

Further information on the algorithms described in this section can be found in the following book,

C. de Boor, A Practical Guide to Splines (1978), Springer-Verlag, ISBN 0-387-90356-9.
Further information of Greville abscissae and B-spline collocation can be found in the following
paper,

Richard W. Johnson, Higher order B-spline collocation at the Greville abscissae. Applied Numerical
Mathematics. vol. 52, 2005, 63–75.

A large collection of B-spline routines is available in the PPPACK library available at
http://www.netlib.org/pppack, which is also part of SLATEC.
!*/

#[cfg(not(feature = "v2"))]
use types::{VectorF64, MatrixF64};
#[cfg(feature = "v2")]
use types::VectorF64;
use ffi;
use enums;

pub struct BSpLineWorkspace {
    w: *mut ffi::gsl_bspline_workspace
}

impl BSpLineWorkspace {
    /// This function allocates a workspace for computing B-splines of order k.
    ///
    /// The number of breakpoints is given by nbreak. This leads to n = nbreak + k - 2 basis
    /// functions.
    ///
    /// Cubic B-splines are specified by k = 4. The size of the workspace is O(5k + nbreak).
    pub fn new(k: usize, nbreak: usize) -> Option<BSpLineWorkspace> {
        let tmp = unsafe { ffi::gsl_bspline_alloc(k, nbreak) };

        if tmp.is_null() {
            None
        } else {
            Some(BSpLineWorkspace {
                w: tmp
            })
        }
    }

    /// This function computes the knots associated with the given breakpoints and stores them
    /// internally in w->knots.
    pub fn knots(&self, breakpts: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_bspline_knots(ffi::FFI::unwrap(breakpts), self.w) }
    }

    /// This function assumes uniformly spaced breakpoints on [a,b] and constructs the corresponding
    /// knot vector using the previously specified nbreak parameter.
    /// The knots are stored in w->knots.
    pub fn knots_uniform(&self, a: f64, b: f64) -> enums::Value {
        unsafe { ffi::gsl_bspline_knots_uniform(a, b, self.w) }
    }

    /// This function evaluates all B-spline basis functions at the position x and stores them in
    /// the vector B, so that the i-th element is B_i(x).
    ///
    /// The vector B must be of length n = nbreak + k - 2. This value may also be obtained by
    /// calling gsl_bspline_ncoeffs.
    ///
    /// Computing all the basis functions at once is more efficient than computing them
    /// individually, due to the nature of the defining recurrence relation.
    pub fn eval(&self, x: f64, B: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_bspline_eval(x, ffi::FFI::unwrap(B), self.w) }
    }

    /// This function evaluates all potentially nonzero B-spline basis functions at the position x
    /// and stores them in the vector Bk, so that the i-th element is B_(istart+i)(x).
    ///
    /// The last element of Bk is B_(iend)(x). The vector Bk must be of length k.
    /// By returning only the nonzero basis functions, this function allows quantities involving
    /// linear combinations of the B_i(x) to be computed without unnecessary terms (such linear
    /// combinations occur, for example, when evaluating an interpolated function).
    pub fn eval_non_zero(&self, x: f64, Bk: &VectorF64, istart: &mut usize,
                         iend: &mut usize) -> enums::Value {
        unsafe { ffi::gsl_bspline_eval_nonzero(x, ffi::FFI::unwrap(Bk), istart, iend, self.w) }
    }

    /// This function returns the number of B-spline coefficients given by n = nbreak + k - 2.
    pub fn ncoeffs(&self) -> usize {
        unsafe { ffi::gsl_bspline_ncoeffs(self.w) }
    }

    /// The Greville abscissae are defined to be the mean location of k-1 consecutive knots in the
    /// knot vector for each basis spline function of order k.
    ///
    /// With the first and last knots in the gsl_bspline_workspace knot vector excluded, there are
    /// gsl_bspline_ncoeffs Greville abscissae for any given B-spline basis.
    ///
    /// These values are often used in B-spline collocation applications and may also be called
    /// Marsden-Schoenberg points.
    ///
    /// Returns the location of the i-th Greville abscissa for the given B-spline basis.
    /// For the ill-defined case when k=1, the implementation chooses to return breakpoint interval
    /// midpoints.
    pub fn greville_abscissa(&self, i: usize) -> f64 {
        unsafe { ffi::gsl_bspline_greville_abscissa(i, self.w) }
    }
}

impl Drop for BSpLineWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_bspline_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_bspline_workspace> for BSpLineWorkspace {
    fn wrap(r: *mut ffi::gsl_bspline_workspace) -> BSpLineWorkspace {
        BSpLineWorkspace {
            w: r
        }
    }

    fn unwrap(bsp: &BSpLineWorkspace) -> *mut ffi::gsl_bspline_workspace {
        bsp.w
    }
}

#[cfg(not(feature = "v2"))]
pub struct BSpLineDerivWorkspace {
    w: *mut ffi::gsl_bspline_deriv_workspace
}

#[cfg(not(feature = "v2"))]
impl BSpLineDerivWorkspace {
    /// This function allocates a workspace for computing the derivatives of a B-spline basis
    /// function of order k.
    ///
    /// The size of the workspace is O(2k^2).
    pub fn new(k: usize) -> Option<BSpLineDerivWorkspace> {
        let tmp = unsafe { ffi::gsl_bspline_deriv_alloc(k) };

        if tmp.is_null() {
            None
        } else {
            Some(BSpLineDerivWorkspace {
                w: tmp
            })
        }
    }

    /// This function evaluates all B-spline basis function derivatives of orders 0 through nderiv
    /// (inclusive) at the position x and stores them in the matrix dB.
    /// The (i,j)-th element of dB is d^jB_i(x)/dx^j. The matrix dB must be of size n = nbreak + k -
    /// 2 by nderiv + 1.
    /// The value n may also be obtained by calling gsl_bspline_ncoeffs. Note that function
    /// evaluations are included as the zeroth order derivatives in dB.
    /// Computing all the basis function derivatives at once is more efficient than computing them
    /// individually, due to the nature of the defining recurrence relation.
    pub fn eval(&self, x: f64, nderiv: usize, dB: &MatrixF64,
                w: &BSpLineWorkspace) -> enums::Value {
        unsafe { ffi::gsl_bspline_deriv_eval(x, nderiv, ffi::FFI::unwrap(dB), ffi::FFI::unwrap(w),
                                             self.w) }
    }

    /// This function evaluates all potentially nonzero B-spline basis function derivatives of
    /// orders 0 through nderiv (inclusive) at the position x and stores them in the matrix dB.
    ///
    /// The (i,j)-th element of dB is d^j/dx^j B_(istart+i)(x). The last row of dB contains d^j/dx^j
    /// B_(iend)(x).
    ///
    /// The matrix dB must be of size k by at least nderiv + 1. Note that function evaluations are
    /// included as the zeroth order derivatives in dB.
    ///
    /// By returning only the nonzero basis functions, this function allows quantities involving
    /// linear combinations of the B_i(x) and their derivatives to be computed without unnecessary
    /// terms.
    pub fn eval_non_zero(&self, x: f64, nderiv: usize, dB: &MatrixF64, istart: &mut usize,
                         iend: &mut usize, w: &BSpLineWorkspace) -> enums::Value {
        unsafe { ffi::gsl_bspline_deriv_eval_nonzero(x, nderiv, ffi::FFI::unwrap(dB), istart, iend,
                                                     ffi::FFI::unwrap(w), self.w) }
    }
}

#[cfg(not(feature = "v2"))]
impl Drop for BSpLineDerivWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_bspline_deriv_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

#[cfg(not(feature = "v2"))]
impl ffi::FFI<ffi::gsl_bspline_deriv_workspace> for BSpLineDerivWorkspace {
    fn wrap(r: *mut ffi::gsl_bspline_deriv_workspace) -> BSpLineDerivWorkspace {
        BSpLineDerivWorkspace {
            w: r
        }
    }

    fn unwrap(bsp: &BSpLineDerivWorkspace) -> *mut ffi::gsl_bspline_deriv_workspace {
        bsp.w
    }
}
