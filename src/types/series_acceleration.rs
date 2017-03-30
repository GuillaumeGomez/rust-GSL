//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Series Acceleration

The functions described in this chapter accelerate the convergence of a series using the Levin u-transform. This method takes a small 
number of terms from the start of a series and uses a systematic approximation to compute an extrapolated value and an estimate of its 
error. The u-transform works for both convergent and divergent series, including asymptotic series.

##Acceleration functions

The following functions compute the full Levin u-transform of a series with its error estimate. The error estimate is computed by 
propagating rounding errors from each term through to the final extrapolation.

These functions are intended for summing analytic series where each term is known to high accuracy, and the rounding errors are assumed 
to originate from finite precision. They are taken to be relative errors of order GSL_DBL_EPSILON for each term.

The calculation of the error in the extrapolated value is an O(N^2) process, which is expensive in time and memory. A faster but less 
reliable method which estimates the error from the convergence of the extrapolated value is described in the next section. For the method 
described here a full table of intermediate values and derivatives through to O(N) must be computed and stored, but this does give a 
reliable error estimate.

##Acceleration functions without error estimation

The functions described in this section compute the Levin u-transform of series and attempt to estimate the error from the “truncation 
error” in the extrapolation, the difference between the final two approximations. Using this method avoids the need to compute an 
intermediate table of derivatives because the error is estimated from the behavior of the extrapolated value itself. Consequently this 
algorithm is an O(N) process and only requires O(N) terms of storage. If the series converges sufficiently fast then this procedure can 
be acceptable. It is appropriate to use this method when there is a need to compute many extrapolations of series with similar convergence 
properties at high-speed. For example, when numerically integrating a function defined by a parameterized series where the parameter 
varies only slightly. A reliable error estimate should be computed first using the full algorithm described above in order to verify the 
consistency of the results.

##References and Further Reading

The algorithms used by these functions are described in the following papers,

T. Fessler, W.F. Ford, D.A. Smith, HURRY: An acceleration algorithm for scalar sequences and series ACM Transactions on Mathematical 
Software, 9(3):346–354, 1983. and Algorithm 602 9(3):355–357, 1983.

The theory of the u-transform was presented by Levin,

D. Levin, Development of Non-Linear Transformations for Improving Convergence of Sequences, Intern. J. Computer Math. B3:371–388, 1973.

A review paper on the Levin Transform is available online,

Herbert H. H. Homeier, Scalar Levin-Type Sequence Transformations, http://arxiv.org/abs/math/0005209.
!*/

use ffi;
use enums;

/// Workspace for Levin U Transform with error estimation
pub struct LevinUWorkspace {
    w: *mut ffi::gsl_sum_levin_u_workspace
}

impl LevinUWorkspace {
    /// This function allocates a workspace for a Levin u-transform of n terms. The size of the workspace is O(2n^2 + 3n).
    pub fn new(n: usize) -> Option<LevinUWorkspace> {
        let tmp = unsafe { ffi::gsl_sum_levin_u_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(LevinUWorkspace {
                w: tmp
            })
        }
    }

    /// This function takes the terms of a series in array of size array_size and computes the extrapolated limit of the series using a
    /// Levin u-transform. Additional working space must be provided in w. The extrapolated sum is stored in sum_accel, with an estimate
    /// of the absolute error stored in abserr. The actual term-by-term sum is returned in w->sum_plain. The algorithm calculates the
    /// truncation error (the difference between two successive extrapolations) and round-off error (propagated from the individual terms)
    /// to choose an optimal number of terms for the extrapolation. All the terms of the series passed in through array should be non-zero.
    pub fn accel(&mut self, array: &[f64], sum_accel: &mut f64, abserr: &mut f64) -> enums::Value {
        unsafe { ffi::gsl_sum_levin_u_accel(array.as_ptr(), array.len() as usize, self.w, sum_accel, abserr) }
    }

    pub fn sum_plain(&self) -> f64 {
        unsafe { (*self.w).sum_plain }
    }

    pub fn terms_used(&self) -> usize {
        unsafe { (*self.w).terms_used }
    }

    pub fn size(&self) -> usize {
        unsafe { (*self.w).size }
    }
}

impl Drop for LevinUWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_sum_levin_u_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_sum_levin_u_workspace> for LevinUWorkspace {
    fn wrap(w: *mut ffi::gsl_sum_levin_u_workspace) -> LevinUWorkspace {
        LevinUWorkspace {
            w: w
        }
    }

    fn soft_wrap(w: *mut ffi::gsl_sum_levin_u_workspace) -> LevinUWorkspace {
        Self::wrap(w)
    }

    fn unwrap_shared(w: &LevinUWorkspace) -> *const ffi::gsl_sum_levin_u_workspace {
        w.w as *const _
    }

    fn unwrap_unique(w: &mut LevinUWorkspace) -> *mut ffi::gsl_sum_levin_u_workspace {
        w.w
    }
}

/// The following functions perform the same calculation without estimating the errors. They require O(N) storage instead of O(N^2).
/// This may be useful for summing many similar series where the size of the error has already been estimated reliably and is not
/// expected to change.
pub struct LevinUTruncWorkspace {
    w: *mut ffi::gsl_sum_levin_utrunc_workspace
}

impl LevinUTruncWorkspace {
    /// This function allocates a workspace for a Levin u-transform of n terms, without error estimation. The size of the workspace is O(3n).
    pub fn new(n: usize) -> Option<LevinUTruncWorkspace> {
        let tmp = unsafe { ffi::gsl_sum_levin_utrunc_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(LevinUTruncWorkspace {
                w: tmp
            })
        }
    }

    /// This function takes the terms of a series in array of size array_size and computes the extrapolated limit of the series using a
    /// Levin u-transform. Additional working space must be provided in w. The extrapolated sum is stored in sum_accel. The actual
    /// term-by-term sum is returned in w->sum_plain. The algorithm terminates when the difference between two successive extrapolations
    /// reaches a minimum or is sufficiently small. The difference between these two values is used as estimate of the error and is stored
    /// in abserr_trunc. To improve the reliability of the algorithm the extrapolated values are replaced by moving averages when
    /// calculating the truncation error, smoothing out any fluctuations.
    pub fn accel(&mut self, array: &[f64], sum_accel: &mut f64, abserr_trunc: &mut f64) -> enums::Value {
        unsafe { ffi::gsl_sum_levin_utrunc_accel(array.as_ptr(), array.len() as usize, self.w, sum_accel, abserr_trunc) }
    }

    pub fn sum_plain(&self) -> f64 {
        unsafe { (*self.w).sum_plain }
    }

    pub fn terms_used(&self) -> usize {
        unsafe { (*self.w).terms_used }
    }

    pub fn size(&self) -> usize {
        unsafe { (*self.w).size }
    }
}

impl Drop for LevinUTruncWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_sum_levin_utrunc_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_sum_levin_utrunc_workspace> for LevinUTruncWorkspace {
    fn wrap(w: *mut ffi::gsl_sum_levin_utrunc_workspace) -> LevinUTruncWorkspace {
        LevinUTruncWorkspace {
            w: w
        }
    }

    fn soft_wrap(w: *mut ffi::gsl_sum_levin_utrunc_workspace) -> LevinUTruncWorkspace {
        Self::wrap(w)
    }

    fn unwrap_shared(w: &LevinUTruncWorkspace) -> *const ffi::gsl_sum_levin_utrunc_workspace {
        w.w as *const _
    }

    fn unwrap_unique(w: &mut LevinUTruncWorkspace) -> *mut ffi::gsl_sum_levin_utrunc_workspace {
        w.w
    }
}
