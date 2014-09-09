//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Interpolation

This chapter describes functions for performing interpolation. The library provides a variety of interpolation methods, including Cubic splines 
and Akima splines. The interpolation types are interchangeable, allowing different methods to be used without recompiling. Interpolations can 
be defined for both normal and periodic boundary conditions. Additional functions are available for computing derivatives and integrals of 
interpolating functions.

These interpolation methods produce curves that pass through each datapoint. To interpolate noisy data with a smoothing curve see Basis Splines.

##Introduction

Given a set of data points (x_1, y_1) \dots (x_n, y_n) the routines described in this section compute a continuous interpolating function 
y(x) such that y(x_i) = y_i. The interpolation is piecewise smooth, and its behavior at the end-points is determined by the type of 
interpolation used.

##Index Look-up and Acceleration

The state of searches can be stored in a gsl_interp_accel object, which is a kind of iterator for interpolation lookups. It caches the previous 
value of an index lookup. When the subsequent interpolation point falls in the same interval its index value can be returned immediately.
!*/

use ffi;
use enums;

/// evaluation accelerator
#[repr(C)]
pub struct InterpAccel {
    /// cache of index
    pub cache: u64,
    /// keep statistics
    pub miss_count: u64,
    pub hit_count: u64
}

impl InterpAccel {
    /// This function returns the index i of the array x_array such that x_array[i] <= x < x_array[i+1]. The index is searched for in the
    /// range [index_lo,index_hi].
    pub fn bsearch(x_array: &[f64], x: f64, index_lo: u64, index_hi: u64) -> u64 {
        unsafe { ffi::gsl_interp_bsearch(x_array.as_ptr(), x, index_lo, index_hi) }
    }

    /// This function returns a pointer to an accelerator object, which is a kind of iterator for interpolation lookups. It tracks the state
    /// of lookups, thus allowing for application of various acceleration strategies.
    pub fn new() -> InterpAccel {
        InterpAccel {
            cache: 0u64,
            miss_count: 0u64,
            hit_count: 0u64
        }
    }

    /// This function reinitializes the accelerator object acc. It should be used when the cached information is no longer applicable—for
    /// example, when switching to a new dataset.
    pub fn reset(&mut self) {
        self.cache = 0u64;
        self.miss_count = 0u64;
        self.hit_count = 0u64;
    }

    /// This function performs a lookup action on the data array x_array of size size, using the given accelerator a. This is how lookups
    /// are performed during evaluation of an interpolation. The function returns an index i such that x_array[i] <= x < x_array[i+1].
    pub fn find(&mut self, x_array: &[f64], x: f64) -> u64 {
        unsafe { ffi::gsl_interp_accel_find(self, x_array.as_ptr(), x_array.len() as u64, x) }
    }
}

pub struct Interp {
    interp: *mut ffi::gsl_interp
}

impl Interp {
    /// This function returns a pointer to a newly allocated interpolation object of type T for size data-points.
    pub fn new(t: &InterpType, size: u64) -> Option<Interp> {
        let tmp = unsafe { ffi::gsl_interp_alloc(t.t, size) };

        if tmp.is_null() {
            None
        } else {
            Some(Interp {
                interp: tmp
            })
        }
    }

    /// This function initializes the interpolation object interp for the data (xa,ya) where xa and ya are arrays of size size. The interpolation
    /// object (gsl_interp) does not save the data arrays xa and ya and only stores the static state computed from the data. The xa data array
    /// is always assumed to be strictly ordered, with increasing x values; the behavior for other arrangements is not defined.
    pub fn init(&self, xa: &[f64], ya: &[f64]) -> enums::Value {
        unsafe { ffi::gsl_interp_init(self.interp, xa.as_ptr(), ya.as_ptr(), xa.len() as u64) }
    }

    /// This function returns the name of the interpolation type used by interp. For example,
    /// 
    /// ```Rust
    /// println!("interp uses '{}' interpolation.", interp.name());
    /// ```
    /// 
    /// would print something like :
    ///
    /// ```Shell
    /// interp uses 'cspline' interpolation.
    /// ```
    pub fn name(&self) -> String {
        let tmp = unsafe { ffi::gsl_interp_name(self.interp as *const ffi::gsl_interp) };

        if tmp.is_null() {
            String::new()
        } else {
            unsafe { ::std::string::raw::from_buf(tmp as *const u8) }
        }
    }

    /// This function returns the minimum number of points required by the interpolation object interp or interpolation type T. For example,
    /// Akima spline interpolation requires a minimum of 5 points.
    pub fn min_size(&self) -> u32 {
        unsafe { ffi::gsl_interp_min_size(self.interp as *const ffi::gsl_interp) }
    }
}

impl Drop for Interp {
    fn drop(&mut self) {
        unsafe { ffi::gsl_interp_free(self.interp) };
        self.interp = ::std::ptr::mut_null();
    }
}

impl ffi::FFI<ffi::gsl_interp> for Interp {
    fn wrap(interp: *mut ffi::gsl_interp) -> Interp {
        Interp {
            interp: interp
        }
    }

    fn unwrap(interp: &Interp) -> *mut ffi::gsl_interp {
        interp.interp
    }
}

pub struct InterpType {
    t: *const ffi::gsl_interp_type
}

impl InterpType {
    /// This function returns the minimum number of points required by the interpolation object interp or interpolation type T. For example,
    /// Akima spline interpolation requires a minimum of 5 points.
    pub fn min_size(&self) -> u32 {
        unsafe { ffi::gsl_interp_type_min_size(self.t) }
    }

    /// Linear interpolation. This interpolation method does not require any additional memory.
    pub fn linear() -> InterpType {
        ffi::FFI::wrap(ffi::gsl_interp_linear as *mut ffi::gsl_interp_type)
    }

    /// Polynomial interpolation. This method should only be used for interpolating small numbers of points because polynomial interpolation
    /// introduces large oscillations, even for well-behaved datasets. The number of terms in the interpolating polynomial is equal to the
    /// number of points.
    pub fn polynomial() -> InterpType {
        ffi::FFI::wrap(ffi::gsl_interp_polynomial as *mut ffi::gsl_interp_type)
    }

    /// Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second
    /// derivatives at the supplied data-points. The second derivative is chosen to be zero at the first point and last point.
    pub fn cspline() -> InterpType {
        ffi::FFI::wrap(ffi::gsl_interp_cspline as *mut ffi::gsl_interp_type)
    }

    /// Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second
    /// derivatives at the supplied data-points. The derivatives at the first and last points are also matched. Note that the last point in
    /// the data must have the same y-value as the first point, otherwise the resulting periodic interpolation will have a discontinuity at
    /// the boundary.
    pub fn cspline_periodic() -> InterpType {
        ffi::FFI::wrap(ffi::gsl_interp_cspline_periodic as *mut ffi::gsl_interp_type)
    }

    /// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.
    pub fn akima() -> InterpType {
        ffi::FFI::wrap(ffi::gsl_interp_akima as *mut ffi::gsl_interp_type)
    }

    /// Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.
    pub fn akima_periodic() -> InterpType {
        ffi::FFI::wrap(ffi::gsl_interp_akima_periodic as *mut ffi::gsl_interp_type)
    }
}

impl ffi::FFI<ffi::gsl_interp_type> for InterpType {
    fn wrap(t: *mut ffi::gsl_interp_type) -> InterpType {
        InterpType {
            t: t as *const ffi::gsl_interp_type
        }
    }

    fn unwrap(t: &InterpType) -> *mut ffi::gsl_interp_type {
        t.t as *mut ffi::gsl_interp_type
    }
}