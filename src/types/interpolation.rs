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

##Higher-level Interface

The functions described in the previous sections required the user to supply pointers to the x and y arrays on each call. The following
functions are equivalent to the corresponding gsl_interp functions but maintain a copy of this data in the gsl_spline object. This removes
the need to pass both xa and ya as arguments on each evaluation.

##References and Further Reading

Descriptions of the interpolation algorithms and further references can be found in the following books:

C.W. Ueberhuber, Numerical Computation (Volume 1), Chapter 9 “Interpolation”, Springer (1997), ISBN 3-540-62058-3.
D.M. Young, R.T. Gregory A Survey of Numerical Mathematics (Volume 1), Chapter 6.8, Dover (1988), ISBN 0-486-65691-8.
!*/

use enums;
use ffi;

/// Evaluation accelerator.
#[derive(Clone)]
pub struct InterpAccel(pub sys::gsl_interp_accel);

impl InterpAccel {
    /// This function returns a pointer to an accelerator object, which is a kind of iterator for
    /// interpolation lookups. It tracks the state of lookups, thus allowing for application of
    /// various acceleration strategies.
    pub fn new() -> InterpAccel {
        InterpAccel(sys::gsl_interp_accel {
            cache: 0,
            miss_count: 0,
            hit_count: 0,
        })
    }

    /// This function reinitializes the accelerator object acc. It should be used when the cached
    /// information is no longer applicable-for example, when switching to a new dataset.
    pub fn reset(&mut self) {
        self.0.cache = 0;
        self.0.miss_count = 0;
        self.0.hit_count = 0;
    }

    /// This function performs a lookup action on the data array x_array of size size, using the
    /// given accelerator a. This is how lookups are performed during evaluation of an
    /// interpolation. The function returns an index i such that `x_array[i] <= x < x_array[i+1]`.
    pub fn find(&mut self, x_array: &[f64], x: f64) -> usize {
        unsafe {
            sys::gsl_interp_accel_find(&mut self.0, x_array.as_ptr(), x_array.len() as usize, x)
        }
    }
}

pub struct Interp {
    interp: *mut sys::gsl_interp,
}

impl Interp {
    /// This function returns a pointer to a newly allocated interpolation object of type T for
    /// size data-points.
    ///
    /// ```
    /// use rgsl::{Interp, InterpType};
    ///
    /// let interp_type = InterpType::linear();
    /// let interp = Interp::new(&interp_type, 2).expect("Failed to initialize `Interp`...");
    /// ```
    pub fn new(t: &InterpType, size: usize) -> Option<Interp> {
        let tmp = unsafe { sys::gsl_interp_alloc(t.t, size) };

        if tmp.is_null() {
            None
        } else {
            Some(Interp { interp: tmp })
        }
    }

    /// This function initializes the interpolation object interp for the data (xa,ya) where xa and
    /// ya are arrays of size size. The interpolation object (gsl_interp) does not save the data
    /// arrays xa and ya and only stores the static state computed from the data. The xa data array
    /// is always assumed to be strictly ordered, with increasing x values; the behavior for other
    /// arrangements is not defined.
    pub fn init(&self, xa: &[f64], ya: &[f64]) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_interp_init(self.interp, xa.as_ptr(), ya.as_ptr(), xa.len() as usize)
        })
    }

    /// This function returns the name of the interpolation type used by interp. For example,
    ///
    /// ```
    /// use rgsl::{Interp, InterpType};
    ///
    /// let interp_type = InterpType::linear();
    /// let interp = Interp::new(&interp_type, 2).expect("Failed to initialize `Interp`...");
    /// println!("interp uses '{}' interpolation.", interp.name());
    /// ```
    ///
    /// would print something like :
    ///
    /// ```Shell
    /// interp uses 'cspline' interpolation.
    /// ```
    pub fn name(&self) -> String {
        let tmp = unsafe { sys::gsl_interp_name(self.interp) };

        if tmp.is_null() {
            String::new()
        } else {
            unsafe {
                String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
            }
        }
    }

    /// This function returns the minimum number of points required by the interpolation object
    /// interp or interpolation type T. For example, Akima spline interpolation requires a minimum
    /// of 5 points.
    pub fn min_size(&self) -> u32 {
        unsafe { sys::gsl_interp_min_size(self.interp) }
    }
}

impl Drop for Interp {
    fn drop(&mut self) {
        unsafe { sys::gsl_interp_free(self.interp) };
        self.interp = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<sys::gsl_interp> for Interp {
    fn wrap(interp: *mut sys::gsl_interp) -> Interp {
        Interp { interp: interp }
    }

    fn soft_wrap(interp: *mut sys::gsl_interp) -> Interp {
        Self::wrap(interp)
    }

    fn unwrap_shared(interp: &Interp) -> *const sys::gsl_interp {
        interp.interp as *const _
    }

    fn unwrap_unique(interp: &mut Interp) -> *mut sys::gsl_interp {
        interp.interp
    }
}

#[derive(Clone, Copy)]
pub struct InterpType {
    t: *const sys::gsl_interp_type,
}

impl InterpType {
    /// This function returns the minimum number of points required by the interpolation object
    /// interp or interpolation type T. For example, Akima spline interpolation requires a minimum
    /// of 5 points.
    pub fn min_size(&self) -> u32 {
        unsafe { sys::gsl_interp_type_min_size(self.t) }
    }

    /// Linear interpolation. This interpolation method does not require any additional memory.
    pub fn linear() -> InterpType {
        ffi_wrap!(gsl_interp_linear, gsl_interp_type)
    }

    /// Polynomial interpolation. This method should only be used for interpolating small numbers
    /// of points because polynomial interpolation introduces large oscillations, even for
    /// well-behaved datasets. The number of terms in the interpolating polynomial is equal to the
    /// number of points.
    pub fn polynomial() -> InterpType {
        ffi_wrap!(gsl_interp_polynomial, gsl_interp_type)
    }

    /// Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on
    /// each interval, with matching first and second derivatives at the supplied data-points. The
    /// second derivative is chosen to be zero at the first point and last point.
    pub fn cspline() -> InterpType {
        ffi_wrap!(gsl_interp_cspline, gsl_interp_type)
    }

    /// Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic on
    /// each interval, with matching first and second derivatives at the supplied data-points. The
    /// derivatives at the first and last points are also matched. Note that the last point in the
    /// data must have the same y-value as the first point, otherwise the resulting periodic
    /// interpolation will have a discontinuity at the boundary.
    pub fn cspline_periodic() -> InterpType {
        ffi_wrap!(gsl_interp_cspline_periodic, gsl_interp_type)
    }

    /// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
    /// corner algorithm of Wodicka.
    pub fn akima() -> InterpType {
        ffi_wrap!(gsl_interp_akima, gsl_interp_type)
    }

    /// Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded
    /// corner algorithm of Wodicka.
    pub fn akima_periodic() -> InterpType {
        ffi_wrap!(gsl_interp_akima_periodic, gsl_interp_type)
    }
}

impl ffi::FFI<sys::gsl_interp_type> for InterpType {
    fn wrap(t: *mut sys::gsl_interp_type) -> InterpType {
        InterpType { t: t }
    }

    fn soft_wrap(t: *mut sys::gsl_interp_type) -> InterpType {
        Self::wrap(t)
    }

    fn unwrap_shared(t: &InterpType) -> *const sys::gsl_interp_type {
        t.t as *const sys::gsl_interp_type
    }

    fn unwrap_unique(t: &mut InterpType) -> *mut sys::gsl_interp_type {
        t.t as *mut sys::gsl_interp_type
    }
}

/// general interpolation object
pub struct Spline {
    spline: *mut sys::gsl_spline,
}

impl Spline {
    pub fn new(t: &InterpType, size: usize) -> Option<Spline> {
        let tmp = unsafe { sys::gsl_spline_alloc(t.t, size) };

        if tmp.is_null() {
            None
        } else {
            Some(Spline { spline: tmp })
        }
    }

    pub fn init(&self, xa: &[f64], ya: &[f64]) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_spline_init(self.spline, xa.as_ptr(), ya.as_ptr(), xa.len() as usize)
        })
    }

    pub fn name(&self) -> String {
        let tmp = unsafe { sys::gsl_spline_name(self.spline) };

        if tmp.is_null() {
            String::new()
        } else {
            unsafe {
                String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
            }
        }
    }

    pub fn min_size(&self) -> u32 {
        unsafe { sys::gsl_spline_min_size(self.spline) }
    }

    pub fn eval(&self, x: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval(self.spline, x, &mut acc.0) }
    }

    pub fn eval_e(&self, x: f64, acc: &mut InterpAccel, y: &mut f64) -> enums::Value {
        enums::Value::from(unsafe { sys::gsl_spline_eval_e(self.spline, x, &mut acc.0, y) })
    }

    pub fn eval_deriv(&self, x: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval_deriv(self.spline, x, &mut acc.0) }
    }

    pub fn eval_deriv_e(&self, x: f64, acc: &mut InterpAccel, d: &mut f64) -> enums::Value {
        enums::Value::from(unsafe { sys::gsl_spline_eval_deriv_e(self.spline, x, &mut acc.0, d) })
    }

    pub fn eval_deriv2(&self, x: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval_deriv2(self.spline, x, &mut acc.0) }
    }

    pub fn eval_deriv2_e(&self, x: f64, acc: &mut InterpAccel, d2: &mut f64) -> enums::Value {
        enums::Value::from(unsafe { sys::gsl_spline_eval_deriv2_e(self.spline, x, &mut acc.0, d2) })
    }

    pub fn eval_integ(&self, a: f64, b: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval_integ(self.spline, a, b, &mut acc.0) }
    }

    pub fn eval_integ_e(
        &self,
        a: f64,
        b: f64,
        acc: &mut InterpAccel,
        result: &mut f64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_spline_eval_integ_e(self.spline, a, b, &mut acc.0, result)
        })
    }
}

impl Drop for Spline {
    fn drop(&mut self) {
        unsafe { sys::gsl_spline_free(self.spline) };
        self.spline = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<sys::gsl_spline> for Spline {
    fn wrap(spline: *mut sys::gsl_spline) -> Spline {
        Spline { spline: spline }
    }

    fn soft_wrap(spline: *mut sys::gsl_spline) -> Spline {
        Self::wrap(spline)
    }

    fn unwrap_shared(spline: &Spline) -> *const sys::gsl_spline {
        spline.spline as *const _
    }

    fn unwrap_unique(spline: &mut Spline) -> *mut sys::gsl_spline {
        spline.spline
    }
}
