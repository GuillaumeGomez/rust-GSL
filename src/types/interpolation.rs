//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Interpolation

This chapter describes functions for performing interpolation. The library provides a variety of interpolation methods, including Cubic splines
and Akima splines. The interpolation types are interchangeable, allowing different methods to be used without recompiling. Interpolations can
be defined for both normal and periodic boundary conditions. Additional functions are available for computing derivatives and integrals of
interpolating functions.

These interpolation methods produce curves that pass through each datapoint. To interpolate noisy data with a smoothing curve see Basis Splines.

## Introduction

Given a set of data points (x_1, y_1) \dots (x_n, y_n) the routines described in this section compute a continuous interpolating function
y(x) such that y(x_i) = y_i. The interpolation is piecewise smooth, and its behavior at the end-points is determined by the type of
interpolation used.

## Index Look-up and Acceleration

The state of searches can be stored in a gsl_interp_accel object, which is a kind of iterator for interpolation lookups. It caches the previous
value of an index lookup. When the subsequent interpolation point falls in the same interval its index value can be returned immediately.

## Higher-level Interface

The functions described in the previous sections required the user to supply pointers to the x and y arrays on each call. The following
functions are equivalent to the corresponding gsl_interp functions but maintain a copy of this data in the gsl_spline object. This removes
the need to pass both xa and ya as arguments on each evaluation.

## References and Further Reading

Descriptions of the interpolation algorithms and further references can be found in the following books:

C.W. Ueberhuber, Numerical Computation (Volume 1), Chapter 9 “Interpolation”, Springer (1997), ISBN 3-540-62058-3.
D.M. Young, R.T. Gregory A Survey of Numerical Mathematics (Volume 1), Chapter 6.8, Dover (1988), ISBN 0-486-65691-8.
!*/

use crate::Value;
use ffi::FFI;

/// Evaluation accelerator.
#[derive(Clone)]
pub struct InterpAccel(pub sys::gsl_interp_accel);

impl InterpAccel {
    /// This function returns a pointer to an accelerator object, which is a kind of iterator for
    /// interpolation lookups. It tracks the state of lookups, thus allowing for application of
    /// various acceleration strategies.
    #[allow(clippy::new_without_default)]
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
    #[doc(alias = "gsl_interp_accel_find")]
    pub fn find(&mut self, x_array: &[f64], x: f64) -> usize {
        unsafe { sys::gsl_interp_accel_find(&mut self.0, x_array.as_ptr(), x_array.len() as _, x) }
    }
}

ffi_wrapper!(Interp, *mut sys::gsl_interp, gsl_interp_free);

impl Interp {
    /// This function returns a pointer to a newly allocated interpolation object of type T for
    /// size data-points.
    ///
    /// ```
    /// use rgsl::{Interp, InterpType};
    ///
    /// let interp_type = InterpType::linear();
    /// let interp = Interp::new(interp_type, 2).expect("Failed to initialize `Interp`...");
    /// ```
    #[doc(alias = "gsl_interp_alloc")]
    pub fn new(t: InterpType, size: usize) -> Option<Interp> {
        let tmp = unsafe { sys::gsl_interp_alloc(t.unwrap_shared(), size) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function initializes the interpolation object interp for the data (xa,ya) where xa and
    /// ya are arrays of size size. The interpolation object (gsl_interp) does not save the data
    /// arrays xa and ya and only stores the static state computed from the data. The xa data array
    /// is always assumed to be strictly ordered, with increasing x values; the behavior for other
    /// arrangements is not defined.
    ///
    /// Asserts that `ya.len() >= xa.len()`.
    #[doc(alias = "gsl_interp_init")]
    pub fn init(&mut self, xa: &[f64], ya: &[f64]) -> Value {
        assert!(ya.len() >= xa.len());
        Value::from(unsafe {
            sys::gsl_interp_init(
                self.unwrap_unique(),
                xa.as_ptr(),
                ya.as_ptr(),
                xa.len() as _,
            )
        })
    }

    /// This function returns the name of the interpolation type used by interp. For example,
    ///
    /// ```
    /// use rgsl::{Interp, InterpType};
    ///
    /// let interp_type = InterpType::linear();
    /// let interp = Interp::new(interp_type, 2).expect("Failed to initialize `Interp`...");
    /// println!("interp uses '{}' interpolation.", interp.name());
    /// ```
    ///
    /// would print something like :
    ///
    /// ```Shell
    /// interp uses 'cspline' interpolation.
    /// ```
    #[doc(alias = "gsl_interp_name")]
    pub fn name(&self) -> String {
        let tmp = unsafe { sys::gsl_interp_name(self.unwrap_shared()) };

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
    #[doc(alias = "gsl_interp_min_size")]
    pub fn min_size(&self) -> u32 {
        unsafe { sys::gsl_interp_min_size(self.unwrap_shared()) }
    }
}

ffi_wrapper!(InterpType, *const sys::gsl_interp_type);

impl InterpType {
    /// This function returns the minimum number of points required by the interpolation object
    /// interp or interpolation type T. For example, Akima spline interpolation requires a minimum
    /// of 5 points.
    #[doc(alias = "gsl_interp_type_min_size")]
    pub fn min_size(&self) -> u32 {
        unsafe { sys::gsl_interp_type_min_size(self.unwrap_shared()) }
    }

    /// Linear interpolation. This interpolation method does not require any additional memory.
    pub fn linear() -> InterpType {
        ffi_wrap!(gsl_interp_linear)
    }

    /// Polynomial interpolation. This method should only be used for interpolating small numbers
    /// of points because polynomial interpolation introduces large oscillations, even for
    /// well-behaved datasets. The number of terms in the interpolating polynomial is equal to the
    /// number of points.
    pub fn polynomial() -> InterpType {
        ffi_wrap!(gsl_interp_polynomial)
    }

    /// Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on
    /// each interval, with matching first and second derivatives at the supplied data-points. The
    /// second derivative is chosen to be zero at the first point and last point.
    pub fn cspline() -> InterpType {
        ffi_wrap!(gsl_interp_cspline)
    }

    /// Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic on
    /// each interval, with matching first and second derivatives at the supplied data-points. The
    /// derivatives at the first and last points are also matched. Note that the last point in the
    /// data must have the same y-value as the first point, otherwise the resulting periodic
    /// interpolation will have a discontinuity at the boundary.
    pub fn cspline_periodic() -> InterpType {
        ffi_wrap!(gsl_interp_cspline_periodic)
    }

    /// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
    /// corner algorithm of Wodicka.
    pub fn akima() -> InterpType {
        ffi_wrap!(gsl_interp_akima)
    }

    /// Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded
    /// corner algorithm of Wodicka.
    pub fn akima_periodic() -> InterpType {
        ffi_wrap!(gsl_interp_akima_periodic)
    }
}

ffi_wrapper!(
    Spline,
    *mut sys::gsl_spline,
    gsl_spline_free,
    "General interpolation object."
);

impl Spline {
    #[doc(alias = "gsl_spline_alloc")]
    pub fn new(t: InterpType, size: usize) -> Option<Spline> {
        let tmp = unsafe { sys::gsl_spline_alloc(t.unwrap_shared(), size) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    #[doc(alias = "gsl_spline_init")]
    pub fn init(&mut self, xa: &[f64], ya: &[f64]) -> Value {
        Value::from(unsafe {
            sys::gsl_spline_init(
                self.unwrap_unique(),
                xa.as_ptr(),
                ya.as_ptr(),
                xa.len() as _,
            )
        })
    }

    #[doc(alias = "gsl_spline_name")]
    pub fn name(&self) -> String {
        let tmp = unsafe { sys::gsl_spline_name(self.unwrap_shared()) };

        if tmp.is_null() {
            String::new()
        } else {
            unsafe {
                String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()
            }
        }
    }

    #[doc(alias = "gsl_spline_min_size")]
    pub fn min_size(&self) -> u32 {
        unsafe { sys::gsl_spline_min_size(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_spline_eval")]
    pub fn eval(&self, x: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval(self.unwrap_shared(), x, &mut acc.0) }
    }

    /// Returns `(Value, y)`.
    #[doc(alias = "gsl_spline_eval_e")]
    pub fn eval_e(&self, x: f64, acc: &mut InterpAccel) -> (Value, f64) {
        let mut y = 0.;
        let ret = unsafe { sys::gsl_spline_eval_e(self.unwrap_shared(), x, &mut acc.0, &mut y) };
        (Value::from(ret), y)
    }

    #[doc(alias = "gsl_spline_eval_deriv")]
    pub fn eval_deriv(&self, x: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval_deriv(self.unwrap_shared(), x, &mut acc.0) }
    }

    /// Returns `(Value, d)`.
    #[doc(alias = "gsl_spline_eval_deriv_e")]
    pub fn eval_deriv_e(&self, x: f64, acc: &mut InterpAccel) -> (Value, f64) {
        let mut d = 0.;
        let ret =
            unsafe { sys::gsl_spline_eval_deriv_e(self.unwrap_shared(), x, &mut acc.0, &mut d) };
        (Value::from(ret), d)
    }

    #[doc(alias = "gsl_spline_eval_deriv2")]
    pub fn eval_deriv2(&self, x: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval_deriv2(self.unwrap_shared(), x, &mut acc.0) }
    }

    /// Returns `(Value, d2)`.
    #[doc(alias = "gsl_spline_eval_deriv2_e")]
    pub fn eval_deriv2_e(&self, x: f64, acc: &mut InterpAccel) -> (Value, f64) {
        let mut d2 = 0.;
        let ret =
            unsafe { sys::gsl_spline_eval_deriv2_e(self.unwrap_shared(), x, &mut acc.0, &mut d2) };
        (Value::from(ret), d2)
    }

    #[doc(alias = "gsl_spline_eval_integ")]
    pub fn eval_integ(&self, a: f64, b: f64, acc: &mut InterpAccel) -> f64 {
        unsafe { sys::gsl_spline_eval_integ(self.unwrap_shared(), a, b, &mut acc.0) }
    }

    /// Returns `(Value, d2)`.
    #[doc(alias = "gsl_spline_eval_integ_e")]
    pub fn eval_integ_e(&self, a: f64, b: f64, acc: &mut InterpAccel) -> (Value, f64) {
        let mut result = 0.;
        let ret = unsafe {
            sys::gsl_spline_eval_integ_e(self.unwrap_shared(), a, b, &mut acc.0, &mut result)
        };
        (Value::from(ret), result)
    }
}
