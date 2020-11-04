//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{FilterEnd, FilterScale, Value, VectorF64, VectorI32};
use ffi::FFI;

ffi_wrapper!(
    FilterGaussianWorkspace,
    *mut sys::gsl_filter_gaussian_workspace,
    gsl_filter_gaussian_free
);

impl FilterGaussianWorkspace {
    pub fn alloc(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_gaussian_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    /// This function applies a Gaussian filter parameterized by `alpha` to the input vector `x`,
    /// storing the output in `y`. The derivative order is specified by `order`, with `0`
    /// corresponding to a Gaussian, `1` corresponding to a first derivative Gaussian, and so on.
    /// The parameter `endtype` specifies how the signal end points are handled. It is allowed for
    /// `x` = `y` for an in-place filter.
    pub fn gaussian(
        &mut self,
        endtype: FilterEnd,
        alpha: f64,
        order: usize,
        x: &VectorF64,
        y: &mut VectorF64,
    ) -> Value {
        Value::from(unsafe {
            sys::gsl_filter_gaussian(
                endtype.into(),
                alpha,
                order,
                x.unwrap_shared(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            )
        })
    }
}

ffi_wrapper!(
    FilterMedianWorkspace,
    *mut sys::gsl_filter_median_workspace,
    gsl_filter_median_free
);

impl FilterMedianWorkspace {
    pub fn alloc(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_median_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    pub fn median(&mut self, endtype: FilterEnd, x: &VectorF64, y: &mut VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_filter_median(
                endtype.into(),
                x.unwrap_shared(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            )
        })
    }
}

ffi_wrapper!(
    FilterRMedianWorkspace,
    *mut sys::gsl_filter_rmedian_workspace,
    gsl_filter_rmedian_free
);

impl FilterRMedianWorkspace {
    pub fn alloc(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_rmedian_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    pub fn rmedian(&mut self, endtype: FilterEnd, x: &VectorF64, y: &mut VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_filter_rmedian(
                endtype.into(),
                x.unwrap_shared(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            )
        })
    }
}

ffi_wrapper!(
    FilterImpulseWorkspace,
    *mut sys::gsl_filter_impulse_workspace,
    gsl_filter_impulse_free
);

impl FilterImpulseWorkspace {
    pub fn alloc(K: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_filter_impulse_alloc(K) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    pub fn impulse(
        &mut self,
        endtype: FilterEnd,
        scale_type: FilterScale,
        t: f64,
        x: &VectorF64,
        y: &mut VectorF64,
        xmedian: &mut VectorF64,
        xsigma: &mut VectorF64,
        noutlier: &mut usize,
        ioutlier: &mut VectorI32,
    ) -> Value {
        Value::from(unsafe {
            sys::gsl_filter_impulse(
                endtype.into(),
                scale_type.into(),
                t,
                x.unwrap_shared(),
                y.unwrap_unique(),
                xmedian.unwrap_unique(),
                xsigma.unwrap_unique(),
                noutlier,
                ioutlier.unwrap_unique(),
                self.unwrap_unique(),
            )
        })
    }
}
