//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::Value;
use ffi::FFI;

ffi_wrapper!(
    RStatQuantileWorkspace,
    *mut sys::gsl_rstat_quantile_workspace,
    gsl_rstat_quantile_free
);

impl RStatQuantileWorkspace {
    #[doc(alias = "gsl_rstat_quantile_alloc")]
    pub fn new(p: f64) -> Option<Self> {
        let s = unsafe { sys::gsl_rstat_quantile_alloc(p) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    #[cfg(feature = "v2_2")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
    #[doc(alias = "gsl_rstat_quantile_reset")]
    pub fn reset(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_rstat_quantile_reset(self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_rstat_quantile_add")]
    pub fn add(&mut self, x: f64) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_rstat_quantile_add(x, self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_rstat_quantile_get")]
    pub fn get(&mut self) -> f64 {
        unsafe { sys::gsl_rstat_quantile_get(self.unwrap_unique()) }
    }
}

ffi_wrapper!(
    RStatWorkspace,
    *mut sys::gsl_rstat_workspace,
    gsl_rstat_free
);

impl RStatWorkspace {
    #[doc(alias = "gsl_rstat_alloc")]
    pub fn new() -> Option<Self> {
        let s = unsafe { sys::gsl_rstat_alloc() };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    #[doc(alias = "gsl_rstat_n")]
    pub fn n(&self) -> usize {
        unsafe { sys::gsl_rstat_n(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_add")]
    pub fn add(&mut self, x: f64) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_rstat_add(x, self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_rstat_min")]
    pub fn min(&self) -> f64 {
        unsafe { sys::gsl_rstat_min(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_max")]
    pub fn max(&self) -> f64 {
        unsafe { sys::gsl_rstat_max(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_mean")]
    pub fn mean(&self) -> f64 {
        unsafe { sys::gsl_rstat_mean(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_variance")]
    pub fn variance(&self) -> f64 {
        unsafe { sys::gsl_rstat_variance(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_sd")]
    pub fn sd(&self) -> f64 {
        unsafe { sys::gsl_rstat_sd(self.unwrap_shared()) }
    }

    #[cfg(feature = "v2_2")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
    #[doc(alias = "gsl_rstat_rms")]
    pub fn rms(&self) -> f64 {
        unsafe { sys::gsl_rstat_rms(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_sd_mean")]
    pub fn sd_mean(&self) -> f64 {
        unsafe { sys::gsl_rstat_sd_mean(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_median")]
    pub fn median(&mut self) -> f64 {
        unsafe { sys::gsl_rstat_median(self.unwrap_unique()) }
    }

    #[doc(alias = "gsl_rstat_skew")]
    pub fn skew(&self) -> f64 {
        unsafe { sys::gsl_rstat_skew(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_kurtosis")]
    pub fn kurtosis(&self) -> f64 {
        unsafe { sys::gsl_rstat_kurtosis(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_rstat_reset")]
    pub fn reset(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_rstat_reset(self.unwrap_unique()) };
        result_handler!(ret, ())
    }
}
