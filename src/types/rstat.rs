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
    pub fn reset(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_rstat_quantile_reset(self.unwrap_unique()) })
    }

    pub fn add(&mut self, x: f64) -> Value {
        Value::from(unsafe { sys::gsl_rstat_quantile_add(x, self.unwrap_unique()) })
    }

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
    pub fn new() -> Option<Self> {
        let s = unsafe { sys::gsl_rstat_alloc() };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    pub fn n(&self) -> usize {
        unsafe { sys::gsl_rstat_n(self.unwrap_shared()) }
    }

    pub fn add(&mut self, x: f64) -> Value {
        Value::from(unsafe { sys::gsl_rstat_add(x, self.unwrap_unique()) })
    }

    pub fn min(&self) -> f64 {
        unsafe { sys::gsl_rstat_min(self.unwrap_shared()) }
    }

    pub fn max(&self) -> f64 {
        unsafe { sys::gsl_rstat_max(self.unwrap_shared()) }
    }

    pub fn mean(&self) -> f64 {
        unsafe { sys::gsl_rstat_mean(self.unwrap_shared()) }
    }

    pub fn variance(&self) -> f64 {
        unsafe { sys::gsl_rstat_variance(self.unwrap_shared()) }
    }

    pub fn sd(&self) -> f64 {
        unsafe { sys::gsl_rstat_sd(self.unwrap_shared()) }
    }

    #[cfg(feature = "v2_2")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
    pub fn rms(&self) -> f64 {
        unsafe { sys::gsl_rstat_rms(self.unwrap_shared()) }
    }

    pub fn sd_mean(&self) -> f64 {
        unsafe { sys::gsl_rstat_sd_mean(self.unwrap_shared()) }
    }

    pub fn median(&mut self) -> f64 {
        unsafe { sys::gsl_rstat_median(self.unwrap_unique()) }
    }

    pub fn skew(&self) -> f64 {
        unsafe { sys::gsl_rstat_skew(self.unwrap_shared()) }
    }

    pub fn kurtosis(&self) -> f64 {
        unsafe { sys::gsl_rstat_kurtosis(self.unwrap_shared()) }
    }

    pub fn reset(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_rstat_reset(self.unwrap_unique()) })
    }
}
