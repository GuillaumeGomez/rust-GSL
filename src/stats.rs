//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[doc(alias = "gsl_stats_wtss")]
pub fn wtss(w: &[f64], wstride: usize, data: &[f64], stride: usize) -> f64 {
    unsafe { sys::gsl_stats_wtss(w.as_ptr(), wstride, data.as_ptr(), stride, data.len() as _) }
}

#[doc(alias = "gsl_stats_wtss_m")]
pub fn wtss_m(w: &[f64], wstride: usize, data: &[f64], stride: usize, wmean: f64) -> f64 {
    unsafe {
        sys::gsl_stats_wtss_m(
            w.as_ptr(),
            wstride,
            data.as_ptr(),
            stride,
            data.len() as _,
            wmean,
        )
    }
}

#[doc(alias = "gsl_stats_wabsdev")]
pub fn wabsdev(w: &[f64], wstride: usize, data: &[f64], stride: usize) -> f64 {
    unsafe { sys::gsl_stats_wabsdev(w.as_ptr(), wstride, data.as_ptr(), stride, data.len() as _) }
}

#[doc(alias = "gsl_stats_wskew")]
pub fn wskew(w: &[f64], wstride: usize, data: &[f64], stride: usize) -> f64 {
    unsafe { sys::gsl_stats_wskew(w.as_ptr(), wstride, data.as_ptr(), stride, data.len() as _) }
}

#[doc(alias = "gsl_stats_wkurtosis")]
pub fn wkurtosis(w: &[f64], wstride: usize, data: &[f64], stride: usize) -> f64 {
    unsafe { sys::gsl_stats_wkurtosis(w.as_ptr(), wstride, data.as_ptr(), stride, data.len() as _) }
}

#[doc(alias = "gsl_stats_wvariance_m")]
pub fn wvariance_m(w: &[f64], wstride: usize, data: &[f64], stride: usize, wmean: f64) -> f64 {
    unsafe {
        sys::gsl_stats_wvariance_m(
            w.as_ptr(),
            wstride,
            data.as_ptr(),
            stride,
            data.len() as _,
            wmean,
        )
    }
}

#[doc(alias = "gsl_stats_wabsdev_m")]
pub fn wabsdev_m(w: &[f64], wstride: usize, data: &[f64], stride: usize, wmean: f64) -> f64 {
    unsafe {
        sys::gsl_stats_wabsdev_m(
            w.as_ptr(),
            wstride,
            data.as_ptr(),
            stride,
            data.len() as _,
            wmean,
        )
    }
}

#[doc(alias = "gsl_stats_wskew_m_sd")]
pub fn wskew_m_sd(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    wmean: f64,
    wsd: f64,
) -> f64 {
    unsafe {
        sys::gsl_stats_wskew_m_sd(
            w.as_ptr(),
            wstride,
            data.as_ptr(),
            stride,
            data.len() as _,
            wmean,
            wsd,
        )
    }
}

#[doc(alias = "gsl_stats_wkurtosis_m_sd")]
pub fn wkurtosis_m_sd(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    wmean: f64,
    wsd: f64,
) -> f64 {
    unsafe {
        sys::gsl_stats_wkurtosis_m_sd(
            w.as_ptr(),
            wstride,
            data.as_ptr(),
            stride,
            data.len() as _,
            wmean,
            wsd,
        )
    }
}

#[doc(alias = "gsl_stats_pvariance")]
pub fn pvariance(data1: &[f64], stride1: usize, data2: &[f64], stride2: usize) -> f64 {
    unsafe {
        sys::gsl_stats_pvariance(
            data1.as_ptr(),
            stride1,
            data1.len() as _,
            data2.as_ptr(),
            stride2,
            data2.len() as _,
        )
    }
}

#[doc(alias = "gsl_stats_ttest")]
pub fn ttest(data1: &[f64], stride1: usize, data2: &[f64], stride2: usize) -> f64 {
    unsafe {
        sys::gsl_stats_ttest(
            data1.as_ptr(),
            stride1,
            data1.len() as _,
            data2.as_ptr(),
            stride2,
            data2.len() as _,
        )
    }
}

#[doc(alias = "gsl_stats_max")]
pub fn max(data: &[f64], stride: usize) -> f64 {
    unsafe { sys::gsl_stats_max(data.as_ptr(), stride, data.len() as _) }
}

#[doc(alias = "gsl_stats_min")]
pub fn min(data: &[f64], stride: usize) -> f64 {
    unsafe { sys::gsl_stats_min(data.as_ptr(), stride, data.len() as _) }
}

/// Returns `(min, max)`.
pub fn gsl_stats_minmax(data: &[f64], stride: usize) -> (f64, f64) {
    let mut min = 0.;
    let mut max = 0.;

    unsafe { sys::gsl_stats_minmax(&mut min, &mut max, data.as_ptr(), stride, data.len() as _) }
    (min, max)
}

#[doc(alias = "gsl_stats_max_index")]
pub fn max_index(data: &[f64], stride: usize) -> usize {
    unsafe { sys::gsl_stats_max_index(data.as_ptr(), stride, data.len() as _) }
}

#[doc(alias = "gsl_stats_min_index")]
pub fn min_index(data: &[f64], stride: usize) -> usize {
    unsafe { sys::gsl_stats_min_index(data.as_ptr(), stride, data.len() as _) }
}

/// Returns `(min, max)`.
pub fn gsl_stats_minmax_index(data: &[f64], stride: usize) -> (usize, usize) {
    let mut min = 0;
    let mut max = 0;

    unsafe {
        sys::gsl_stats_minmax_index(&mut min, &mut max, data.as_ptr(), stride, data.len() as _)
    }
    (min, max)
}

#[cfg(feature = "v2_5")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_5")))]
#[doc(alias = "gsl_stats_select")]
pub fn select(data: &mut [f64], stride: usize, k: usize) -> f64 {
    unsafe { sys::gsl_stats_select(data.as_mut_ptr(), stride, data.len() as _, k) }
}

#[cfg(feature = "v2_5")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_5")))]
#[doc(alias = "gsl_stats_median")]
pub fn median(data: &mut [f64], stride: usize) -> f64 {
    unsafe { sys::gsl_stats_median(data.as_mut_ptr(), stride, data.len() as _) }
}
