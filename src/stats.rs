//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::vector::{self, Vector};

#[cfg(feature = "v2_5")]
use crate::vector::VectorMut;

// FIXME: Many functions are missing.

/// # Weighted Samples
///
/// The functions described in this section allow the computation of
/// statistics for weighted samples.  The functions accept a vector of
/// samples, xᵢ, with associated weights, wᵢ.  Each sample xᵢ is
/// considered as having been drawn from a Gaussian distribution with
/// variance σᵢ².  The sample weight wᵢ is defined as the reciprocal
/// of this variance, wᵢ = 1/σᵢ².  Setting a weight to zero
/// corresponds to removing a sample from a dataset.

/// Return the weighted mean of the dataset `data` using the set of
/// weights `w`. The weighted mean is defined as,
/// ̂μ = (∑ wᵢ xᵢ) / (∑ wᵢ).
///
/// # Example
///
/// ```
/// use rgsl::{stats::wmean, vector::Vector};
/// let m = wmean(&[1., 1.], &[1., 1.]);
/// assert_eq!(m, 1.);
/// ```
#[doc(alias = "gsl_stats_wmean")]
pub fn wmean<T>(w: &T, data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wmean: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wmean(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
}

/// Returns the estimated variance of the weighted dataset `data`
/// using the set of weights `w`.  The estimated variance of a
/// weighted dataset is calculated as,
/// ̂σ² = (∑ wᵢ) / ((∑ wᵢ)² - ∑ wᵢ²) · ∑ wᵢ (xᵢ - ̂μ)².
#[doc(alias = "gsl_stats_wvariance")]
pub fn wvariance<T>(w: &T, data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wvariance: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wvariance(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
}

/// Returns the estimated variance of the weighted dataset `data`
/// using the given weighted mean `wmean`.
#[doc(alias = "gsl_stats_wvariance_m")]
pub fn wvariance_m<T>(w: &T, data: &T, wmean: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wvariance_m: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wvariance_m(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            wmean,
        )
    }
}

/// Return the standard deviation is defined as the square root of the
/// variance.
#[doc(alias = "gsl_stats_wsd")]
pub fn wsd<T>(w: &T, data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wsd: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wsd(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
}

/// Return the standard deviation is defined as the square root of the
/// variance using the given weighted mean `wmean`.
#[doc(alias = "gsl_stats_wsd_m")]
pub fn wsd_m<T>(w: &T, data: &T, wmean: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wsd_m: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wsd_m(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            wmean,
        )
    }
}

/// Return an unbiased estimate of the variance of the weighted
/// dataset `data` when the population mean `mean` of the underlying
/// distribution is known a priori.  In this case the estimator for
/// the variance replaces the sample mean ̂μ by the known population
/// mean μ:
/// σ² = ∑ wᵢ (xᵢ - μ)² / (∑ wᵢ)..
#[doc(alias = "gsl_stats_wvariance_with_fixed_mean")]
pub fn wvariance_with_fixed_mean<T>(w: &T, data: &T, mean: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wvariance_with_fixed_mean: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wvariance_with_fixed_mean(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            mean,
        )
    }
}

/// Return the standard deviation which is defined as the square root
/// of the variance computed by [`wvariance_with_fixed_mean`].
#[doc(alias = "gsl_stats_wsd_with_fixed_mean")]
pub fn wsd_with_fixed_mean<T>(w: &T, data: &T, mean: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wsd_with_fixed_mean: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wsd_with_fixed_mean(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            mean,
        )
    }
}

/// Return the weighted total sum of squares (TSS) of data about the
/// weighted mean.  TSS = ∑ wᵢ (xᵢ - wmean)² where the weighted mean
/// wmean is computed internally.
#[doc(alias = "gsl_stats_wtss")]
pub fn wtss<T>(w: &T, data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wtss: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wtss(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
}

/// Return the weighted total sum of squares (TSS) of data about the
/// weighted mean.  TSS = ∑ wᵢ (xᵢ - `wmean`)².
#[doc(alias = "gsl_stats_wtss_m")]
pub fn wtss_m<T>(w: &T, data: &T, wmean: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wtss_m: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wtss_m(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            wmean,
        )
    }
}

/// Return the weighted absolute deviation from the weighted mean of
/// data.  The absolute deviation from the mean is defined as,
/// absdev = (∑ wᵢ |xᵢ - ̂μ|) / (∑ wᵢ)
#[doc(alias = "gsl_stats_wabsdev")]
pub fn wabsdev<T>(w: &T, data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wabsdev: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wabsdev(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
}

/// Return the absolute deviation of the weighted dataset data about
/// the given weighted mean `wmean`.
#[doc(alias = "gsl_stats_wabsdev_m")]
pub fn wabsdev_m<T>(w: &T, data: &T, wmean: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wabsdev_m: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wabsdev_m(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            wmean,
        )
    }
}

/// Return the weighted skewness of the dataset `data`.
/// skew = (∑ wᵢ ((xᵢ - ̂x) / ̂σ)³) / (∑ wᵢ)
#[doc(alias = "gsl_stats_wskew")]
pub fn wskew<T>(w: &T, data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wskew: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wskew(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
}

/// Return the weighted skewness of the dataset `data` using the given
/// values of the weighted mean and weighted standard deviation,
/// `wmean` and `wsd`.
#[doc(alias = "gsl_stats_wskew_m_sd")]
pub fn wskew_m_sd<T>(w: &T, data: &T, wmean: f64, wsd: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wskew_m_sd: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wskew_m_sd(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            wmean,
            wsd,
        )
    }
}

/// Return the weighted kurtosis of the dataset `data`.
/// kurtosis = (∑ wᵢ ((xᵢ - ̂x) / ̂σ)⁴) / (∑ wᵢ) - 3
#[doc(alias = "gsl_stats_wkurtosis")]
pub fn wkurtosis<T>(w: &T, data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wkurtosis: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wkurtosis(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
}

/// Return the weighted kurtosis of the dataset `data` using the given
/// values of the weighted mean and weighted standard deviation,
/// `wmean` and `wsd`.
#[doc(alias = "gsl_stats_wkurtosis_m_sd")]
pub fn wkurtosis_m_sd<T>(w: &T, data: &T, wmean: f64, wsd: f64) -> f64
where
    T: Vector<f64> + ?Sized,
{
    if T::len(w) != T::len(data) {
        panic!("rgsl::stats::wkurtosis_m_sd: the size of w and data must be the same");
    }
    unsafe {
        sys::gsl_stats_wkurtosis_m_sd(
            vector::as_ptr(w),
            T::stride(w),
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
            wmean,
            wsd,
        )
    }
}

#[doc(alias = "gsl_stats_pvariance")]
pub fn pvariance<T>(data1: &T, data2: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    unsafe {
        sys::gsl_stats_pvariance(
            vector::as_ptr(data1),
            T::stride(data1),
            T::len(data1),
            vector::as_ptr(data2),
            T::stride(data2),
            T::len(data2),
        )
    }
}

#[doc(alias = "gsl_stats_ttest")]
pub fn ttest<T>(data1: &T, data2: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    unsafe {
        sys::gsl_stats_ttest(
            vector::as_ptr(data1),
            T::stride(data1),
            T::len(data1),
            vector::as_ptr(data2),
            T::stride(data2),
            T::len(data2),
        )
    }
}

#[doc(alias = "gsl_stats_max")]
pub fn max<T>(data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    unsafe { sys::gsl_stats_max(vector::as_ptr(data), T::stride(data), T::len(data)) }
}

#[doc(alias = "gsl_stats_min")]
pub fn min<T>(data: &T) -> f64
where
    T: Vector<f64> + ?Sized,
{
    unsafe { sys::gsl_stats_min(vector::as_ptr(data), T::stride(data), T::len(data)) }
}

/// Returns `(min, max)`.
#[doc(alias = "gsl_stats_minmax")]
pub fn stats_minmax<T>(data: &T) -> (f64, f64)
where
    T: Vector<f64> + ?Sized,
{
    let mut min = 0.;
    let mut max = 0.;

    unsafe {
        sys::gsl_stats_minmax(
            &mut min,
            &mut max,
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
    (min, max)
}

#[doc(alias = "gsl_stats_max_index")]
pub fn max_index<T>(data: &T) -> usize
where
    T: Vector<f64> + ?Sized,
{
    unsafe { sys::gsl_stats_max_index(vector::as_ptr(data), T::stride(data), T::len(data)) }
}

#[doc(alias = "gsl_stats_min_index")]
pub fn min_index<T>(data: &T) -> usize
where
    T: Vector<f64> + ?Sized,
{
    unsafe { sys::gsl_stats_min_index(vector::as_ptr(data), T::stride(data), T::len(data)) }
}

/// Returns `(min, max)`.
#[doc(alias = "gsl_stats_minmax_index")]
pub fn stats_minmax_index<T>(data: &T) -> (usize, usize)
where
    T: Vector<f64> + ?Sized,
{
    let mut min = 0;
    let mut max = 0;

    unsafe {
        sys::gsl_stats_minmax_index(
            &mut min,
            &mut max,
            vector::as_ptr(data),
            T::stride(data),
            T::len(data),
        )
    }
    (min, max)
}

#[cfg(feature = "v2_5")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_5")))]
#[doc(alias = "gsl_stats_select")]
pub fn select<T>(data: &mut T, k: usize) -> f64
where
    T: VectorMut<f64> + ?Sized,
{
    unsafe { sys::gsl_stats_select(vector::as_mut_ptr(data), T::stride(data), T::len(data), k) }
}

#[cfg(feature = "v2_5")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_5")))]
#[doc(alias = "gsl_stats_median")]
pub fn median<T>(data: &mut T) -> f64
where
    T: VectorMut<f64> + ?Sized,
{
    unsafe { sys::gsl_stats_median(vector::as_mut_ptr(data), T::stride(data), T::len(data)) }
}
