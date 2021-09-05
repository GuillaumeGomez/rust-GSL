//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Statistics

This chapter describes the statistical functions in the library. The basic statistical functions include routines to compute the mean,
variance and standard deviation. More advanced functions allow you to calculate absolute deviations, skewness, and kurtosis as well as the
median and arbitrary percentiles. The algorithms use recurrence relations to compute average quantities in a stable way, without large
intermediate values that might overflow.

## Weighted Samples

The functions described in this section allow the computation of statistics for weighted samples. The functions accept an array of
samples, x_i, with associated weights, w_i. Each sample x_i is considered as having been drawn from a Gaussian distribution with variance
\sigma_i^2. The sample weight w_i is defined as the reciprocal of this variance, w_i = 1/\sigma_i^2. Setting a weight to zero corresponds
to removing a sample from a dataset.

## Maximum and Minimum values

The following functions find the maximum and minimum values of a dataset (or their indices). If the data contains NaNs then a NaN will be
returned, since the maximum or minimum value is undefined. For functions which return an index, the location of the first NaN in the array is returned.

## Median and Percentiles

The median and percentile functions described in this section operate on sorted data. For convenience we use quantiles, measured on a
scale of 0 to 1, instead of percentiles (which use a scale of 0 to 100).

## References and Further Reading

The standard reference for almost any topic in statistics is the multi-volume Advanced Theory of Statistics by Kendall and Stuart.

Maurice Kendall, Alan Stuart, and J. Keith Ord. The Advanced Theory of Statistics (multiple volumes) reprinted as Kendall’s Advanced
Theory of Statistics. Wiley, ISBN 047023380X.
Many statistical concepts can be more easily understood by a Bayesian approach. The following book by Gelman, Carlin, Stern and Rubin
gives a comprehensive coverage of the subject.

Andrew Gelman, John B. Carlin, Hal S. Stern, Donald B. Rubin. Bayesian Data Analysis. Chapman & Hall, ISBN 0412039915.
For physicists the Particle Data Group provides useful reviews of Probability and Statistics in the “Mathematical Tools” section of its
Annual Review of Particle Physics.

Review of Particle Properties R.M. Barnett et al., Physical Review D54, 1 (1996)
The Review of Particle Physics is available online at the website http://pdg.lbl.gov/.
!*/

/// This function returns the arithmetic mean of data, a dataset of length n with stride stride. The
/// arithmetic mean, or sample mean, is denoted by \Hat\mu and defined as,
///
/// \Hat\mu = (1/N) \sum x_i
///
/// where x_i are the elements of the dataset data. For samples drawn from a gaussian distribution
/// the variance of \Hat\mu is \sigma^2 / N.
#[doc(alias = "gsl_stats_mean")]
pub fn mean(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_mean(data.as_ptr(), stride, n) }
}

/// This function returns the estimated, or sample, variance of data, a dataset of length n with
/// stride stride. The estimated variance is denoted by \Hat\sigma^2 and is defined by,
///
/// \Hat\sigma^2 = (1/(N-1)) \sum (x_i - \Hat\mu)^2
///
/// where x_i are the elements of the dataset data. Note that the normalization factor of 1/(N-1)
/// results from the derivation of \Hat\sigma^2 as an unbiased estimator of the population variance
/// \sigma^2. For samples drawn from a Gaussian distribution the variance of \Hat\sigma^2 itself is
/// 2 \sigma^4 / N.
///
/// This function computes the mean via a call to gsl_stats_mean. If you have already computed the
/// mean then you can pass it directly to gsl_stats_variance_m.
#[doc(alias = "gsl_stats_variance")]
pub fn variance(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_variance(data.as_ptr(), stride, n) }
}

/// This function returns the sample variance of data relative to the given value of mean. The
/// function is computed with \Hat\mu replaced by the value of mean that you supply,
///
/// \Hat\sigma^2 = (1/(N-1)) \sum (x_i - mean)^2
#[doc(alias = "gsl_stats_variance_m")]
pub fn variance_m(data: &[f64], stride: usize, n: usize, mean: f64) -> f64 {
    unsafe { sys::gsl_stats_variance_m(data.as_ptr(), stride, n, mean) }
}

/// The standard deviation is defined as the square root of the variance. This function returns the
/// square root of the corresponding variance functions above.
#[doc(alias = "gsl_stats_sd")]
pub fn sd(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_sd(data.as_ptr(), stride, n) }
}

/// The standard deviation is defined as the square root of the variance. This function returns the
/// square root of the corresponding variance functions above.
#[doc(alias = "gsl_stats_sd_m")]
pub fn sd_m(data: &[f64], stride: usize, n: usize, mean: f64) -> f64 {
    unsafe { sys::gsl_stats_sd_m(data.as_ptr(), stride, n, mean) }
}

/// This function returns the total sum of squares (TSS) of data about the mean. For gsl_stats_tss_m
/// the user-supplied value of mean is used, and for gsl_stats_tss it is computed using
/// gsl_stats_mean.
///
/// TSS =  \sum (x_i - mean)^2
#[doc(alias = "gsl_stats_tss")]
pub fn tss(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_tss(data.as_ptr(), stride, n) }
}

/// This function returns the total sum of squares (TSS) of data about the mean. For gsl_stats_tss_m
/// the user-supplied value of mean is used, and for gsl_stats_tss it is computed using
/// gsl_stats_mean.
///
/// TSS =  \sum (x_i - mean)^2
#[doc(alias = "gsl_stats_tss_m")]
pub fn tss_m(data: &[f64], stride: usize, n: usize, mean: f64) -> f64 {
    unsafe { sys::gsl_stats_tss_m(data.as_ptr(), stride, n, mean) }
}

/// This function computes an unbiased estimate of the variance of data when the population mean
/// mean of the underlying distribution is known a priori. In this case the estimator for the
/// variance uses the factor 1/N and the sample mean \Hat\mu is replaced by the known population
/// mean \mu,
///
/// \Hat\sigma^2 = (1/N) \sum (x_i - \mu)^2
#[doc(alias = "gsl_stats_variance_with_fixed_mean")]
pub fn variance_with_fixed_mean(data: &[f64], stride: usize, n: usize, mean: f64) -> f64 {
    unsafe { sys::gsl_stats_variance_with_fixed_mean(data.as_ptr(), stride, n, mean) }
}

/// This function calculates the standard deviation of data for a fixed population mean mean. The
/// result is the square root of the corresponding variance function.
#[doc(alias = "gsl_stats_sd_with_fixed_mean")]
pub fn sd_with_fixed_mean(data: &[f64], stride: usize, n: usize, mean: f64) -> f64 {
    unsafe { sys::gsl_stats_sd_with_fixed_mean(data.as_ptr(), stride, n, mean) }
}

/// This function computes the absolute deviation from the mean of data, a dataset of length n with
/// stride stride. The absolute deviation from the mean is defined as,
///
/// absdev  = (1/N) \sum |x_i - \Hat\mu|
///
/// where x_i are the elements of the dataset data. The absolute deviation from the mean provides a
/// more robust measure of the width of a distribution than the variance. This function computes the
/// mean of data via a call to gsl_stats_mean.
#[doc(alias = "gsl_stats_absdev")]
pub fn absdev(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_absdev(data.as_ptr(), stride, n) }
}

/// This function computes the absolute deviation of the dataset data relative to the given value of
/// mean,
///
/// absdev  = (1/N) \sum |x_i - mean|
///
/// This function is useful if you have already computed the mean of data (and want to avoid
/// recomputing it), or wish to calculate the absolute deviation relative to another value (such as
/// zero, or the median).
#[doc(alias = "gsl_stats_absdev_m")]
pub fn absdev_m(data: &[f64], stride: usize, n: usize, mean: f64) -> f64 {
    unsafe { sys::gsl_stats_absdev_m(data.as_ptr(), stride, n, mean) }
}

/// This function computes the skewness of data, a dataset of length n with stride stride. The
/// skewness is defined as,
///
/// skew = (1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^3
///
/// where x_i are the elements of the dataset data. The skewness measures the asymmetry of the tails
/// of a distribution.
///
/// The function computes the mean and estimated standard deviation of data via calls to [`mean`]
/// and [`sd`].
#[doc(alias = "gsl_stats_skew")]
pub fn skew(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_skew(data.as_ptr(), stride, n) }
}

/// This function computes the skewness of the dataset data using the given values of the mean mean
/// and standard deviation sd,
///
/// skew = (1/N) \sum ((x_i - mean)/sd)^3
///
/// These functions are useful if you have already computed the mean and standard deviation of data
/// and want to avoid recomputing them.
#[doc(alias = "gsl_stats_skew_m_sd")]
pub fn skew_m_sd(data: &[f64], stride: usize, n: usize, mean: f64, sd: f64) -> f64 {
    unsafe { sys::gsl_stats_skew_m_sd(data.as_ptr(), stride, n, mean, sd) }
}

/// This function computes the kurtosis of data, a dataset of length n with stride stride. The
/// kurtosis is defined as,
///
/// kurtosis = ((1/N) \sum ((x_i - \Hat\mu)/\Hat\sigma)^4)  - 3
///
/// The kurtosis measures how sharply peaked a distribution is, relative to its width. The kurtosis
/// is normalized to zero for a Gaussian distribution.
#[doc(alias = "gsl_stats_kurtosis")]
pub fn kurtosis(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_kurtosis(data.as_ptr(), stride, n) }
}

/// This function computes the kurtosis of the dataset data using the given values of the mean mean
/// and standard deviation sd,
///
/// kurtosis = ((1/N) \sum ((x_i - mean)/sd)^4) - 3
///
/// This function is useful if you have already computed the mean and standard deviation of data and
/// want to avoid recomputing them.
#[doc(alias = "gsl_stats_kurtosis_m_sd")]
pub fn kurtosis_m_sd(data: &[f64], stride: usize, n: usize, mean: f64, sd: f64) -> f64 {
    unsafe { sys::gsl_stats_kurtosis_m_sd(data.as_ptr(), stride, n, mean, sd) }
}

/// This function computes the lag-1 autocorrelation of the dataset data.
///
/// a_1 = {\sum_{i = 1}^{n} (x_{i} - \Hat\mu) (x_{i-1} - \Hat\mu)
///        \over
///        \sum_{i = 1}^{n} (x_{i} - \Hat\mu) (x_{i} - \Hat\mu)}
#[doc(alias = "gsl_stats_lag1_autocorrelation")]
pub fn lag1_autocorrelation(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_lag1_autocorrelation(data.as_ptr(), stride, n) }
}

/// This function computes the lag-1 autocorrelation of the dataset data using the given value of
/// the mean mean.
#[doc(alias = "gsl_stats_lag1_autocorrelation_m")]
pub fn lag1_autocorrelation_m(data: &[f64], stride: usize, n: usize, mean: f64) -> f64 {
    unsafe { sys::gsl_stats_lag1_autocorrelation_m(data.as_ptr(), stride, n, mean) }
}

/// This function computes the covariance of the datasets data1 and data2 which must both be of the
/// same length n.
///
/// covar = (1/(n - 1)) \sum_{i = 1}^{n} (x_i - \Hat x) (y_i - \Hat y)
#[doc(alias = "gsl_stats_covariance")]
pub fn covariance(data1: &[f64], stride1: usize, data2: &[f64], stride2: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_covariance(data1.as_ptr(), stride1, data2.as_ptr(), stride2, n) }
}

/// This function computes the covariance of the datasets data1 and data2 using the given values of
/// the means, mean1 and mean2. This is useful if you have already computed the means of data1 and
/// data2 and want to avoid recomputing them.
#[doc(alias = "gsl_stats_covariance_m")]
pub fn covariance_m(
    data1: &[f64],
    stride1: usize,
    data2: &[f64],
    stride2: usize,
    n: usize,
    mean1: f64,
    mean2: f64,
) -> f64 {
    unsafe {
        sys::gsl_stats_covariance_m(
            data1.as_ptr(),
            stride1,
            data2.as_ptr(),
            stride2,
            n,
            mean1,
            mean2,
        )
    }
}

/// This function efficiently computes the Pearson correlation coefficient between the datasets
/// data1 and data2 which must both be of the same length n.
///
/// r = cov(x, y) / (\Hat\sigma_x \Hat\sigma_y)
///   = {1/(n-1) \sum (x_i - \Hat x) (y_i - \Hat y)
///      \over
///      \sqrt{1/(n-1) \sum (x_i - \Hat x)^2} \sqrt{1/(n-1) \sum (y_i - \Hat y)^2}
///     }
#[doc(alias = "gsl_stats_correlation")]
pub fn correlation(data1: &[f64], stride1: usize, data2: &[f64], stride2: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_correlation(data1.as_ptr(), stride1, data2.as_ptr(), stride2, n) }
}

/// This function computes the Spearman rank correlation coefficient between the datasets data1 and
/// data2 which must both be of the same length n. Additional workspace of size 2*n is required in
/// work. The Spearman rank correlation between vectors x and y is equivalent to the Pearson
/// correlation between the ranked vectors x_R and y_R, where ranks are defined to be the average of
/// the positions of an element in the ascending order of the values.
#[doc(alias = "gsl_stats_spearman")]
pub fn spearman(
    data1: &[f64],
    stride1: usize,
    data2: &[f64],
    stride2: usize,
    n: usize,
    work: &mut [f64],
) -> f64 {
    unsafe {
        sys::gsl_stats_spearman(
            data1.as_ptr(),
            stride1,
            data2.as_ptr(),
            stride2,
            n,
            work.as_mut_ptr(),
        )
    }
}

/// This function returns the weighted mean of the dataset data with stride stride and length n,
/// using the set of weights w with stride wstride and length n. The weighted mean is defined as,
///
/// \Hat\mu = (\sum w_i x_i) / (\sum w_i)
#[doc(alias = "gsl_stats_wmean")]
pub fn wmean(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_wmean(w.as_ptr(), wstride, data.as_ptr(), stride, n) }
}

/// This function returns the estimated variance of the dataset data with stride stride and length
/// n, using the set of weights w with
/// stride wstride and length n. The estimated variance of a weighted dataset is calculated as,
///
/// \Hat\sigma^2 = ((\sum w_i)/((\sum w_i)^2 - \sum (w_i^2)))
///                 \sum w_i (x_i - \Hat\mu)^2
///
/// Note that this expression reduces to an unweighted variance with the familiar 1/(N-1) factor
/// when there are N equal non-zero weights.
#[doc(alias = "gsl_stats_wvariance")]
pub fn wvariance(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_wvariance(w.as_ptr(), wstride, data.as_ptr(), stride, n) }
}

/// This function returns the estimated variance of the weighted dataset data using the given
/// weighted mean wmean.
#[doc(alias = "gsl_stats_wvariance_m")]
pub fn wvariance_m(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    n: usize,
    wmean: f64,
) -> f64 {
    unsafe { sys::gsl_stats_wvariance_m(w.as_ptr(), wstride, data.as_ptr(), stride, n, wmean) }
}

/// The standard deviation is defined as the square root of the variance. This function returns the
/// square root of the corresponding variance function [`wvariance`] above.
#[doc(alias = "gsl_stats_wsd")]
pub fn wsd(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_wsd(w.as_ptr(), wstride, data.as_ptr(), stride, n) }
}

/// This function returns the square root of the corresponding variance function
/// [`wvariance_m`] above.
#[doc(alias = "gsl_stats_wsd_m")]
pub fn wsd_m(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize, wmean: f64) -> f64 {
    unsafe { sys::gsl_stats_wsd_m(w.as_ptr(), wstride, data.as_ptr(), stride, n, wmean) }
}

/// This function computes an unbiased estimate of the variance of the weighted dataset data when
/// the population mean mean of the underlying distribution is known a priori. In this case the
/// estimator for the variance replaces the sample mean \Hat\mu by the known population mean \mu,
///
/// \Hat\sigma^2 = (\sum w_i (x_i - \mu)^2) / (\sum w_i)
#[doc(alias = "gsl_stats_wvariance_with_fixed_mean")]
pub fn wvariance_with_fixed_mean(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    n: usize,
    mean: f64,
) -> f64 {
    unsafe {
        sys::gsl_stats_wvariance_with_fixed_mean(
            w.as_ptr(),
            wstride,
            data.as_ptr(),
            stride,
            n,
            mean,
        )
    }
}

/// The standard deviation is defined as the square root of the variance. This function returns the
/// square root of the corresponding variance function above.
#[doc(alias = "gsl_stats_wsd_with_fixed_mean")]
pub fn wsd_with_fixed_mean(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    n: usize,
    mean: f64,
) -> f64 {
    unsafe {
        sys::gsl_stats_wsd_with_fixed_mean(w.as_ptr(), wstride, data.as_ptr(), stride, n, mean)
    }
}

/// This function returns the weighted total sum of squares (TSS) of data about the weighted mean.
/// For gsl_stats_wtss_m the user-supplied value of wmean is used, and for gsl_stats_wtss it is
/// computed using gsl_stats_wmean.
///
/// TSS =  \sum w_i (x_i - wmean)^2
#[doc(alias = "gsl_stats_wtss")]
pub fn wtss(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_wtss(w.as_ptr(), wstride, data.as_ptr(), stride, n) }
}

/// This function returns the weighted total sum of squares (TSS) of data about the weighted mean.
/// For gsl_stats_wtss_m the user-supplied value of wmean is used, and for gsl_stats_wtss it is
/// computed using gsl_stats_wmean.
///
/// TSS =  \sum w_i (x_i - wmean)^2
#[doc(alias = "gsl_stats_wtss_m")]
pub fn wtss_m(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize, wmean: f64) -> f64 {
    unsafe { sys::gsl_stats_wtss_m(w.as_ptr(), wstride, data.as_ptr(), stride, n, wmean) }
}

/// This function computes the weighted absolute deviation from the weighted mean of data. The absolute deviation from the mean is defined
/// as,
///
/// absdev = (\sum w_i |x_i - \Hat\mu|) / (\sum w_i)
#[doc(alias = "gsl_stats_wabsdev")]
pub fn wabsdev(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_wabsdev(w.as_ptr(), wstride, data.as_ptr(), stride, n) }
}

/// This function computes the absolute deviation of the weighted dataset data about the given weighted mean wmean.
#[doc(alias = "gsl_stats_wabsdev_m")]
pub fn wabsdev_m(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    n: usize,
    wmean: f64,
) -> f64 {
    unsafe { sys::gsl_stats_wabsdev_m(w.as_ptr(), wstride, data.as_ptr(), stride, n, wmean) }
}

/// This function computes the weighted skewness of the dataset data.
///
/// skew = (\sum w_i ((x_i - \Hat x)/\Hat \sigma)^3) / (\sum w_i)
#[doc(alias = "gsl_stats_wskew")]
pub fn wskew(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_wskew(w.as_ptr(), wstride, data.as_ptr(), stride, n) }
}

/// This function computes the weighted skewness of the dataset data using the given values of the
/// weighted mean and weighted standard deviation, wmean and wsd.
#[doc(alias = "gsl_stats_wskew_m_sd")]
pub fn wskew_m_sd(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    n: usize,
    wmean: f64,
    wsd: f64,
) -> f64 {
    unsafe { sys::gsl_stats_wskew_m_sd(w.as_ptr(), wstride, data.as_ptr(), stride, n, wmean, wsd) }
}

/// This function computes the weighted kurtosis of the dataset data.
///
/// kurtosis = ((\sum w_i ((x_i - \Hat x)/\Hat \sigma)^4) / (\sum w_i)) - 3
#[doc(alias = "gsl_stats_wkurtosis")]
pub fn wkurtosis(w: &[f64], wstride: usize, data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_wkurtosis(w.as_ptr(), wstride, data.as_ptr(), stride, n) }
}

/// This function computes the weighted kurtosis of the dataset data using the given values of the
/// weighted mean and weighted standard deviation, wmean and wsd.
#[doc(alias = "gsl_stats_wkurtosis_m_sd")]
pub fn wkurtosis_m_sd(
    w: &[f64],
    wstride: usize,
    data: &[f64],
    stride: usize,
    n: usize,
    wmean: f64,
    wsd: f64,
) -> f64 {
    unsafe {
        sys::gsl_stats_wkurtosis_m_sd(w.as_ptr(), wstride, data.as_ptr(), stride, n, wmean, wsd)
    }
}

/// This function returns the maximum value in data, a dataset of length n with stride stride. The
/// maximum value is defined as the value of the element x_i which satisfies x_i >= x_j for all j.
///
/// If you want instead to find the element with the largest absolute magnitude you will need to
/// apply fabs or abs to your data before calling this function.
#[doc(alias = "gsl_stats_max")]
pub fn max(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_max(data.as_ptr(), stride, n) }
}

/// This function returns the minimum value in data, a dataset of length n with stride stride. The
/// minimum value is defined as the value of the element x_i which satisfies x_i <= x_j for all j.
///
/// If you want instead to find the element with the smallest absolute magnitude you will need to
/// apply fabs or abs to your data before calling this function.
#[doc(alias = "gsl_stats_min")]
pub fn min(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_min(data.as_ptr(), stride, n) }
}

/// This function finds both the minimum and maximum values min, max in data in a single pass.
///
/// Returns `(min, max)`.
#[doc(alias = "gsl_stats_minmax")]
pub fn minmax(data: &[f64], stride: usize, n: usize) -> (f64, f64) {
    let mut min = 0.;
    let mut max = 0.;
    unsafe { sys::gsl_stats_minmax(&mut min, &mut max, data.as_ptr(), stride, n) };
    (min, max)
}

/// This function returns the index of the maximum value in data, a dataset of length n with stride
/// stride. The maximum value is defined as the value of the element x_i which satisfies x_i >= x_j
/// for all j. When there are several equal maximum elements then the first one is chosen.
#[doc(alias = "gsl_stats_max_index")]
pub fn max_index(data: &[f64], stride: usize, n: usize) -> usize {
    unsafe { sys::gsl_stats_max_index(data.as_ptr(), stride, n) }
}

/// This function returns the index of the minimum value in data, a dataset of length n with stride
/// stride. The minimum value is defined as the value of the element x_i which satisfies x_i >= x_j
/// for all j. When there are several equal minimum elements then the first one is chosen.
#[doc(alias = "gsl_stats_min_index")]
pub fn min_index(data: &[f64], stride: usize, n: usize) -> usize {
    unsafe { sys::gsl_stats_min_index(data.as_ptr(), stride, n) }
}

/// This function returns the indexes min_index, max_index of the minimum and maximum values in data
/// in a single pass.
///
/// Returns `(min_index, max_index)`.
#[doc(alias = "gsl_stats_minmax_index")]
pub fn minmax_index(data: &[f64], stride: usize, n: usize) -> (usize, usize) {
    let mut min_index = 0;
    let mut max_index = 0;
    unsafe {
        sys::gsl_stats_minmax_index(&mut min_index, &mut max_index, data.as_ptr(), stride, n)
    };
    (min_index, max_index)
}

/// This function returns the median value of sorted_data, a dataset of length n with stride stride.
/// The elements of the array must be in ascending numerical order. There are no checks to see
/// whether the data are sorted, so the function gsl_sort should always be used first.
///
/// When the dataset has an odd number of elements the median is the value of element (n-1)/2. When
/// the dataset has an even number of elements the median is the mean of the two nearest middle
/// values, elements (n-1)/2 and n/2. Since the algorithm for computing the median involves
/// interpolation this function always returns a floating-point number, even for integer data types.
#[doc(alias = "gsl_stats_median_from_sorted_data")]
pub fn median_from_sorted_data(data: &[f64], stride: usize, n: usize) -> f64 {
    unsafe { sys::gsl_stats_median_from_sorted_data(data.as_ptr(), stride, n) }
}

/// This function returns a quantile value of sorted_data, a double-precision array of length n with
/// stride stride. The elements of the array must be in ascending numerical order. The quantile is
/// determined by the f, a fraction between 0 and 1. For example, to compute the value of the 75th
/// percentile f should have the value 0.75.
///
/// There are no checks to see whether the data are sorted, so the function gsl_sort should always
/// be used first.
///
/// The quantile is found by interpolation, using the formula
///
/// quantile = (1 - \delta) x_i + \delta x_{i+1}
///
/// where i is floor((n - 1)f) and \delta is (n-1)f - i.
///
/// Thus the minimum value of the array (data[0*stride]) is given by f equal to zero, the maximum
/// value (data[(n-1)*stride]) is given by f equal to one and the median value is given by f equal
/// to 0.5. Since the algorithm for computing quantiles involves interpolation this function always
/// returns a floating-point number, even for integer data types.
#[doc(alias = "gsl_stats_quantile_from_sorted_data")]
pub fn quantile_from_sorted_data(data: &[f64], stride: usize, n: usize, f: f64) -> f64 {
    unsafe { sys::gsl_stats_quantile_from_sorted_data(data.as_ptr(), stride, n, f) }
}
