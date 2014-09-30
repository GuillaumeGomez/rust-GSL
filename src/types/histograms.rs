//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Histograms

This chapter describes functions for creating histograms. Histograms provide a convenient way of summarizing the distribution of a set of data. 
A histogram consists of a set of bins which count the number of events falling into a given range of a continuous variable x. In GSL the bins 
of a histogram contain floating-point numbers, so they can be used to record both integer and non-integer distributions. The bins can use arbitrary 
sets of ranges (uniformly spaced bins are the default). Both one and two-dimensional histograms are supported.

Once a histogram has been created it can also be converted into a probability distribution function. The library provides efficient routines 
for selecting random samples from probability distributions. This can be useful for generating simulations based on real data.

##Resampling from histograms

A histogram made by counting events can be regarded as a measurement of a probability distribution. Allowing for statistical error, the height 
of each bin represents the probability of an event where the value of x falls in the range of that bin. The probability distribution function 
has the one-dimensional form p(x)dx where,

p(x) = n_i/ (N w_i)
In this equation n_i is the number of events in the bin which contains x, w_i is the width of the bin and N is the total number of events. 
The distribution of events within each bin is assumed to be uniform.
!*/

use ffi;
use enums;

pub struct Histogram {
    h: *mut ffi::gsl_histogram
}

impl Histogram {
    /// This function allocates memory for a histogram with n bins, and returns a pointer to a newly created gsl_histogram struct. If insufficient
    /// memory is available a null pointer is returned and the error handler is invoked with an error code of Value::NoMem. The bins and ranges are
    /// not initialized, and should be prepared using one of the range-setting functions below in order to make the histogram ready for use.
    pub fn new(n: u64) -> Option<Histogram> {
        let tmp = unsafe { ffi::gsl_histogram_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Histogram {
                h: tmp
            })
        }
    }

    /// This function sets the ranges of the existing histogram h using the array range of size size. The values of the histogram bins are reset
    /// to zero. The range array should contain the desired bin limits. The ranges can be arbitrary, subject to the restriction that they are
    /// monotonically increasing.
    /// 
    /// The following example shows how to create a histogram with logarithmic bins with ranges [1,10), [10,100) and [100,1000).
    ///
    /// ```C
    /// gsl_histogram * h = gsl_histogram_alloc (3);
    /// 
    /// /* bin[0] covers the range 1 <= x < 10 */
    /// /* bin[1] covers the range 10 <= x < 100 */
    /// /* bin[2] covers the range 100 <= x < 1000 */
    /// 
    /// double range[4] = { 1.0, 10.0, 100.0, 1000.0 };
    /// 
    /// gsl_histogram_set_ranges (h, range, 4);
    /// ```
    /// 
    /// Note that the size of the range array should be defined to be one element bigger than the number of bins. The additional element is
    /// required for the upper value of the final bin.
    pub fn set_ranges(&self, range: &[f64]) -> enums::Value {
        unsafe { ffi::gsl_histogram_set_ranges(self.h, range.as_ptr(), range.len() as u64) }
    }

    /// This function sets the ranges of the existing histogram h to cover the range xmin to xmax uniformly. The values of the histogram bins
    /// are reset to zero. The bin ranges are shown in the table below,
    /// 
    /// bin[0] corresponds to xmin <= x < xmin + d
    /// bin[1] corresponds to xmin + d <= x < xmin + 2 d
    /// ......
    /// bin[n-1] corresponds to xmin + (n-1)d <= x < xmax
    /// where d is the bin spacing, d = (xmax-xmin)/n.
    pub fn set_ranges_uniform(&self, xmin: f64, xmax: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram_set_ranges_uniform(self.h, xmin, xmax) }
    }

    /// This function copies the self histogram into the pre-existing histogram dest, making dest into an exact copy of self. The two histograms
    /// must be of the same size.
    pub fn copy(&self, dest: &Histogram) -> enums::Value {
        unsafe { ffi::gsl_histogram_memcpy(dest.h, self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns a pointer to a newly created histogram which is an exact copy of the self histogram.
    pub fn clone(&self) -> Option<Histogram> {
        let tmp = unsafe { ffi::gsl_histogram_clone(self.h as *const ffi::gsl_histogram) };

        if tmp.is_null() {
            None
        } else {
            Some(Histogram {
                h: tmp
            })
        }
    }

    /// This function updates the self histogram by adding one (1.0) to the bin whose range contains the coordinate x.
    /// 
    /// If x lies in the valid range of the histogram then the function returns zero to indicate success. If x is less than the lower limit of
    /// the histogram then the function returns Value::Dom, and none of bins are modified. Similarly, if the value of x is greater than or equal
    /// to the upper limit of the histogram then the function returns Value::Dom, and none of the bins are modified. The error handler is not
    /// called, however, since it is often necessary to compute histograms for a small range of a larger dataset, ignoring the values outside
    /// the range of interest.
    pub fn increment(&self, x: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram_increment(self.h, x) }
    }

    /// This function is similar to gsl_histogram_increment but increases the value of the appropriate bin in the histogram h by the floating-point
    /// number weight.
    pub fn accumulate(&self, x: f64, weight: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram_accumulate(self.h, x, weight) }
    }

    /// This function returns the contents of the i-th bin of the histogram h. If i lies outside the valid range of indices for the histogram then
    /// the error handler is called with an error code of Value::Dom and the function returns 0.
    pub fn get(&self, i: u64) -> f64 {
        unsafe { ffi::gsl_histogram_get(self.h as *const ffi::gsl_histogram, i) }
    }

    /// This function finds the upper and lower range limits of the i-th bin of the self histogram. If the index i is valid then the corresponding
    /// range limits are stored in lower and upper. The lower limit is inclusive (i.e. events with this coordinate are included in the bin) and
    /// the upper limit is exclusive (i.e. events with the coordinate of the upper limit are excluded and fall in the neighboring higher bin,
    /// if it exists). The function returns 0 to indicate success. If i lies outside the valid range of indices for the histogram then
    /// the error handler is called and the function returns an error code of Value::Dom.
    pub fn get_range(&self, i: u64, lower: &mut f64, upper: &mut f64) -> enums::Value {
        unsafe { ffi::gsl_histogram_get_range(self.h as *const ffi::gsl_histogram, i, lower, upper) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins of the self histogram. They provide a way
    /// of determining these values without accessing the gsl_histogram struct directly.
    pub fn max(&self) -> f64 {
        unsafe { ffi::gsl_histogram_max(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins of the self histogram. They provide a way
    /// of determining these values without accessing the gsl_histogram struct directly.
    pub fn min(&self) -> f64 {
        unsafe { ffi::gsl_histogram_min(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins of the self histogram. They provide a way
    /// of determining these values without accessing the gsl_histogram struct directly.
    pub fn bins(&self) -> u64 {
        unsafe { ffi::gsl_histogram_bins(self.h as *const ffi::gsl_histogram) }
    }

    /// This function resets all the bins in the self histogram to zero.
    pub fn reset(&self) {
        unsafe { ffi::gsl_histogram_reset(self.h) }
    }

    /// This function finds and sets the index i to the bin number which covers the coordinate x in the self histogram. The bin is located using
    /// a binary search. The search includes an optimization for histograms with uniform range, and will return the correct bin immediately in
    /// this case. If x is found in the range of the histogram then the function sets the index i and returns Value::Success. If x lies outside
    /// the valid range of the histogram then the function returns Value::Dom and the error handler is invoked.
    pub fn find(&self, x: f64, i: &mut u64) -> enums::Value {
        unsafe { ffi::gsl_histogram_find(self.h as *const ffi::gsl_histogram, x, i) }
    }

    /// This function returns the maximum value contained in the histogram bins.
    pub fn max_val(&self) -> f64 {
        unsafe { ffi::gsl_histogram_max_val(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the index of the bin containing the maximum value. In the case where several bins contain the same maximum value
    /// the smallest index is returned.
    pub fn max_bin(&self) -> u64 {
        unsafe { ffi::gsl_histogram_max_bin(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the minimum value contained in the histogram bins.
    pub fn min_val(&self) -> f64 {
        unsafe { ffi::gsl_histogram_min_val(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the index of the bin containing the minimum value. In the case where several bins contain the same maximum value
    /// the smallest index is returned.
    pub fn min_bin(&self) -> u64 {
        unsafe { ffi::gsl_histogram_min_bin(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the mean of the histogrammed variable, where the histogram is regarded as a probability distribution. Negative
    /// bin values are ignored for the purposes of this calculation. The accuracy of the result is limited by the bin width.
    pub fn mean(&self) -> f64 {
        unsafe { ffi::gsl_histogram_mean(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the standard deviation of the histogrammed variable, where the histogram is regarded as a probability distribution.
    /// Negative bin values are ignored for the purposes of this calculation. The accuracy of the result is limited by the bin width.
    pub fn sigma(&self) -> f64 {
        unsafe { ffi::gsl_histogram_sigma(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns the sum of all bin values. Negative bin values are included in the sum.
    pub fn sum(&self) -> f64 {
        unsafe { ffi::gsl_histogram_sum(self.h as *const ffi::gsl_histogram) }
    }

    /// This function returns true if the all of the individual bin ranges of the two histograms are identical, and false otherwise.
    pub fn equal_bins_p(&self, other: &Histogram) -> bool {
        match unsafe { ffi::gsl_histogram_equal_bins_p(self.h as *const ffi::gsl_histogram, other.h as *const ffi::gsl_histogram) } {
            0i32 => false,
            _ => true
        }
    }

    /// This function adds the contents of the bins in histogram other to the corresponding bins of self histogram, i.e. h'_1(i) = h_1(i) + h_2(i).
    /// The two histograms must have identical bin ranges.
    pub fn add(&self, other: &Histogram) -> enums::Value {
        unsafe { ffi::gsl_histogram_add(self.h, other.h as *const ffi::gsl_histogram) }
    }

    /// This function subtracts the contents of the bins in histogram other from the corresponding bins of self histogram, i.e. h'_1(i) = h_1(i) - h_2(i).
    /// The two histograms must have identical bin ranges.
    pub fn sub(&self, other: &Histogram) -> enums::Value {
        unsafe { ffi::gsl_histogram_sub(self.h, other.h as *const ffi::gsl_histogram) }
    }

    /// This function multiplies the contents of the bins of self histogram by the contents of the corresponding bins in other histogram, i.e. h'_1(i) =
    /// h_1(i) * h_2(i). The two histograms must have identical bin ranges.
    pub fn mul(&self, other: &Histogram) -> enums::Value {
        unsafe { ffi::gsl_histogram_mul(self.h, other.h as *const ffi::gsl_histogram) }
    }

    /// This function divides the contents of the bins of self histogram by the contents of the corresponding bins in other histogram, i.e. h'_1(i) = h_1(i)
    /// / h_2(i). The two histograms must have identical bin ranges.
    pub fn div(&self, other: &Histogram) -> enums::Value {
        unsafe { ffi::gsl_histogram_div(self.h, other.h as *const ffi::gsl_histogram) }
    }

    /// This function multiplies the contents of the bins of self histogram by the constant scale, i.e. h'_1(i) = h_1(i) * scale.
    pub fn scale(&self, scale: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram_scale(self.h, scale) }
    }

    /// This function shifts the contents of the bins of self histogram by the constant offset, i.e. h'_1(i) = h_1(i) + offset.
    pub fn shift(&self, offset: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram_shift(self.h, offset) }
    }

    #[allow(unused_must_use)]
    pub fn print(&self, stream: &mut Writer/*, range_format: &str, bin_format: &str*/) -> enums::Value {
        unsafe {
            let n = (*self.h).n as int;

            for i in range(0i, n) {
                write!(stream, "{}", *(*self.h).range.offset(i));
                write!(stream, " ");
                write!(stream, "{}", *(*self.h).range.offset(i + 1));
                write!(stream, " ");
                write!(stream, "{}", *(*self.h).range.offset(i));
                write!(stream, "\n");
            }

            enums::Success
        }
    }
}

impl Drop for Histogram {
    fn drop(&mut self) {
        unsafe { ffi::gsl_histogram_free(self.h) };
        self.h = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_histogram> for Histogram {
    fn wrap(h: *mut ffi::gsl_histogram) -> Histogram {
        Histogram {
            h: h
        }
    }

    fn unwrap(h: &Histogram) -> *mut ffi::gsl_histogram {
        h.h
    }
}

/// The probability distribution function for a histogram consists of a set of bins which measure the probability of an event falling into a
/// given range of a continuous variable x. A probability distribution function is defined by the following struct, which actually stores the
/// cumulative probability distribution function. This is the natural quantity for generating samples via the inverse transform method, because
/// there is a one-to-one mapping between the cumulative probability distribution and the range [0,1]. It can be shown that by taking a uniform
/// random number in this range and finding its corresponding coordinate in the cumulative probability distribution we obtain samples with the
/// desired probability distribution.
pub struct HistogramPdf {
    h: *mut ffi::gsl_histogram_pdf
}

impl HistogramPdf {
    /// This function allocates memory for a probability distribution with n bins and returns a pointer to a newly initialized gsl_histogram_pdf
    /// struct. If insufficient memory is available a null pointer is returned and the error handler is invoked with an error code of Value::NoMem.
    pub fn new(n: u64) -> Option<HistogramPdf> {
        let tmp = unsafe { ffi::gsl_histogram_pdf_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(HistogramPdf {
                h: tmp
            })
        }
    }

    /// This function initializes the probability distribution self with the contents of the histogram h. If any of the bins of h are negative then
    /// the error handler is invoked with an error code of Value::Dom because a probability distribution cannot contain negative values.
    pub fn init(&self, h: &Histogram) -> enums::Value {
        unsafe { ffi::gsl_histogram_pdf_init(self.h, h.h as *const ffi::gsl_histogram) }
    }

    /// This function uses r, a uniform random number between zero and one, to compute a single random sample from the probability distribution
    /// self. The algorithm used to compute the sample s is given by the following formula,
    /// 
    /// s = range[i] + delta * (range[i+1] - range[i])
    /// 
    /// where i is the index which satisfies sum[i] <= r < sum[i+1] and delta is (r - sum[i])/(sum[i+1] - sum[i]).
    pub fn sample(&self, r: f64) -> f64 {
        unsafe { ffi::gsl_histogram_pdf_sample(self.h as *const ffi::gsl_histogram_pdf, r) }
    }
}

impl Drop for HistogramPdf {
    fn drop(&mut self) {
        unsafe { ffi::gsl_histogram_pdf_free(self.h) };
        self.h = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_histogram_pdf> for HistogramPdf {
    fn wrap(h: *mut ffi::gsl_histogram_pdf) -> HistogramPdf {
        HistogramPdf {
            h: h
        }
    }

    fn unwrap(h: &HistogramPdf) -> *mut ffi::gsl_histogram_pdf {
        h.h
    }
}

/// A two dimensional histogram consists of a set of bins which count the number of events falling in a given area of the (x,y) plane. The simplest
/// way to use a two dimensional histogram is to record two-dimensional position information, n(x,y). Another possibility is to form a joint
/// distribution by recording related variables. For example a detector might record both the position of an event (x) and the amount of energy
/// it deposited E. These could be histogrammed as the joint distribution n(x,E).
pub struct Histogram2D {
    h: *mut ffi::gsl_histogram2d
}

impl Histogram2D {
    /// This function allocates memory for a two-dimensional histogram with nx bins in the x direction and ny bins in the y direction. The
    /// function returns a pointer to a newly created gsl_histogram2d struct. If insufficient memory is available a null pointer is returned
    /// and the error handler is invoked with an error code of Value::NoMem. The bins and ranges must be initialized with one of the
    /// functions below before the histogram is ready for use.
    pub fn new(nx: u64, ny: u64) -> Option<Histogram2D> {
        let tmp = unsafe { ffi::gsl_histogram2d_alloc(nx, ny) };

        if tmp.is_null() {
            None
        } else {
            Some(Histogram2D {
                h: tmp
            })
        }
    }

    /// This function sets the ranges of the existing histogram h using the arrays xrange and yrange of size xsize and ysize respectively.
    /// The values of the histogram bins are reset to zero.
    pub fn set_ranges(&self, xrange: &[f64], yrange: &[f64]) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_set_ranges(self.h, xrange.as_ptr(), xrange.len() as u64, yrange.as_ptr(), yrange.len() as u64) }
    }

    /// This function sets the ranges of the existing histogram h to cover the ranges xmin to xmax and ymin to ymax uniformly. The values
    /// of the histogram bins are reset to zero.
    pub fn set_ranges_uniform(&self, xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_set_ranges_uniform(self.h, xmin, xmax, ymin, ymax) }
    }

    /// This function copies the histogram src into the pre-existing histogram dest, making dest into an exact copy of src. The two histograms
    /// must be of the same size.
    pub fn copy(&self, dest: &Histogram2D) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_memcpy(dest.h, self.h as *const ffi::gsl_histogram2d) }
    }

    /// his function returns a pointer to a newly created histogram which is an exact copy of the histogram self.
    pub fn clone(&self) -> Option<Histogram2D> {
        let tmp = unsafe { ffi::gsl_histogram2d_clone(self.h as *const ffi::gsl_histogram2d) };

        if tmp.is_null() {
            None
        } else {
            Some(Histogram2D {
                h: tmp
            })
        }
    }

    /// This function updates the histogram h by adding one (1.0) to the bin whose x and y ranges contain the coordinates (x,y).
    /// 
    /// If the point (x,y) lies inside the valid ranges of the histogram then the function returns zero to indicate success. If (x,y) lies
    /// outside the limits of the histogram then the function returns Value::Dom, and none of the bins are modified. The error handler is not
    /// called, since it is often necessary to compute histograms for a small range of a larger dataset, ignoring any coordinates outside the
    /// range of interest.
    pub fn increment(&self, x: f64, y: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_increment(self.h, x, y) }
    }

    /// This function is similar to gsl_histogram2d_increment but increases the value of the appropriate bin in the histogram h by the floating-point
    /// number weight.
    pub fn accumulate(&self, x: f64, y: f64, weight: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_accumulate(self.h, x, y, weight) }
    }

    /// This function returns the contents of the (i,j)-th bin of the histogram h. If (i,j) lies outside the valid range of indices for the
    /// histogram then the error handler is called with an error code of Value::Dom and the function returns 0.
    pub fn get(&self, i: u64, j: u64) -> f64 {
        unsafe { ffi::gsl_histogram2d_get(self.h as *const ffi::gsl_histogram2d, i, j) }
    }

    /// This function finds the upper and lower range limits of the i-th and j-th bins in the x and y directions of the histogram h. The range
    /// limits are stored in xlower and xupper or ylower and yupper. The lower limits are inclusive (i.e. events with these coordinates are included
    /// in the bin) and the upper limits are exclusive (i.e. events with the value of the upper limit are not included and fall in the neighboring
    /// higher bin, if it exists). The functions return 0 to indicate success. If i or j lies outside the valid range of indices for the histogram
    /// then the error handler is called with an error code of Value::Dom.
    pub fn get_xrange(&self, i: u64, xlower: &mut f64, xupper: &mut f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_get_xrange(self.h as *const ffi::gsl_histogram2d, i, xlower, xupper) }
    }

    /// This function finds the upper and lower range limits of the i-th and j-th bins in the x and y directions of the histogram h. The range
    /// limits are stored in xlower and xupper or ylower and yupper. The lower limits are inclusive (i.e. events with these coordinates are included
    /// in the bin) and the upper limits are exclusive (i.e. events with the value of the upper limit are not included and fall in the neighboring
    /// higher bin, if it exists). The functions return 0 to indicate success. If i or j lies outside the valid range of indices for the histogram
    /// then the error handler is called with an error code of Value::Dom.
    pub fn get_yrange(&self, j: u64, ylower: &mut f64, yupper: &mut f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_get_yrange(self.h as *const ffi::gsl_histogram2d, j, ylower, yupper) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins for the x and y directions of the histogram h.
    /// They provide a way of determining these values without accessing the gsl_histogram2d struct directly.
    pub fn xmax(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_xmax(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins for the x and y directions of the histogram h.
    /// They provide a way of determining these values without accessing the gsl_histogram2d struct directly.
    pub fn xmin(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_xmin(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins for the x and y directions of the histogram h.
    /// They provide a way of determining these values without accessing the gsl_histogram2d struct directly.
    pub fn nx(&self) -> u64 {
        unsafe { ffi::gsl_histogram2d_nx(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins for the x and y directions of the histogram h.
    /// They provide a way of determining these values without accessing the gsl_histogram2d struct directly.
    pub fn ymax(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_ymax(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins for the x and y directions of the histogram h.
    /// They provide a way of determining these values without accessing the gsl_histogram2d struct directly.
    pub fn ymin(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_ymin(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the maximum upper and minimum lower range limits and the number of bins for the x and y directions of the histogram h.
    /// They provide a way of determining these values without accessing the gsl_histogram2d struct directly.
    pub fn ny(&self) -> u64 {
        unsafe { ffi::gsl_histogram2d_ny(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function resets all the bins of the histogram h to zero.
    pub fn reset(&self) {
        unsafe { ffi::gsl_histogram2d_reset(self.h) }
    }

    /// This function finds and sets the indices i and j to the bin which covers the coordinates (x,y). The bin is located using a binary search.
    /// The search includes an optimization for histograms with uniform ranges, and will return the correct bin immediately in this case. If
    /// (x,y) is found then the function sets the indices (i,j) and returns Value::Success. If (x,y) lies outside the valid range of the histogram
    /// then the function returns Value::Dom and the error handler is invoked.
    pub fn find(&self, x: f64, y: f64, i: &mut u64, j: &mut u64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_find(self.h as *const ffi::gsl_histogram2d, x, y, i, j) }
    }

    /// This function returns the maximum value contained in the histogram bins.
    pub fn max_val(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_max_val(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function finds the indices of the bin containing the maximum value in the histogram h and stores the result in (i,j). In the case
    /// where several bins contain the same maximum value the first bin found is returned.
    pub fn max_bin(&self, i: &mut u64, j: &mut u64) {
        unsafe { ffi::gsl_histogram2d_max_bin(self.h as *const ffi::gsl_histogram2d, i, j) }
    }

    /// This function returns the minimum value contained in the histogram bins.
    pub fn min_val(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_min_val(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function finds the indices of the bin containing the minimum value in the histogram h and stores the result in (i,j). In the case
    /// where several bins contain the same maximum value the first bin found is returned.
    pub fn min_bin(&self, i: &mut u64, j: &mut u64) {
        unsafe { ffi::gsl_histogram2d_min_bin(self.h as *const ffi::gsl_histogram2d, i, j) }
    }

    /// This function returns the mean of the histogrammed x variable, where the histogram is regarded as a probability distribution. Negative
    /// bin values are ignored for the purposes of this calculation.
    pub fn xmean(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_xmean(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the mean of the histogrammed y variable, where the histogram is regarded as a probability distribution. Negative
    /// bin values are ignored for the purposes of this calculation.
    pub fn ymean(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_ymean(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the standard deviation of the histogrammed x variable, where the histogram is regarded as a probability
    /// distribution. Negative bin values are ignored for the purposes of this calculation.
    pub fn xsigma(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_xsigma(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the standard deviation of the histogrammed y variable, where the histogram is regarded as a probability
    /// distribution. Negative bin values are ignored for the purposes of this calculation.
    pub fn ysigma(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_ysigma(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the covariance of the histogrammed x and y variables, where the histogram is regarded as a probability
    /// distribution. Negative bin values are ignored for the purposes of this calculation.
    pub fn cov(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_cov(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns the sum of all bin values. Negative bin values are included in the sum.
    pub fn sum(&self) -> f64 {
        unsafe { ffi::gsl_histogram2d_sum(self.h as *const ffi::gsl_histogram2d) }
    }

    /// This function returns 1 if all the individual bin ranges of the two histograms are identical, and 0 otherwise.
    pub fn equal_bins_p(&self, other: &Histogram2D) -> bool {
        match unsafe { ffi::gsl_histogram2d_equal_bins_p(self.h as *const ffi::gsl_histogram2d, other.h as *const ffi::gsl_histogram2d) } {
            0 => false,
            _ => true
        }
    }

    /// This function adds the contents of the bins in histogram h2 to the corresponding bins of histogram h1, i.e. h'_1(i,j) = h_1(i,j)
    /// + h_2(i,j). The two histograms must have identical bin ranges.
    pub fn add(&self, other: &Histogram2D) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_add(self.h, other.h as *const ffi::gsl_histogram2d) }
    }

    /// This function subtracts the contents of the bins in histogram h2 from the corresponding bins of histogram h1, i.e. h'_1(i,j) = h_1(i,j)
    /// - h_2(i,j). The two histograms must have identical bin ranges.
    pub fn sub(&self, other: &Histogram2D) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_sub(self.h, other.h as *const ffi::gsl_histogram2d) }
    }

    /// This function multiplies the contents of the bins of histogram h1 by the contents of the corresponding bins in histogram h2, i.e. h'_1(i,j)
    /// = h_1(i,j) * h_2(i,j). The two histograms must have identical bin ranges.
    pub fn mul(&self, other: &Histogram2D) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_mul(self.h, other.h as *const ffi::gsl_histogram2d) }
    }

    /// This function divides the contents of the bins of histogram h1 by the contents of the corresponding bins in histogram h2, i.e. h'_1(i,j) =
    /// h_1(i,j) / h_2(i,j). The two histograms must have identical bin ranges.
    pub fn div(&self, other: &Histogram2D) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_div(self.h, other.h as *const ffi::gsl_histogram2d) }
    }

    /// This function multiplies the contents of the bins of histogram h by the constant scale, i.e. h'_1(i,j) = h_1(i,j) scale.
    pub fn scale(&self, scale: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_scale(self.h, scale) }
    }

    /// This function shifts the contents of the bins of histogram h by the constant offset, i.e. h'_1(i,j) = h_1(i,j) + offset.
    pub fn shift(&self, offset: f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_shift(self.h, offset) }
    }
}

impl Drop for Histogram2D {
    fn drop(&mut self) {
        unsafe { ffi::gsl_histogram2d_free(self.h) };
        self.h = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_histogram2d> for Histogram2D {
    fn wrap(h: *mut ffi::gsl_histogram2d) -> Histogram2D {
        Histogram2D {
            h: h
        }
    }

    fn unwrap(h: &Histogram2D) -> *mut ffi::gsl_histogram2d {
        h.h
    }
}



/// As in the one-dimensional case, a two-dimensional histogram made by counting events can be regarded as a measurement of a probability distribution.
/// Allowing for statistical error, the height of each bin represents the probability of an event where (x,y) falls in the range of that bin. For a
/// two-dimensional histogram the probability distribution takes the form p(x,y) dx dy where,
/// 
/// p(x,y) = n_{ij}/ (N A_{ij})
/// 
/// In this equation n_{ij} is the number of events in the bin which contains (x,y), A_{ij} is the area of the bin and N is the total number of
/// events. The distribution of events within each bin is assumed to be uniform.
pub struct Histogram2DPdf {
    h: *mut ffi::gsl_histogram2d_pdf
}

impl Histogram2DPdf {
    /// This function allocates memory for a two-dimensional probability distribution of size nx-by-ny and returns a pointer to a newly initialized
    /// gsl_histogram2d_pdf struct. If insufficient memory is available a null pointer is returned and the error handler is invoked with an error
    /// code of Value::NoMem.
    pub fn new(nx: u64, ny: u64) -> Option<Histogram2DPdf> {
        let tmp = unsafe { ffi::gsl_histogram2d_pdf_alloc(nx, ny) };

        if tmp.is_null() {
            None
        } else {
            Some(Histogram2DPdf {
                h: tmp
            })
        }
    }

    ///This function initializes the two-dimensional probability distribution calculated p from the histogram h. If any of the bins of h are
    /// negative then the error handler is invoked with an error code of GSL_EDOM because a probability distribution cannot contain negative
    /// values.
    pub fn init(&self, h: &Histogram2D) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_pdf_init(self.h, h.h as *const ffi::gsl_histogram2d) }
    }

    /// This function uses two uniform random numbers between zero and one, r1 and r2, to compute a single random sample from the two-dimensional
    /// probability distribution p.
    pub fn sample(&self, r1: f64, r2: f64, x: &mut f64, y: &mut f64) -> enums::Value {
        unsafe { ffi::gsl_histogram2d_pdf_sample(self.h as *const ffi::gsl_histogram2d_pdf, r1, r2, x, y) }
    }
}

impl Drop for Histogram2DPdf {
    fn drop(&mut self) {
        unsafe { ffi::gsl_histogram2d_pdf_free(self.h) };
        self.h = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_histogram2d_pdf> for Histogram2DPdf {
    fn wrap(h: *mut ffi::gsl_histogram2d_pdf) -> Histogram2DPdf {
        Histogram2DPdf {
            h: h
        }
    }

    fn unwrap(h: &Histogram2DPdf) -> *mut ffi::gsl_histogram2d_pdf {
        h.h
    }
}