//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
Linear Regression

The functions described in this section can be used to perform least-squares fits to a straight line model, Y(c,x) = c_0 + c_1 x.
!*/

use crate::{
    vector::{check_equal_len, Vector},
    Value,
};

/// This function computes the best-fit linear regression coefficients
/// (c0, c1) of the model Y = c_0 + c_1 X for the dataset (`x`, `y`),
/// two vectors of the same length (possibly with strides).
///
/// The errors on `y` are assumed unknown so the variance-covariance
/// matrix for the parameters (c0, c1) is estimated from the scatter
/// of the points around the best-fit line and returned via the
/// parameters (cov00, cov01, cov11).
///
/// The sum of squares of the residuals from the best-fit line is
/// returned in sumsq. Note: the correlation coefficient of the data
/// can be computed using gsl_stats_correlation (see
/// [`Correlation`](http://www.gnu.org/software/gsl/manual/html_node/Correlation.html#Correlation)),
/// it does not depend on the fit.
///
/// Returns `(c0, c1, cov00, cov01, cov11, sumsq)`.
///
/// # Example
///
/// ```
/// use rgsl::fit;
/// let (c0, c1, _, _, _, _) = fit::linear(&[0., 1.], &[0., 1.])?;
/// assert_eq!(c0, 0.);
/// assert_eq!(c1, 1.);
/// # Ok::<(), rgsl::Value>(())
/// ```
#[doc(alias = "gsl_fit_linear")]
pub fn linear<T>(x: &T, y: &T) -> Result<(f64, f64, f64, f64, f64, f64), Value>
where
    T: Vector<f64> + ?Sized,
{
    check_equal_len(x, y)?;
    let mut c0 = 0.;
    let mut c1 = 0.;
    let mut cov00 = 0.;
    let mut cov01 = 0.;
    let mut cov11 = 0.;
    let mut sumsq = 0.;
    let ret = unsafe {
        ::sys::gsl_fit_linear(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c0,
            &mut c1,
            &mut cov00,
            &mut cov01,
            &mut cov11,
            &mut sumsq,
        )
    };
    result_handler!(ret, (c0, c1, cov00, cov01, cov11, sumsq))
}

/// This function computes the best-fit linear regression coefficients
/// (c0, c1) of the model Y = c_0 + c_1 X for the weighted dataset
/// (`x`, `y`), two vectors of the same length (possibly with strides).
///
/// The vector `w`, of the same length as `x` and `y`, specifies the
/// weight of each datapoint.
///
/// The weight is the reciprocal of the variance for each datapoint in y.
///
/// The covariance matrix for the parameters (c0, c1) is computed using the weights and returned via
/// the parameters (cov00, cov01, cov11).
/// The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq.
///
/// Returns `(c0, c1, cov00, cov01, cov11, chisq)`.
#[doc(alias = "gsl_fit_wlinear")]
pub fn wlinear<T: Vector<f64> + ?Sized>(
    x: &T,
    w: &T,
    y: &T,
) -> Result<(f64, f64, f64, f64, f64, f64), Value> {
    check_equal_len(x, y)?;
    check_equal_len(x, w)?;
    let mut c0 = 0.;
    let mut c1 = 0.;
    let mut cov00 = 0.;
    let mut cov01 = 0.;
    let mut cov11 = 0.;
    let mut chisq = 0.;
    let ret = unsafe {
        ::sys::gsl_fit_wlinear(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(w).as_ptr(),
            T::stride(w),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c0,
            &mut c1,
            &mut cov00,
            &mut cov01,
            &mut cov11,
            &mut chisq,
        )
    };
    result_handler!(ret, (c0, c1, cov00, cov01, cov11, chisq))
}

/// This function uses the best-fit linear regression coefficients c0, c1 and their covariance
/// cov00, cov01, cov11 to compute the fitted function y and its standard deviation y_err for the
/// model Y = c_0 + c_1 X at the point x.
///
/// Returns `(y, y_err)`.
#[doc(alias = "gsl_fit_linear_est")]
pub fn linear_est(
    x: f64,
    c0: f64,
    c1: f64,
    cov00: f64,
    cov01: f64,
    cov11: f64,
) -> Result<(f64, f64), Value> {
    let mut y = 0.;
    let mut y_err = 0.;
    let ret =
        unsafe { sys::gsl_fit_linear_est(x, c0, c1, cov00, cov01, cov11, &mut y, &mut y_err) };
    result_handler!(ret, (y, y_err))
}

/// This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for
/// the datasets (x, y), two vectors of length n with strides xstride and ystride.
/// The errors on y are assumed unknown so the variance of the parameter c1 is estimated from the
/// scatter of the points around the best-fit line and returned via the parameter cov11.
/// The sum of squares of the residuals from the best-fit line is returned in sumsq.
///
/// Returns `(c1, cov11, sumsq)`.
#[doc(alias = "gsl_fit_mul")]
pub fn mul<T: Vector<f64> + ?Sized>(x: &T, y: &T) -> Result<(f64, f64, f64), Value> {
    check_equal_len(x, y)?;
    let mut c1 = 0.;
    let mut cov11 = 0.;
    let mut sumsq = 0.;
    let ret = unsafe {
        sys::gsl_fit_mul(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c1,
            &mut cov11,
            &mut sumsq,
        )
    };
    result_handler!(ret, (c1, cov11, sumsq))
}

/// Returns `(c1, cov11, sumsq)`.
#[doc(alias = "gsl_fit_wmul")]
pub fn wmul<T: Vector<f64> + ?Sized>(x: &T, w: &T, y: &T) -> Result<(f64, f64, f64), Value> {
    check_equal_len(x, y)?;
    check_equal_len(x, w)?;
    let mut c1 = 0.;
    let mut cov11 = 0.;
    let mut sumsq = 0.;
    let ret = unsafe {
        sys::gsl_fit_wmul(
            T::as_slice(x).as_ptr(),
            T::stride(x),
            T::as_slice(w).as_ptr(),
            T::stride(w),
            T::as_slice(y).as_ptr(),
            T::stride(y),
            T::len(x),
            &mut c1,
            &mut cov11,
            &mut sumsq,
        )
    };
    result_handler!(ret, (c1, cov11, sumsq))
}

/// This function uses the best-fit linear regression coefficient c1 and its covariance cov11 to
/// compute the fitted function y and its standard deviation y_err for the model Y = c_1 X at the
/// point x.
///
/// Returns `(y, y_err)`.
#[doc(alias = "gsl_fit_mul_est")]
pub fn mul_est(x: f64, c1: f64, cov11: f64) -> Result<(f64, f64), Value> {
    let mut y = 0.;
    let mut y_err = 0.;
    let ret = unsafe { sys::gsl_fit_mul_est(x, c1, cov11, &mut y, &mut y_err) };
    result_handler!(ret, (y, y_err))
}
