//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
Linear Regression

The functions described in this section can be used to perform least-squares fits to a straight line model, Y(c,x) = c_0 + c_1 x.
!*/

use enums;

/// This function computes the best-fit linear regression coefficients (c0,c1) of the model Y = c_0 + c_1 X for the dataset (x, y), two vectors of length n with strides xstride and ystride.
/// The errors on y are assumed unknown so the variance-covariance matrix for the parameters (c0, c1) is estimated from the scatter of the points around the best-fit line and returned via the parameters (cov00, cov01, cov11).
/// The sum of squares of the residuals from the best-fit line is returned in sumsq. Note: the correlation coefficient of the data can be computed using gsl_stats_correlation (see [`Correlation`](http://www.gnu.org/software/gsl/manual/html_node/Correlation.html#Correlation)), it does not depend on the fit.
pub fn linear(
    x: &[f64],
    xstride: usize,
    y: &[f64],
    ystride: usize,
    n: usize,
    c0: &mut f64,
    c1: &mut f64,
    cov00: &mut f64,
    cov01: &mut f64,
    cov11: &mut f64,
    sumsq: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        ::ffi::gsl_fit_linear(
            x.as_ptr(),
            xstride,
            y.as_ptr(),
            ystride,
            n,
            c0,
            c1,
            cov00,
            cov01,
            cov11,
            sumsq,
        )
    })
}

/// This function computes the best-fit linear regression coefficients (c0,c1) of the model Y = c_0 + c_1 X for the weighted dataset (x, y), two vectors of length n with strides xstride and ystride.
/// The vector w, of length n and stride wstride, specifies the weight of each datapoint.
/// The weight is the reciprocal of the variance for each datapoint in y.
///
/// The covariance matrix for the parameters (c0, c1) is computed using the weights and returned via the parameters (cov00, cov01, cov11).
/// The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq.
pub fn wlinear(
    x: &[f64],
    xstride: usize,
    w: &[f64],
    wstride: usize,
    y: &[f64],
    ystride: usize,
    n: usize,
    c0: &mut f64,
    c1: &mut f64,
    cov00: &mut f64,
    cov01: &mut f64,
    cov11: &mut f64,
    chisq: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        ::ffi::gsl_fit_wlinear(
            x.as_ptr(),
            xstride,
            w.as_ptr(),
            wstride,
            y.as_ptr(),
            ystride,
            n,
            c0,
            c1,
            cov00,
            cov01,
            cov11,
            chisq,
        )
    })
}

/// This function uses the best-fit linear regression coefficients c0, c1 and their covariance cov00, cov01, cov11 to compute the fitted function y and its standard deviation y_err for the model Y = c_0 + c_1 X at the point x.
pub fn linear_est(
    x: f64,
    c0: f64,
    c1: f64,
    cov00: f64,
    cov01: f64,
    cov11: f64,
    y: &mut f64,
    y_err: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        ::ffi::gsl_fit_linear_est(x, c0, c1, cov00, cov01, cov11, y, y_err)
    })
}

/// This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for the datasets (x, y), two vectors of length n with strides xstride and ystride.
/// The errors on y are assumed unknown so the variance of the parameter c1 is estimated from the scatter of the points around the best-fit line and returned via the parameter cov11.
/// The sum of squares of the residuals from the best-fit line is returned in sumsq.
pub fn mul(
    x: &[f64],
    xstride: usize,
    y: &[f64],
    ystride: usize,
    n: usize,
    c1: &mut f64,
    cov11: &mut f64,
    sumsq: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        ::ffi::gsl_fit_mul(
            x.as_ptr(),
            xstride,
            y.as_ptr(),
            ystride,
            n,
            c1,
            cov11,
            sumsq,
        )
    })
}

/// This function computes the best-fit linear regression coefficient c1 of the model Y = c_1 X for the weighted datasets (x, y), two vectors of length n with strides xstride and ystride.
/// The vector w, of length n and stride wstride, specifies the weight of each datapoint.
/// The weight is the reciprocal of the variance for each datapoint in y.
///
/// The variance of the parameter c1 is computed using the weights and returned via the parameter cov11.
/// The weighted sum of squares of the residuals from the best-fit line, \chi^2, is returned in chisq.
pub fn wmul(
    x: &[f64],
    xstride: usize,
    w: &[f64],
    wstride: usize,
    y: &[f64],
    ystride: usize,
    n: usize,
    c1: &mut f64,
    cov11: &mut f64,
    sumsq: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        ::ffi::gsl_fit_wmul(
            x.as_ptr(),
            xstride,
            w.as_ptr(),
            wstride,
            y.as_ptr(),
            ystride,
            n,
            c1,
            cov11,
            sumsq,
        )
    })
}

/// This function uses the best-fit linear regression coefficient c1 and its covariance cov11 to compute the fitted function y and its standard deviation y_err for the model Y = c_1 X at the point x.
pub fn mul_est(x: f64, c1: f64, cov11: f64, y: &mut f64, y_err: &mut f64) -> enums::Value {
    enums::Value::from(unsafe { ::ffi::gsl_fit_mul_est(x, c1, cov11, y, y_err) })
}
