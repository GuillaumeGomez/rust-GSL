//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Numerical Differentiation

The functions described in this chapter compute numerical derivatives by finite differencing. An adaptive algorithm is used to find the
best choice of finite difference and to estimate the error in the derivative.

##References and Further Reading

The algorithms used by these functions are described in the following sources:

Abramowitz and Stegun, Handbook of Mathematical Functions, Section 25.3.4, and Table 25.5 (Coefficients for Differentiation).
S.D. Conte and Carl de Boor, Elementary Numerical Analysis: An Algorithmic Approach, McGraw-Hill, 1972.
!*/

use enums;

fn my_max(x1: f64, x2: f64) -> f64 {
    if x1 > x2 {
        x1
    } else {
        x2
    }
}

fn central_deriv<T>(
    f: ::function<T>,
    param: &mut T,
    x: f64,
    h: f64,
    result: &mut f64,
    abs_err_round: &mut f64,
    abs_err_trunc: &mut f64,
) {
    let fm1 = f(x - h, param);
    let fp1 = f(x + h, param);

    let fmh = f(x - h / 2f64, param);
    let fph = f(x + h / 2f64, param);

    let r3 = 0.5f64 * (fp1 - fm1);
    let r5 = (4f64 / 3f64) * (fph - fmh) - (1f64 / 3f64) * r3;

    let e3 = unsafe { (fp1.abs() + fm1.abs()) * ::DBL_EPSILON };
    let e5 = unsafe { 2f64 * (fph.abs() + fmh.abs()) * ::DBL_EPSILON + e3 };

    let dy = unsafe { my_max((r3 / h).abs(), (r5 / h).abs()) * (x.abs() / h) * ::DBL_EPSILON };

    *result = r5 / h;
    *abs_err_trunc = unsafe { ((r5 - r3) / h).abs() };
    *abs_err_round = unsafe { (e5 / h).abs() + dy };
}

/// This function computes the numerical derivative of the function f at the point x using an adaptive central difference algorithm with a step-size
/// of h. The derivative is returned in result and an estimate of its absolute error is returned in abserr.
///
/// The initial value of h is used to estimate an optimal step-size, based on the scaling of the truncation error and round-off error in the
/// derivative calculation. The derivative is computed using a 5-point rule for equally spaced abscissae at x-h, x-h/2, x, x+h/2, x+h, with
/// an error estimate taken from the difference between the 5-point rule and the corresponding 3-point rule x-h, x, x+h. Note that the value
/// of the function at x does not contribute to the derivative calculation, so only 4-points are actually used.
pub fn deriv_central<T>(
    f: ::function<T>,
    param: &mut T,
    x: f64,
    h: f64,
    result: &mut f64,
    abs_err: &mut f64,
) -> enums::Value {
    let mut r_0 = 0f64;
    let mut round = 0f64;
    let mut trunc = 0f64;

    central_deriv(f, param, x, h, &mut r_0, &mut round, &mut trunc);
    let mut error = round + trunc;

    if round < trunc && (round > 0f64 && trunc > 0f64) {
        let mut r_opt = 0f64;
        let mut round_opt = 0f64;
        let mut trunc_opt = 0f64;

        let h_opt = unsafe { h * (round / (2f64 * trunc)).powf(1f64 / 3f64) };
        central_deriv(
            f,
            param,
            x,
            h_opt,
            &mut r_opt,
            &mut round_opt,
            &mut trunc_opt,
        );
        let error_opt = round_opt + trunc_opt;

        if error_opt < error && unsafe { (r_opt - r_0).abs() } < 4f64 * error {
            r_0 = r_opt;
            error = error_opt;
        }
    }
    *result = r_0;
    *abs_err = error;
    ::Value::Success
}

fn forward_deriv<T>(
    f: ::function<T>,
    param: &mut T,
    x: f64,
    h: f64,
    result: &mut f64,
    abs_err_round: &mut f64,
    abs_err_trunc: &mut f64,
) {
    let f1 = f(x + h / 4f64, param);
    let f2 = f(x + h / 2f64, param);
    let f3 = f(x + (3f64 / 4f64) * h, param);
    let f4 = f(x + h, param);

    let r2 = 2f64 * (f4 - f2);
    let r4 = (22f64 / 3f64) * (f4 - f3) - (62f64 / 3f64) * (f3 - f2) + (52f64 / 3f64) * (f2 - f1);

    let e4 =
        unsafe { 2f64 * 20.67f64 * (f4.abs() + f3.abs() + f2.abs() + f1.abs()) * ::DBL_EPSILON };

    let dy = unsafe { my_max((r2 / h).abs(), (r4 / h).abs()) * (x.abs() / h) * ::DBL_EPSILON };

    *result = r4 / h;
    *abs_err_trunc = unsafe { ((r4 - r2) / h).abs() };
    *abs_err_round = unsafe { (e4 / h).abs() + dy };
}

/// This function computes the numerical derivative of the function f at the point x using an adaptive forward difference algorithm with a step-size
/// of h. The function is evaluated only at points greater than x, and never at x itself. The derivative is returned in result and an estimate
/// of its absolute error is returned in abserr. This function should be used if f(x) has a discontinuity at x, or is undefined for values less
/// than x.
///
/// The initial value of h is used to estimate an optimal step-size, based on the scaling of the truncation error and round-off error in the
/// derivative calculation. The derivative at x is computed using an “open” 4-point rule for equally spaced abscissae at x+h/4, x+h/2, x+3h/4,
/// x+h, with an error estimate taken from the difference between the 4-point rule and the corresponding 2-point rule x+h/2, x+h.
pub fn deriv_forward<T>(
    f: ::function<T>,
    param: &mut T,
    x: f64,
    h: f64,
    result: &mut f64,
    abs_err: &mut f64,
) -> enums::Value {
    let mut r_0 = 0f64;
    let mut round = 0f64;
    let mut trunc = 0f64;

    forward_deriv(f, param, x, h, &mut r_0, &mut round, &mut trunc);
    let mut error = round + trunc;

    if round < trunc && (round > 0f64 && trunc > 0f64) {
        let mut r_opt = 0f64;
        let mut round_opt = 0f64;
        let mut trunc_opt = 0f64;

        let h_opt = unsafe { h * (round / trunc).powf(1f64 / 2f64) };
        forward_deriv(
            f,
            param,
            x,
            h_opt,
            &mut r_opt,
            &mut round_opt,
            &mut trunc_opt,
        );
        let error_opt = round_opt + trunc_opt;

        if error_opt < error && unsafe { (r_opt - r_0).abs() } < 4f64 * error {
            r_0 = r_opt;
            error = error_opt;
        }
    }
    *result = r_0;
    *abs_err = error;
    ::Value::Success
}

/// This function computes the numerical derivative of the function f at the point x using an adaptive backward difference algorithm with a
/// step-size of h. The function is evaluated only at points less than x, and never at x itself. The derivative is returned in result and an
/// estimate of its absolute error is returned in abserr. This function should be used if f(x) has a discontinuity at x, or is undefined for
/// values greater than x.
///
/// This function is equivalent to calling gsl_deriv_forward with a negative step-size.
pub fn deriv_backward<T>(
    f: ::function<T>,
    param: &mut T,
    x: f64,
    h: f64,
    result: &mut f64,
    abs_err: &mut f64,
) -> enums::Value {
    deriv_forward(f, param, x, -h, result, abs_err)
}
