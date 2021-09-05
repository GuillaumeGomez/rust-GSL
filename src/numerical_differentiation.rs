//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Numerical Differentiation

The functions described in this chapter compute numerical derivatives by finite differencing. An adaptive algorithm is used to find the
best choice of finite difference and to estimate the error in the derivative.

## References and Further Reading

The algorithms used by these functions are described in the following sources:

Abramowitz and Stegun, Handbook of Mathematical Functions, Section 25.3.4, and Table 25.5 (Coefficients for Differentiation).
S.D. Conte and Carl de Boor, Elementary Numerical Analysis: An Algorithmic Approach, McGraw-Hill, 1972.
!*/

use crate::Value;

/// This function computes the numerical derivative of the function f at the point x using an
/// adaptive central difference algorithm with a step-size of h. The derivative is returned in
/// result and an estimate of its absolute error is returned in abserr.
///
/// The initial value of h is used to estimate an optimal step-size, based on the scaling of the
/// truncation error and round-off error in the derivative calculation. The derivative is computed
/// using a 5-point rule for equally spaced abscissae at x-h, x-h/2, x, x+h/2, x+h, with an error
/// estimate taken from the difference between the 5-point rule and the corresponding 3-point rule
/// x-h, x, x+h. Note that the value of the function at x does not contribute to the derivative
/// calculation, so only 4-points are actually used.
///
/// Returns `(result, abs_err)`.
#[doc(alias = "gsl_deriv_central")]
pub fn deriv_central<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> (Value, f64, f64) {
    let mut result = 0.;
    let mut abs_err = 0.;
    let function = wrap_callback!(f, F);

    let ret = unsafe { sys::gsl_deriv_central(&function, x, h, &mut result, &mut abs_err) };
    (::Value::from(ret), result, abs_err)
}

/// This function computes the numerical derivative of the function f at the point x using an
/// adaptive forward difference algorithm with a step-size of h. The function is evaluated only at
/// points greater than x, and never at x itself. The derivative is returned in result and an
/// estimate of its absolute error is returned in abserr. This function should be used if f(x) has a
/// discontinuity at x, or is undefined for values less than x.
///
/// The initial value of h is used to estimate an optimal step-size, based on the scaling of the
/// truncation error and round-off error in the derivative calculation. The derivative at x is
/// computed using an “open” 4-point rule for equally spaced abscissae at x+h/4, x+h/2, x+3h/4,
/// x+h, with an error estimate taken from the difference between the 4-point rule and the
/// corresponding 2-point rule x+h/2, x+h.
///
/// Returns `(result, abs_err)`.
#[doc(alias = "gsl_deriv_forward")]
pub fn deriv_forward<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> (Value, f64, f64) {
    let mut result = 0.;
    let mut abs_err = 0.;
    let function = wrap_callback!(f, F);

    let ret = unsafe { sys::gsl_deriv_forward(&function, x, h, &mut result, &mut abs_err) };
    (::Value::from(ret), result, abs_err)
}

/// This function computes the numerical derivative of the function f at the point x using an
/// adaptive backward difference algorithm with a step-size of h. The function is evaluated only at
/// points less than x, and never at x itself. The derivative is returned in result and an estimate
/// of its absolute error is returned in abserr. This function should be used if f(x) has a
/// discontinuity at x, or is undefined for values greater than x.
///
/// This function is equivalent to calling gsl_deriv_forward with a negative step-size.
///
/// Returns `(result, abs_err)`.
#[doc(alias = "gsl_deriv_backward")]
pub fn deriv_backward<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> (Value, f64, f64) {
    let mut result = 0.;
    let mut abs_err = 0.;
    let function = wrap_callback!(f, F);

    let ret = unsafe { sys::gsl_deriv_backward(&function, x, h, &mut result, &mut abs_err) };
    (::Value::from(ret), result, abs_err)
}
