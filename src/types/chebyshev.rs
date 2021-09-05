//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Chebyshev Approximations

This chapter describes routines for computing Chebyshev approximations to univariate functions. A
Chebyshev approximation is a truncation of the series f(x) = \sum c_n T_n(x), where the Chebyshev
polynomials T_n(x) = \cos(n \arccos x) provide an orthogonal basis of polynomials on the interval
[-1,1] with the weight function 1 / \sqrt{1-x^2}. The first few Chebyshev polynomials are,
T_0(x) = 1, T_1(x) = x, T_2(x) = 2 x^2 - 1.

For further information see Abramowitz & Stegun, Chapter 22.

## Definitions

The approximation is made over the range [a,b] using order+1 terms, including the coefficient
`c[0]`. The series is computed using the following convention,

f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)

which is needed when accessing the coefficients directly.

## References and Further Reading

The following paper describes the use of Chebyshev series,

R. Broucke, `Ten Subroutines for the Manipulation of Chebyshev Series [C1] (Algorithm 446)`.
Communications of the ACM 16(4), 254–256 (1973)
!*/

use crate::Value;
use ffi::FFI;

ffi_wrapper!(ChebSeries, *mut sys::gsl_cheb_series, gsl_cheb_free);

impl ChebSeries {
    #[doc(alias = "gsl_cheb_alloc")]
    pub fn new(n: usize) -> Option<Self> {
        let tmp = unsafe { sys::gsl_cheb_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function computes the Chebyshev approximation cs for the function f over the range
    /// (a,b) to the previously specified order. The computation of the Chebyshev approximation is
    /// an O(n^2) process, and requires n function evaluations.
    #[doc(alias = "gsl_cheb_init")]
    pub fn init<F: Fn(f64) -> f64>(&mut self, f: F, a: f64, b: f64) -> Value {
        let function = wrap_callback!(f, F);

        Value::from(unsafe { sys::gsl_cheb_init(self.unwrap_unique(), &function, a, b) })
    }

    /// This function returns the order of Chebyshev series cs.
    #[doc(alias = "gsl_cheb_order")]
    pub fn order(&self) -> usize {
        unsafe { sys::gsl_cheb_order(self.unwrap_shared()) }
    }

    /// This function returns the size of the Chebyshev coefficient array c[] for the Chebyshev
    /// series cs.
    #[doc(alias = "gsl_cheb_size")]
    pub fn size(&self) -> usize {
        unsafe { sys::gsl_cheb_size(self.unwrap_shared()) }
    }

    /// This function evaluates the Chebyshev series cs at a given point x.
    #[doc(alias = "gsl_cheb_eval")]
    pub fn eval(&self, x: f64) -> f64 {
        unsafe { sys::gsl_cheb_eval(self.unwrap_shared(), x) }
    }

    /// This function computes the Chebyshev series cs at a given point x, estimating both the
    /// series result and its absolute error abserr. The error estimate is made from the first
    /// neglected term in the series.
    ///
    /// Returns `(result, abs_err)`.
    #[doc(alias = "gsl_cheb_eval_err")]
    pub fn eval_err(&self, x: f64) -> (Value, f64, f64) {
        let mut result = 0.;
        let mut abs_err = 0.;

        let ret =
            unsafe { sys::gsl_cheb_eval_err(self.unwrap_shared(), x, &mut result, &mut abs_err) };
        (::Value::from(ret), result, abs_err)
    }

    /// This function evaluates the Chebyshev series cs at a given point x, to (at most) the given
    /// order order.
    #[doc(alias = "gsl_cheb_eval_n")]
    pub fn eval_n(&self, order: usize, x: f64) -> f64 {
        unsafe { sys::gsl_cheb_eval_n(self.unwrap_shared(), order, x) }
    }

    /// This function evaluates a Chebyshev series cs at a given point x, estimating both the series
    /// result and its absolute error abserr, to (at most) the given order order. The error estimate
    /// is made from the first neglected term in the series.
    ///
    /// Returns `(result, abs_err)`.
    #[doc(alias = "gsl_cheb_eval_n_err")]
    pub fn eval_n_err(&self, order: usize, x: f64) -> (Value, f64, f64) {
        let mut result = 0.;
        let mut abs_err = 0.;

        let ret = unsafe {
            sys::gsl_cheb_eval_n_err(self.unwrap_shared(), order, x, &mut result, &mut abs_err)
        };
        (::Value::from(ret), result, abs_err)
    }

    /// This function computes the derivative of the series cs, storing the derivative coefficients
    /// in the previously allocated deriv. The two series cs and deriv must have been allocated with
    /// the same order.
    #[doc(alias = "gsl_cheb_calc_deriv")]
    pub fn calc_deriv(&self, deriv: &mut ChebSeries) -> Value {
        Value::from(unsafe {
            sys::gsl_cheb_calc_deriv(deriv.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function computes the integral of the series cs, storing the integral coefficients in
    /// the previously allocated integ. The two series cs and integ must have been allocated with
    /// the same order. The lower limit of the integration is taken to be the left hand end of the
    /// range a.
    #[doc(alias = "gsl_cheb_calc_integ")]
    pub fn calc_integ(&self, integ: &mut ChebSeries) -> Value {
        Value::from(unsafe {
            sys::gsl_cheb_calc_integ(integ.unwrap_unique(), self.unwrap_shared())
        })
    }
}
