//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Chebyshev Approximations

This chapter describes routines for computing Chebyshev approximations to univariate functions. A Chebyshev approximation is a truncation of the 
series f(x) = \sum c_n T_n(x), where the Chebyshev polynomials T_n(x) = \cos(n \arccos x) provide an orthogonal basis of polynomials on the 
interval [-1,1] with the weight function 1 / \sqrt{1-x^2}. The first few Chebyshev polynomials are, T_0(x) = 1, T_1(x) = x, T_2(x) = 2 x^2 - 1. 
For further information see Abramowitz & Stegun, Chapter 22.

##Definitions

The approximation is made over the range [a,b] using order+1 terms, including the coefficient c[0]. The series is computed using the following convention,

f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)
which is needed when accessing the coefficients directly.

##References and Further Reading

The following paper describes the use of Chebyshev series,

R. Broucke, “Ten Subroutines for the Manipulation of Chebyshev Series [C1] (Algorithm 446)”. Communications of the ACM 16(4), 254–256 (1973)
!*/

use ffi;
use enums;
use std::f64::consts::PI;
use std::num::Float;
use c_vec::CVec;

pub struct ChebSeries<'a> {
    c: *mut ffi::gsl_cheb_series,
    data: CVec<'a, f64>
}

impl<'a> ChebSeries<'a> {
    pub fn new(n: u64) -> Option<ChebSeries<'a>> {
        let tmp = unsafe { ffi::gsl_cheb_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                Some(ChebSeries {
                    c: tmp,
                    data: CVec::new((*tmp).c, (*tmp).order as usize + 1)
                })
            }
        }
    }

    /// This function computes the Chebyshev approximation cs for the function f over the range (a,b) to the previously specified order. The
    /// computation of the Chebyshev approximation is an O(n^2) process, and requires n function evaluations.
    pub fn init<T>(&mut self, func: ::function<T>, a: f64, b: f64, param: &mut T) -> enums::value::Value {
        if a >= b {
            rgsl_error!("null function interval [a,b]", ::Value::Dom);
            ::Value::Failure
        } else {
            unsafe {
                (*self.c).a = a;
                (*self.c).b = b;

                let bma = 0.5 * (b - a);
                let bpa = 0.5 * (b + a);
                let fac = 2.0 / ((*self.c).order as f64 + 1.0);

                let mut tmp_vec = CVec::new((*self.c).f, (*self.c).order_sp as usize + 1);
                for k in range(0, (*self.c).order + 1) {
                    let y = (PI * (k as f64 + 0.5) / ((*self.c).order as f64 + 1f64)).cos();
                    tmp_vec.as_mut_slice()[k as usize] = func(y * bma + bpa, param);
                }

                for j in range(0, (*self.c).order + 1) {
                    let mut sum = 0f64;

                    for k in range(0, (*self.c).order + 1) {
                        sum += tmp_vec.as_slice()[k as usize] * (PI * j as f64 * (k as f64 + 0.5) / ((*self.c).order as f64 + 1f64)).cos();
                    }
                    self.data.as_mut_slice()[j as usize] = fac * sum;
                }
            }
            ::Value::Success
        }
    }

    /// This function returns the order of Chebyshev series cs.
    pub fn order(&self) -> u64 {
        unsafe { ffi::gsl_cheb_order(self.c as *const ffi::gsl_cheb_series) }
    }

    /// This function returns the size of the Chebyshev coefficient array c[] for the Chebyshev series cs.
    pub fn size(&self) -> u64 {
        unsafe { ffi::gsl_cheb_size(self.c as *const ffi::gsl_cheb_series) }
    }

    /// This function returns a pointer to the coefficient array c[] location in memory for the Chebyshev series cs.
    pub fn as_slice<'r>(&'r self) -> &'r [f64] {
        self.data.as_slice()
    }

    /// This function returns a pointer to the coefficient array c[] location in memory for the Chebyshev series cs.
    pub fn as_mut_slice<'r>(&'r mut self) -> &'r mut [f64] {
        self.data.as_mut_slice()
    }

    /// This function evaluates the Chebyshev series cs at a given point x.
    pub fn eval(&self, x: f64) -> f64 {
        unsafe { ffi::gsl_cheb_eval(self.c as *const ffi::gsl_cheb_series, x) }
    }

    /// This function computes the Chebyshev series cs at a given point x, estimating both the series result and its absolute error abserr.
    /// The error estimate is made from the first neglected term in the series.
    pub fn eval_err(&self, x: f64, result: &mut f64, abs_err: &mut f64) -> enums::value::Value {
        unsafe { ffi::gsl_cheb_eval_err(self.c as *const ffi::gsl_cheb_series, x, result, abs_err) }
    }

    /// This function evaluates the Chebyshev series cs at a given point x, to (at most) the given order order.
    pub fn eval_n(&self, order: u64, x: f64) -> f64 {
        unsafe { ffi::gsl_cheb_eval_n(self.c as *const ffi::gsl_cheb_series, order, x) }
    }

    /// This function evaluates a Chebyshev series cs at a given point x, estimating both the series result and its absolute error abserr, to
    /// (at most) the given order order. The error estimate is made from the first neglected term in the series.
    pub fn eval_n_err(&self, order: u64, x: f64, result: &mut f64, abs_err: &mut f64) -> enums::value::Value {
        unsafe { ffi::gsl_cheb_eval_n_err(self.c as *const ffi::gsl_cheb_series, order, x, result, abs_err) }
    }

    /// This function computes the derivative of the series cs, storing the derivative coefficients in the previously allocated deriv. The
    /// two series cs and deriv must have been allocated with the same order.
    pub fn calc_deriv(&self, deriv: &ChebSeries) -> enums::value::Value {
        unsafe { ffi::gsl_cheb_calc_deriv(deriv.c, self.c as *const ffi::gsl_cheb_series) }
    }

    /// This function computes the integral of the series cs, storing the integral coefficients in the previously allocated integ. The two series
    /// cs and integ must have been allocated with the same order. The lower limit of the integration is taken to be the left hand end of the range a.
    pub fn calc_integ(&self, integ: &ChebSeries) -> enums::value::Value {
        unsafe { ffi::gsl_cheb_calc_integ(integ.c, self.c as *const ffi::gsl_cheb_series) }
    }
}

#[unsafe_destructor]
impl<'a> Drop for ChebSeries<'a> {
    fn drop(&mut self) {
        unsafe { ffi::gsl_cheb_free(self.c) };
        self.c = ::std::ptr::null_mut();
    }
}

impl<'a> ffi::FFI<ffi::gsl_cheb_series> for ChebSeries<'a> {
    fn wrap(c: *mut ffi::gsl_cheb_series) -> ChebSeries<'a> {
        unsafe {
            ChebSeries {
                c: c,
                data: CVec::new((*c).c, (*c).order as usize + 1)
            }
        }
    }

    fn unwrap(c: &ChebSeries) -> *mut ffi::gsl_cheb_series {
        c.c
    }
}
