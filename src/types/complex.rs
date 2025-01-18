//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use num_complex::Complex;

#[deprecated(since = "8.0.0", note = "use `Complex<f64>` instead")]
pub type ComplexF64 = Complex<f64>;
#[deprecated(since = "8.0.0", note = "use `Complex<f32>` instead")]
pub type ComplexF32 = Complex<f32>;

pub(crate) trait ToC<T> {
    fn unwrap(self) -> T;
}

pub(crate) trait FromC<T> {
    fn wrap(self) -> T;
}

impl ToC<sys::gsl_complex> for Complex<f64> {
    fn unwrap(self) -> sys::gsl_complex {
        // Complex<T> is memory layout compatible with [T; 2]
        unsafe { std::mem::transmute(self) }
    }
}

impl FromC<Complex<f64>> for sys::gsl_complex {
    fn wrap(self) -> Complex<f64> {
        unsafe { std::mem::transmute(self) }
    }
}

/// Define the capabilities provided for Complex numbers by the GSL.
pub trait ComplexOps<T> {
    /// This function uses the rectangular Cartesian components (x,y)
    /// to return the complex number z = x + i y.
    #[doc(alias = "gsl_complex_rect")]
    fn rect(x: T, y: T) -> Complex<T>;

    /// This function returns the complex number z = r \exp(iθ) = r
    /// (\cos(θ) + i \sin(θ)) from the polar representation (`r`, `theta`).
    #[doc(alias = "gsl_complex_polar")]
    fn polar(r: T, theta: T) -> Complex<T>;

    /// This function returns the magnitude of the complex number z, |z|.
    #[doc(alias = "gsl_complex_abs")]
    #[deprecated(since = "8.0.0", note = "please use `.norm()` instead")]
    fn abs(&self) -> T;

    /// This function returns the squared magnitude of the complex
    /// number z = `self`, |z|².
    #[doc(alias = "gsl_complex_abs2")]
    #[deprecated(since = "8.0.0", note = "please use `.norm_sqr()` instead")]
    fn abs2(&self) -> T;

    /// This function returns the natural logarithm of the magnitude
    /// of the complex number z = `self`, log|z|.
    ///
    /// It allows an accurate evaluation of \log|z| when |z| is close to one.
    ///
    /// The direct evaluation of log(gsl_complex_abs(z)) would lead to
    /// a loss of precision in this case.
    #[doc(alias = "gsl_complex_logabs")]
    fn logabs(&self) -> T;

    /// This function returns the sum of the complex numbers a and b, z=a+b.
    #[doc(alias = "gsl_complex_add")]
    #[deprecated(since = "8.0.0", note = "please use `+` instead")]
    fn add(&self, other: &Complex<T>) -> Complex<T>;

    /// This function returns the difference of the complex numbers a
    /// and b, z=a-b.
    #[doc(alias = "gsl_complex_sub")]
    #[deprecated(since = "8.0.0", note = "please use `-` instead")]
    fn sub(&self, other: &Complex<T>) -> Complex<T>;

    /// This function returns the product of the complex numbers a and b, z=ab.
    #[doc(alias = "gsl_complex_mul")]
    #[deprecated(since = "8.0.0", note = "please use `*` instead")]
    fn mul(&self, other: &Complex<T>) -> Complex<T>;

    /// This function returns the quotient of the complex numbers a
    /// and b, z=a/b.
    #[doc(alias = "gsl_complex_div")]
    #[deprecated(since = "8.0.0", note = "please use `/` of `fdiv` instead")]
    fn div(&self, other: &Complex<T>) -> Complex<T>;

    /// This function returns the sum of the complex number a and the
    /// real number x, z = a + x.
    #[doc(alias = "gsl_complex_add_real")]
    #[deprecated(since = "8.0.0", note = "please use `+` instead")]
    fn add_real(&self, x: T) -> Complex<T>;

    /// This function returns the difference of the complex number a
    /// and the real number x, z=a-x.
    #[doc(alias = "gsl_complex_sub_real")]
    #[deprecated(since = "8.0.0", note = "please use `-` instead")]
    fn sub_real(&self, x: T) -> Complex<T>;

    /// This function returns the product of the complex number a and
    /// the real number x, z=ax.
    #[doc(alias = "gsl_complex_mul_real")]
    #[deprecated(since = "8.0.0", note = "please use `*` instead")]
    fn mul_real(&self, x: T) -> Complex<T>;

    /// This function returns the quotient of the complex number a and
    /// the real number x, z=a/x.
    #[doc(alias = "gsl_complex_div_real")]
    #[deprecated(since = "8.0.0", note = "please use `/` instead")]
    fn div_real(&self, x: T) -> Complex<T>;

    /// This function returns the sum of the complex number a and the
    /// imaginary number iy, z=a+iy.
    #[doc(alias = "gsl_complex_add_imag")]
    #[deprecated(since = "8.0.0", note = "please use `self + x * Complex::I` instead")]
    fn add_imag(&self, x: T) -> Complex<T>;

    /// This function returns the difference of the complex number a
    /// and the imaginary number iy, z=a-iy.
    #[doc(alias = "gsl_complex_sub_imag")]
    #[deprecated(since = "8.0.0", note = "please use `self - x * Complex::I` instead")]
    fn sub_imag(&self, x: T) -> Complex<T>;

    /// This function returns the product of the complex number a and
    /// the imaginary number iy, z=a*(iy).
    #[doc(alias = "gsl_complex_mul_imag")]
    #[deprecated(since = "8.0.0", note = "please use `self * x * Complex::I` instead")]
    fn mul_imag(&self, x: T) -> Complex<T>;

    /// This function returns the quotient of the complex number a and
    /// the imaginary number iy, z=a/(iy).
    #[doc(alias = "gsl_complex_div_imag")]
    #[deprecated(since = "8.0.0", note = "please use `self / (x * Complex::I)` instead")]
    fn div_imag(&self, x: T) -> Complex<T>;

    /// This function returns the complex conjugate of the complex
    /// number z, z^* = x - i y.
    #[doc(alias = "gsl_complex_conjugate")]
    #[deprecated(since = "8.0.0", note = "please use `.conj()` instead")]
    fn conjugate(&self) -> Complex<T>;

    /// This function returns the inverse, or reciprocal, of the
    /// complex number z, 1/z = (x - i y)/ (x^2 + y^2).
    #[doc(alias = "gsl_complex_inverse")]
    #[deprecated(since = "8.0.0", note = "please use `.inv()` instead")]
    fn inverse(&self) -> Complex<T>;

    /// This function returns the negative of the complex number z, -z
    /// = (-x) + i(-y).
    #[doc(alias = "gsl_complex_negative")]
    #[deprecated(since = "8.0.0", note = "please use the unary `-` instead")]
    fn negative(&self) -> Complex<T>;

    /// This function returns the complex square root of the real
    /// number x, where x may be negative.
    #[doc(alias = "gsl_complex_sqrt_real")]
    fn sqrt_real(x: T) -> Complex<T>;

    /// The function returns the complex number z raised to the
    /// complex power a z^a.  This is computed as \exp(\log(z)*a)
    /// using complex logarithms and complex exponentials.
    #[doc(alias = "gsl_complex_pow")]
    #[deprecated(since = "8.0.0", note = "please use the unary `-` instead")]
    fn pow(&self, other: &Complex<T>) -> Complex<T>;

    /// This function returns the complex number z raised to the real
    /// power x, z^x.
    #[doc(alias = "gsl_complex_pow_real")]
    #[deprecated(since = "8.0.0", note = "please use `.powf(x)` instead")]
    fn pow_real(&self, x: T) -> Complex<T>;

    /// This function returns the complex base-b logarithm of the
    /// complex number z, \log_b(z).  This quantity is computed as the
    /// ratio \log(z)/\log(b).
    #[doc(alias = "gsl_complex_log_b")]
    #[deprecated(since = "8.0.0", note = "please use `.log(base)` instead")]
    fn log_b(&self, base: &Complex<T>) -> Complex<T>;

    /// This function returns the complex secant of the complex number
    /// z, \sec(z) = 1/\cos(z).
    #[doc(alias = "gsl_complex_sec")]
    fn sec(&self) -> Complex<T>;

    /// This function returns the complex cosecant of the complex
    /// number z, \csc(z) = 1/\sin(z).
    #[doc(alias = "gsl_complex_csc")]
    fn csc(&self) -> Complex<T>;

    /// This function returns the complex cotangent of the complex
    /// number z, \cot(z) = 1/\tan(z).
    #[doc(alias = "gsl_complex_cot")]
    fn cot(&self) -> Complex<T>;

    /// This function returns the complex arcsine of the complex
    /// number z, \arcsin(z).  The branch cuts are on the real axis,
    /// less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arcsin")]
    #[deprecated(since = "8.0.0", note = "please use `.asin()` instead")]
    fn arcsin(&self) -> Complex<T>;

    /// This function returns the complex arcsine of the real number
    /// z, \arcsin(z).
    ///
    /// * For z between -1 and 1, the function returns a real value in
    ///   the range \[-π/2,π/2\].
    /// * For z less than -1 the result has a real part of -π/2 and
    ///   a positive imaginary part.
    /// * For z greater than 1 the result has a real part of π/2 and
    ///   a negative imaginary part.
    #[doc(alias = "gsl_complex_arcsin_real")]
    fn arcsin_real(z: T) -> Complex<T>;

    /// This function returns the complex arccosine of the complex
    /// number z, \arccos(z).  The branch cuts are on the real axis,
    /// less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arccos")]
    #[deprecated(since = "8.0.0", note = "please use `.acos()` instead")]
    fn arccos(&self) -> Complex<T>;

    /// This function returns the complex arccosine of the real number
    /// z, \arccos(z).
    ///
    /// * For z between -1 and 1, the function returns a real value in
    ///   the range \[0,π\].
    /// * For z less than -1 the result has a real part of \pi and a
    ///   negative imaginary part.
    /// * For z greater than 1 the result is purely imaginary and positive.
    #[doc(alias = "gsl_complex_arccos_real")]
    fn arccos_real(z: T) -> Complex<T>;

    /// This function returns the complex arctangent of the complex
    /// number z, \arctan(z).  The branch cuts are on the imaginary
    /// axis, below -i and above i.
    #[doc(alias = "gsl_complex_arctan")]
    #[deprecated(since = "8.0.0", note = "please use `.atan()` instead")]
    fn arctan(&self) -> Complex<T>;

    /// This function returns the complex arcsecant of the complex
    /// number z, \arcsec(z) = \arccos(1/z).
    #[doc(alias = "gsl_complex_arcsec")]
    fn arcsec(&self) -> Complex<T>;

    /// This function returns the complex arcsecant of the real number
    /// z, \arcsec(z) = \arccos(1/z).
    #[doc(alias = "gsl_complex_arcsec_real")]
    fn arcsec_real(z: T) -> Complex<T>;

    /// This function returns the complex arccosecant of the complex
    /// number z, \arccsc(z) = \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsc")]
    fn arccsc(&self) -> Complex<T>;

    /// This function returns the complex arccosecant of the real
    /// number z, \arccsc(z) = \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsc_real")]
    fn arccsc_real(z: T) -> Complex<T>;

    /// This function returns the complex arccotangent of the complex
    /// number z, \arccot(z) = \arctan(1/z).
    #[doc(alias = "gsl_complex_arccot")]
    fn arccot(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic secant of the
    /// complex number z, \sech(z) = 1/\cosh(z).
    #[doc(alias = "gsl_complex_sech")]
    fn sech(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic cosecant of the
    /// complex number z, \csch(z) = 1/\sinh(z).
    #[doc(alias = "gsl_complex_csch")]
    fn csch(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic cotangent of the
    /// complex number z, \coth(z) = 1/\tanh(z).
    #[doc(alias = "gsl_complex_coth")]
    fn coth(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic arcsine of the
    /// complex number z, \arcsinh(z).  The branch cuts are on the
    /// imaginary axis, below -i and above i.
    #[doc(alias = "gsl_complex_arcsinh")]
    #[deprecated(since = "8.0.0", note = "please use `.asinh()` instead")]
    fn arcsinh(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic arccosine of the
    /// complex number z, \arccosh(z).  The branch cut is on the real
    /// axis, less than 1.  Note that in this case we use the negative
    /// square root in formula 4.6.21 of Abramowitz & Stegun giving
    /// \arccosh(z)=\log(z-\sqrt{z^2-1}).
    #[doc(alias = "gsl_complex_arccosh")]
    #[deprecated(since = "8.0.0", note = "please use `.acosh()` instead")]
    fn arccosh(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic arccosine of the
    /// real number z, \arccosh(z).
    #[doc(alias = "gsl_complex_arccosh_real")]
    fn arccosh_real(z: T) -> Complex<T>;

    /// This function returns the complex hyperbolic arctangent of the
    /// complex number z, \arctanh(z).
    ///
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arctanh")]
    #[deprecated(since = "8.0.0", note = "please use `.atanh()` instead")]
    fn arctanh(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic arctangent of the
    /// real number z, \arctanh(z).
    #[doc(alias = "gsl_complex_arctanh_real")]
    fn arctanh_real(z: T) -> Complex<T>;

    /// This function returns the complex hyperbolic arcsecant of the
    /// complex number z, \arcsech(z) = \arccosh(1/z).
    #[doc(alias = "gsl_complex_arcsech")]
    fn arcsech(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic arccosecant of
    /// the complex number z, \arccsch(z) = \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsch")]
    fn arccsch(&self) -> Complex<T>;

    /// This function returns the complex hyperbolic arccotangent of
    /// the complex number z, \arccoth(z) = \arctanh(1/z).
    #[doc(alias = "gsl_complex_arccoth")]
    fn arccoth(&self) -> Complex<T>;

    #[deprecated(since = "8.0.0", note = "please use `.re` instead")]
    fn real(&self) -> T;

    #[deprecated(since = "8.0.0", note = "please use `.im` instead")]
    fn imaginary(&self) -> T;
}

impl ComplexOps<f64> for Complex<f64> {
    #[doc(alias = "gsl_complex_rect")]
    fn rect(x: f64, y: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_rect(x, y).wrap() }
    }

    #[doc(alias = "gsl_complex_polar")]
    fn polar(r: f64, theta: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_polar(r, theta).wrap() }
    }

    #[doc(alias = "gsl_complex_abs")]
    fn abs(&self) -> f64 {
        unsafe { sys::gsl_complex_abs(self.unwrap()) }
    }

    #[doc(alias = "gsl_complex_abs2")]
    fn abs2(&self) -> f64 {
        unsafe { sys::gsl_complex_abs2(self.unwrap()) }
    }

    #[doc(alias = "gsl_complex_logabs")]
    fn logabs(&self) -> f64 {
        unsafe { sys::gsl_complex_logabs(self.unwrap()) }
    }

    #[doc(alias = "gsl_complex_add")]
    fn add(&self, other: &Complex<f64>) -> Complex<f64> {
        unsafe { sys::gsl_complex_add(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sub")]
    fn sub(&self, other: &Complex<f64>) -> Complex<f64> {
        unsafe { sys::gsl_complex_sub(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_mul")]
    fn mul(&self, other: &Complex<f64>) -> Complex<f64> {
        unsafe { sys::gsl_complex_mul(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_div")]
    fn div(&self, other: &Complex<f64>) -> Complex<f64> {
        unsafe { sys::gsl_complex_div(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_add_real")]
    fn add_real(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_add_real(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_sub_real")]
    fn sub_real(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_sub_real(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_mul_real")]
    fn mul_real(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_mul_real(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_div_real")]
    fn div_real(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_div_real(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_add_imag")]
    fn add_imag(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_add_imag(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_sub_imag")]
    fn sub_imag(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_sub_imag(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_mul_imag")]
    fn mul_imag(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_mul_imag(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_div_imag")]
    fn div_imag(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_div_imag(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_conjugate")]
    fn conjugate(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_conjugate(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_inverse")]
    fn inverse(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_inverse(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_negative")]
    fn negative(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_negative(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sqrt_real")]
    fn sqrt_real(x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_sqrt_real(x).wrap() }
    }

    #[doc(alias = "gsl_complex_pow")]
    fn pow(&self, other: &Complex<f64>) -> Complex<f64> {
        unsafe { sys::gsl_complex_pow(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_pow_real")]
    fn pow_real(&self, x: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_pow_real(self.unwrap(), x).wrap() }
    }

    #[doc(alias = "gsl_complex_log_b")]
    fn log_b(&self, other: &Complex<f64>) -> Complex<f64> {
        unsafe { sys::gsl_complex_log_b(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sec")]
    fn sec(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_sec(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_csc")]
    fn csc(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_csc(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_cot")]
    fn cot(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_cot(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsin")]
    fn arcsin(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arcsin(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsin_real")]
    fn arcsin_real(z: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_arcsin_real(z).wrap() }
    }

    #[doc(alias = "gsl_complex_arccos")]
    fn arccos(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccos(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccos_real")]
    fn arccos_real(z: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccos_real(z).wrap() }
    }

    #[doc(alias = "gsl_complex_arctan")]
    fn arctan(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arctan(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsec")]
    fn arcsec(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arcsec(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsec_real")]
    fn arcsec_real(z: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_arcsec_real(z).wrap() }
    }

    #[doc(alias = "gsl_complex_arccsc")]
    fn arccsc(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccsc(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccsc_real")]
    fn arccsc_real(z: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccsc_real(z).wrap() }
    }

    #[doc(alias = "gsl_complex_arccot")]
    fn arccot(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccot(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sech")]
    fn sech(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_sech(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_csch")]
    fn csch(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_csch(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_coth")]
    fn coth(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_coth(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsinh")]
    fn arcsinh(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arcsinh(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccosh")]
    fn arccosh(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccosh(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccosh_real")]
    fn arccosh_real(z: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccosh_real(z).wrap() }
    }

    #[doc(alias = "gsl_complex_arctanh")]
    fn arctanh(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arctanh(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arctanh_real")]
    fn arctanh_real(z: f64) -> Complex<f64> {
        unsafe { sys::gsl_complex_arctanh_real(z).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsech")]
    fn arcsech(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arcsech(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccsch")]
    fn arccsch(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccsch(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccoth")]
    fn arccoth(&self) -> Complex<f64> {
        unsafe { sys::gsl_complex_arccoth(self.unwrap()).wrap() }
    }

    fn real(&self) -> f64 {
        self.re
    }

    fn imaginary(&self) -> f64 {
        self.im
    }
}

// The GLS Complex module does not support `f32` operations.  Thus we
// convert back and forth to `f64`.
impl ToC<sys::gsl_complex> for Complex<f32> {
    fn unwrap(self) -> sys::gsl_complex {
        sys::gsl_complex {
            dat: [self.re as f64, self.im as f64],
        }
    }
}

// For use by other modules (e.g. `blas`).
impl FromC<Complex<f32>> for sys::gsl_complex_float {
    fn wrap(self) -> Complex<f32> {
        // Complex<T> is memory layout compatible with [T; 2]
        unsafe { std::mem::transmute(self) }
    }
}

impl ToC<sys::gsl_complex_float> for Complex<f32> {
    fn unwrap(self) -> sys::gsl_complex_float {
        unsafe { std::mem::transmute(self) }
    }
}

impl FromC<Complex<f32>> for sys::gsl_complex {
    fn wrap(self) -> Complex<f32> {
        let [re, im] = self.dat;
        Complex {
            re: re as f32,
            im: im as f32,
        }
    }
}

impl ComplexOps<f32> for Complex<f32> {
    #[doc(alias = "gsl_complex_rect")]
    fn rect(x: f32, y: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_rect(x as f64, y as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_polar")]
    fn polar(r: f32, theta: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_polar(r as f64, theta as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_abs")]
    fn abs(&self) -> f32 {
        unsafe { sys::gsl_complex_abs(self.unwrap()) as f32 }
    }

    #[doc(alias = "gsl_complex_abs2")]
    fn abs2(&self) -> f32 {
        unsafe { sys::gsl_complex_abs2(self.unwrap()) as f32 }
    }

    #[doc(alias = "gsl_complex_logabs")]
    fn logabs(&self) -> f32 {
        unsafe { sys::gsl_complex_logabs(self.unwrap()) as f32 }
    }

    #[doc(alias = "gsl_complex_add")]
    fn add(&self, other: &Complex<f32>) -> Complex<f32> {
        unsafe { sys::gsl_complex_add(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sub")]
    fn sub(&self, other: &Complex<f32>) -> Complex<f32> {
        unsafe { sys::gsl_complex_sub(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_mul")]
    fn mul(&self, other: &Complex<f32>) -> Complex<f32> {
        unsafe { sys::gsl_complex_mul(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_div")]
    fn div(&self, other: &Complex<f32>) -> Complex<f32> {
        unsafe { sys::gsl_complex_div(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_add_real")]
    fn add_real(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_add_real(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_sub_real")]
    fn sub_real(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_sub_real(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_mul_real")]
    fn mul_real(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_mul_real(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_div_real")]
    fn div_real(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_div_real(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_add_imag")]
    fn add_imag(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_add_imag(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_sub_imag")]
    fn sub_imag(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_sub_imag(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_mul_imag")]
    fn mul_imag(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_mul_imag(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_div_imag")]
    fn div_imag(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_div_imag(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_conjugate")]
    fn conjugate(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_conjugate(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_inverse")]
    fn inverse(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_inverse(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_negative")]
    fn negative(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_negative(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sqrt_real")]
    fn sqrt_real(x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_sqrt_real(x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_pow")]
    fn pow(&self, other: &Complex<f32>) -> Complex<f32> {
        unsafe { sys::gsl_complex_pow(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_pow_real")]
    fn pow_real(&self, x: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_pow_real(self.unwrap(), x as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_log_b")]
    fn log_b(&self, other: &Complex<f32>) -> Complex<f32> {
        unsafe { sys::gsl_complex_log_b(self.unwrap(), other.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sec")]
    fn sec(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_sec(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_csc")]
    fn csc(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_csc(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_cot")]
    fn cot(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_cot(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsin")]
    fn arcsin(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arcsin(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsin_real")]
    fn arcsin_real(z: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_arcsin_real(z as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_arccos")]
    fn arccos(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccos(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccos_real")]
    fn arccos_real(z: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccos_real(z as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_arctan")]
    fn arctan(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arctan(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsec")]
    fn arcsec(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arcsec(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsec_real")]
    fn arcsec_real(z: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_arcsec_real(z as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_arccsc")]
    fn arccsc(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccsc(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccsc_real")]
    fn arccsc_real(z: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccsc_real(z as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_arccot")]
    fn arccot(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccot(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_sech")]
    fn sech(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_sech(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_csch")]
    fn csch(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_csch(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_coth")]
    fn coth(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_coth(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsinh")]
    fn arcsinh(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arcsinh(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccosh")]
    fn arccosh(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccosh(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccosh_real")]
    fn arccosh_real(z: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccosh_real(z as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_arctanh")]
    fn arctanh(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arctanh(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arctanh_real")]
    fn arctanh_real(z: f32) -> Complex<f32> {
        unsafe { sys::gsl_complex_arctanh_real(z as f64).wrap() }
    }

    #[doc(alias = "gsl_complex_arcsech")]
    fn arcsech(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arcsech(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccsch")]
    fn arccsch(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccsch(self.unwrap()).wrap() }
    }

    #[doc(alias = "gsl_complex_arccoth")]
    fn arccoth(&self) -> Complex<f32> {
        unsafe { sys::gsl_complex_arccoth(self.unwrap()).wrap() }
    }

    fn real(&self) -> f32 {
        self.re
    }

    fn imaginary(&self) -> f32 {
        self.im
    }
}

#[cfg(test)]
mod tests {
    // All these tests have been tested against the following C code:
    //
    // ```ignore
    // #include <gsl/gsl_complex.h>
    // #include <gsl/gsl_complex_math.h>
    // #include <gsl/gsl_inline.h>
    // #include <gsl/gsl_complex.h>
    //
    // void print_complex(gsl_complex *c, const char *text) {
    //   printf("%s: %f %f\n", text, c->dat[0], c->dat[1]);
    // }
    //
    // int main (void)
    // {
    //   gsl_complex c = gsl_complex_rect(10., 10.);
    //   gsl_complex c2 = gsl_complex_rect(1., -1.);
    //   print_complex(&c, "rect");
    //   print_complex(&c2, "rect");
    //   gsl_complex c3 = gsl_complex_polar(5., 7.);
    //   print_complex(&c3, "polar");
    //
    //   printf("-> %f\n", gsl_complex_arg(c3));
    //   printf("-> %f\n", gsl_complex_abs(c3));
    //   printf("-> %f\n", gsl_complex_abs2(c3));
    //   printf("-> %f\n", gsl_complex_logabs(c3));
    //
    //   {
    //     gsl_complex c4 = gsl_complex_add(c3, c2);
    //     print_complex(&c4, "\nadd");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sub(c3, c2);
    //     print_complex(&c4, "sub");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_mul(c3, c2);
    //     print_complex(&c4, "mul");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_div(c3, c2);
    //     print_complex(&c4, "div");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_add_real(c3, 5.);
    //     print_complex(&c4, "add_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sub_real(c3, 5.);
    //     print_complex(&c4, "sub_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_mul_real(c3, 5.);
    //     print_complex(&c4, "mul_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_div_real(c3, 5.);
    //     print_complex(&c4, "div_real");
    //   }
    //
    //
    //   {
    //     gsl_complex c4 = gsl_complex_add_imag(c3, 5.);
    //     print_complex(&c4, "\nadd_imag");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sub_imag(c3, 5.);
    //     print_complex(&c4, "sub_imag");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_mul_imag(c3, 5.);
    //     print_complex(&c4, "mul_imag");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_div_imag(c3, 5.);
    //     print_complex(&c4, "div_imag");
    //   }
    //
    //
    //   {
    //     gsl_complex c4 = gsl_complex_conjugate(c3);
    //     print_complex(&c4, "\nconjugate");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_inverse(c3);
    //     print_complex(&c4, "inverse");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_negative(c3);
    //     print_complex(&c4, "negative");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sqrt(c3);
    //     print_complex(&c4, "sqrt");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sqrt_real(5.);
    //     print_complex(&c4, "sqrt_real");
    //   }
    //
    //
    //   {
    //     gsl_complex c4 = gsl_complex_pow(c3, c2);
    //     print_complex(&c4, "\npow");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_pow_real(c3, 5.);
    //     print_complex(&c4, "pow_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_exp(c3);
    //     print_complex(&c4, "exp");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_log(c3);
    //     print_complex(&c4, "log");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_log10(c3);
    //     print_complex(&c4, "log10");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_log_b(c3, c2);
    //     print_complex(&c4, "log_b");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sin(c3);
    //     print_complex(&c4, "sin");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_cos(c3);
    //     print_complex(&c4, "cos");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_tan(c3);
    //     print_complex(&c4, "tan");
    //   }
    //
    //
    //   {
    //     gsl_complex c4 = gsl_complex_sec(c3);
    //     print_complex(&c4, "\nsec");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_csc(c3);
    //     print_complex(&c4, "csc");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_cot(c3);
    //     print_complex(&c4, "cot");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arcsin(c3);
    //     print_complex(&c4, "arcsin");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arcsin_real(5.);
    //     print_complex(&c4, "arcsin_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccos(c3);
    //     print_complex(&c4, "arccos");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccos_real(5.);
    //     print_complex(&c4, "arccos_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arctan(c3);
    //     print_complex(&c4, "arctan");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arcsec(c3);
    //     print_complex(&c4, "arcsec");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arcsec_real(5.);
    //     print_complex(&c4, "arcsec_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccsc(c3);
    //     print_complex(&c4, "arccsc");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccsc_real(5.);
    //     print_complex(&c4, "arccsc_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccot(c3);
    //     print_complex(&c4, "arccot");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sinh(c3);
    //     print_complex(&c4, "sinh");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_cosh(c3);
    //     print_complex(&c4, "cosh");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_tanh(c3);
    //     print_complex(&c4, "tanh");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_sech(c3);
    //     print_complex(&c4, "sech");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_csch(c3);
    //     print_complex(&c4, "csch");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_coth(c3);
    //     print_complex(&c4, "coth");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arcsinh(c3);
    //     print_complex(&c4, "arcsinh");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccosh(c3);
    //     print_complex(&c4, "arccosh");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccosh_real(5.);
    //     print_complex(&c4, "arccosh_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arctanh(c3);
    //     print_complex(&c4, "arctanh");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arctanh_real(5.);
    //     print_complex(&c4, "arctanh_real");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arcsech(c3);
    //     print_complex(&c4, "arcsech");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccsch(c3);
    //     print_complex(&c4, "arccsch");
    //   }
    //   {
    //     gsl_complex c4 = gsl_complex_arccoth(c3);
    //     print_complex(&c4, "arccoth");
    //   }
    //   return 0;
    // }
    // ```

    #[test]
    #[allow(deprecated)]
    fn complex_f64() {
        type C = num_complex::Complex<f64>;
        use crate::complex::ComplexOps;

        let v = C::rect(10., 10.);
        assert_eq!(v, C { re: 10., im: 10. });
        let v2 = C::rect(1., -1.);
        assert_eq!(v2, C { re: 1., im: -1. });
        let v = C::polar(5., 7.);
        assert_eq!(
            format!("{:.4} {:.4}", v.re, v.im),
            "3.7695 3.2849".to_owned()
        );

        let arg = v.arg();
        assert_eq!(format!("{:.4}", arg), "0.7168".to_owned());
        let arg = v.abs();
        assert_eq!(format!("{:.3}", arg), "5.000".to_owned());
        let arg = v.abs2();
        assert_eq!(format!("{:.3}", arg), "25.000".to_owned());
        let arg = v.logabs();
        assert_eq!(format!("{:.4}", arg), "1.6094".to_owned());

        let v3 = v.add(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "4.7695 2.2849".to_owned()
        );
        let v3 = v.sub(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.7695 4.2849".to_owned()
        );
        let v3 = v.mul(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "7.0544 -0.4846".to_owned()
        );
        let v3 = v.div(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.2423 3.5272".to_owned()
        );
        let v3 = v.add_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "8.7695 3.2849".to_owned()
        );
        let v3 = v.sub_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-1.2305 3.2849".to_owned()
        );
        let v3 = v.mul_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "18.8476 16.4247".to_owned()
        );
        let v3 = v.div_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.7539 0.6570".to_owned()
        );

        let v3 = v.add_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "3.7695 8.2849".to_owned()
        );
        let v3 = v.sub_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "3.7695 -1.7151".to_owned()
        );
        let v3 = v.mul_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-16.4247 18.8476".to_owned()
        );
        let v3 = v.div_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.6570 -0.7539".to_owned()
        );

        let v3 = v.conjugate();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "3.7695 -3.2849".to_owned()
        );
        let v3 = v.inverse();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1508 -0.1314".to_owned()
        );
        let v3 = v.negative();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-3.7695 -3.2849".to_owned()
        );
        let v3 = v.sqrt();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.0940 0.7844".to_owned()
        );
        let v3 = C::sqrt_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.2361 0.0000".to_owned()
        );

        let v3 = v.pow(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "6.4240 -7.9737".to_owned()
        );
        let v3 = v.pow_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-2824.0381 -1338.0708".to_owned()
        );
        let v3 = v.exp();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-42.9142 -6.1938".to_owned()
        );
        let v3 = v.ln();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.6094 0.7168".to_owned()
        );
        let v3 = v.log10();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.6990 0.3113".to_owned()
        );
        let v3 = v.log_b(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0071 2.0523".to_owned()
        );
        let v3 = v.sin();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-7.8557 -10.7913".to_owned()
        );
        let v3 = v.cos();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-10.8216 7.8337".to_owned()
        );
        let v3 = v.tan();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.0027 0.9991".to_owned()
        );

        let v3 = v.sec();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0606 -0.0439".to_owned()
        );
        let v3 = v.csc();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0441 0.0606".to_owned()
        );
        let v3 = v.cot();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.0027 -1.0009".to_owned()
        );
        let v3 = v.arcsin();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.8440 2.3014".to_owned()
        );
        let v3 = C::arcsin_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.5708 -2.2924".to_owned()
        );
        let v3 = v.arccos();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.7268 -2.3014".to_owned()
        );
        let v3 = C::arccos_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.0000 2.2924".to_owned()
        );
        let v3 = v.arctan();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.4186 0.1291".to_owned()
        );
        let v3 = v.arcsec();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.4208 0.1325".to_owned()
        );
        let v3 = C::arcsec_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.3694 0.0000".to_owned()
        );
        let v3 = v.arccsc();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1500 -0.1325".to_owned()
        );
        let v3 = C::arccsc_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.2014 0.0000".to_owned()
        );
        let v3 = v.arccot();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1522 -0.1291".to_owned()
        );
        let v3 = v.sinh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-21.4457 -3.0986".to_owned()
        );
        let v3 = v.cosh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-21.4685 -3.0953".to_owned()
        );
        let v3 = v.tanh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.9990 0.0003".to_owned()
        );
        let v3 = v.sech();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0456 0.0066".to_owned()
        );
        let v3 = v.csch();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0457 0.0066".to_owned()
        );
        let v3 = v.coth();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.0010 -0.0003".to_owned()
        );
        let v3 = v.arcsinh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.3041 0.7070".to_owned()
        );
        let v3 = v.arccosh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.3014 0.7268".to_owned()
        );
        let v3 = C::arccosh_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.2924 0.0000".to_owned()
        );
        let v3 = v.arctanh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1493 1.4372".to_owned()
        );
        let v3 = C::arctanh_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.2027 -1.5708".to_owned()
        );
        let v3 = v.arcsech();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1325 -1.4208".to_owned()
        );
        let v3 = v.arccsch();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1515 -0.1303".to_owned()
        );
        let v3 = v.arccoth();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1493 -0.1336".to_owned()
        );
    }

    #[test]
    #[allow(deprecated)]
    fn complex_f32() {
        type C = num_complex::Complex<f32>;
        use crate::complex::ComplexOps;

        let v = C::rect(10., 10.);
        assert_eq!(v, C { re: 10., im: 10. });
        let v2 = C::rect(1., -1.);
        assert_eq!(v2, C { re: 1., im: -1. });
        let v = C::polar(5., 7.);
        assert_eq!(
            format!("{:.4} {:.4}", v.re, v.im),
            "3.7695 3.2849".to_owned()
        );

        let arg = v.arg();
        assert_eq!(format!("{:.4}", arg), "0.7168".to_owned());
        let arg = v.abs();
        assert_eq!(format!("{:.3}", arg), "5.000".to_owned());
        let arg = v.abs2();
        assert_eq!(format!("{:.3}", arg), "25.000".to_owned());
        let arg = v.logabs();
        assert_eq!(format!("{:.4}", arg), "1.6094".to_owned());

        let v3 = v.add(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "4.7695 2.2849".to_owned()
        );
        let v3 = v.sub(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.7695 4.2849".to_owned()
        );
        let v3 = v.mul(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "7.0544 -0.4846".to_owned()
        );
        let v3 = v.div(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.2423 3.5272".to_owned()
        );
        let v3 = v.add_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "8.7695 3.2849".to_owned()
        );
        let v3 = v.sub_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-1.2305 3.2849".to_owned()
        );
        let v3 = v.mul_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "18.8476 16.4247".to_owned()
        );
        let v3 = v.div_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.7539 0.6570".to_owned()
        );

        let v3 = v.add_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "3.7695 8.2849".to_owned()
        );
        let v3 = v.sub_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "3.7695 -1.7151".to_owned()
        );
        let v3 = v.mul_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-16.4247 18.8476".to_owned()
        );
        let v3 = v.div_imag(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.6570 -0.7539".to_owned()
        );

        let v3 = v.conjugate();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "3.7695 -3.2849".to_owned()
        );
        let v3 = v.inverse();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1508 -0.1314".to_owned()
        );
        let v3 = v.negative();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-3.7695 -3.2849".to_owned()
        );
        let v3 = v.sqrt();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.0940 0.7844".to_owned()
        );
        let v3 = C::sqrt_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.2361 0.0000".to_owned()
        );

        let v3 = v.pow(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "6.4240 -7.9737".to_owned()
        );
        let v3 = v.pow_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-2824.0381 -1338.0712".to_owned()
        );
        let v3 = v.exp();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-42.9142 -6.1938".to_owned()
        );
        let v3 = v.ln();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.6094 0.7168".to_owned()
        );
        let v3 = v.log10();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.6990 0.3113".to_owned()
        );
        let v3 = v.log_b(&v2);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0071 2.0523".to_owned()
        );
        let v3 = v.sin();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-7.8557 -10.7913".to_owned()
        );
        let v3 = v.cos();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-10.8216 7.8337".to_owned()
        );
        let v3 = v.tan();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.0027 0.9991".to_owned()
        );

        let v3 = v.sec();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0606 -0.0439".to_owned()
        );
        let v3 = v.csc();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0441 0.0606".to_owned()
        );
        let v3 = v.cot();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.0027 -1.0009".to_owned()
        );
        let v3 = v.arcsin();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.8440 2.3014".to_owned()
        );
        let v3 = C::arcsin_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.5708 -2.2924".to_owned()
        );
        let v3 = v.arccos();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.7268 -2.3014".to_owned()
        );
        let v3 = C::arccos_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.0000 2.2924".to_owned()
        );
        let v3 = v.arctan();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.4186 0.1291".to_owned()
        );
        let v3 = v.arcsec();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.4208 0.1325".to_owned()
        );
        let v3 = C::arcsec_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.3694 0.0000".to_owned()
        );
        let v3 = v.arccsc();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1500 -0.1325".to_owned()
        );
        let v3 = C::arccsc_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.2014 0.0000".to_owned()
        );
        let v3 = v.arccot();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1522 -0.1291".to_owned()
        );
        let v3 = v.sinh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-21.4457 -3.0986".to_owned()
        );
        let v3 = v.cosh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-21.4685 -3.0953".to_owned()
        );
        let v3 = v.tanh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.9990 0.0003".to_owned()
        );
        let v3 = v.sech();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0456 0.0066".to_owned()
        );
        let v3 = v.csch();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "-0.0457 0.0066".to_owned()
        );
        let v3 = v.coth();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "1.0010 -0.0003".to_owned()
        );
        let v3 = v.arcsinh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.3041 0.7070".to_owned()
        );
        let v3 = v.arccosh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.3014 0.7268".to_owned()
        );
        let v3 = C::arccosh_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "2.2924 0.0000".to_owned()
        );
        let v3 = v.arctanh();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1493 1.4372".to_owned()
        );
        let v3 = C::arctanh_real(5.);
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.2027 -1.5708".to_owned()
        );
        let v3 = v.arcsech();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1325 -1.4208".to_owned()
        );
        let v3 = v.arccsch();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1515 -0.1303".to_owned()
        );
        let v3 = v.arccoth();
        assert_eq!(
            format!("{:.4} {:.4}", v3.re, v3.im),
            "0.1493 -0.1336".to_owned()
        );
    }
}
