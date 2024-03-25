//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// TODO : port to Rust type : http://doc.rust-lang.org/num/complex/struct.Complex.html

use std::fmt::{self, Debug, Formatter};

#[doc(hidden)]
#[allow(clippy::upper_case_acronyms)]
pub trait CFFI<T> {
    fn wrap(s: T) -> Self;
    fn unwrap(self) -> T;
}

#[doc(hidden)]
#[allow(clippy::upper_case_acronyms)]
pub trait FFFI<T> {
    fn wrap(self) -> T;
    fn unwrap(t: T) -> Self;
}

//#[deprecated(note = "Use `Complex64` from the `num_complex` create instead")]
#[repr(C)]
#[derive(Clone, Copy, PartialEq)]
pub struct ComplexF64 {
    pub dat: [f64; 2],
}

impl ComplexF64 {
    /// This function uses the rectangular Cartesian components (x,y) to return the complex number
    /// z = x + i y.
    #[doc(alias = "gsl_complex_rect")]
    pub fn rect(x: f64, y: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_rect(x, y).wrap() }
    }

    /// This function returns the complex number z = r \exp(i \theta) = r (\cos(\theta) + i
    /// \sin(\theta)) from the polar representation (r,theta).
    #[doc(alias = "gsl_complex_polar")]
    pub fn polar(r: f64, theta: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_polar(r, theta).wrap() }
    }

    /// This function returns the argument of the complex number z, \arg(z), where -\pi < \arg(z)
    /// <= \pi.
    #[doc(alias = "gsl_complex_arg")]
    pub fn arg(&self) -> f64 {
        unsafe { sys::gsl_complex_arg(self.unwrap()) }
    }

    /// This function returns the magnitude of the complex number z, |z|.
    #[doc(alias = "gsl_complex_abs")]
    pub fn abs(&self) -> f64 {
        unsafe { sys::gsl_complex_abs(self.unwrap()) }
    }

    /// This function returns the squared magnitude of the complex number z, |z|^2.
    #[doc(alias = "gsl_complex_abs2")]
    pub fn abs2(&self) -> f64 {
        unsafe { sys::gsl_complex_abs2(self.unwrap()) }
    }

    /// This function returns the natural logarithm of the magnitude of the complex number z,
    /// \log|z|.
    ///
    /// It allows an accurate evaluation of \log|z| when |z| is close to one.
    ///
    /// The direct evaluation of log(gsl_complex_abs(z)) would lead to a loss of precision in this
    /// case.
    #[doc(alias = "gsl_complex_logabs")]
    pub fn logabs(&self) -> f64 {
        unsafe { sys::gsl_complex_logabs(self.unwrap()) }
    }

    /// This function returns the sum of the complex numbers a and b, z=a+b.
    #[doc(alias = "gsl_complex_add")]
    pub fn add(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { sys::gsl_complex_add(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the difference of the complex numbers a and b, z=a-b.
    #[doc(alias = "gsl_complex_sub")]
    pub fn sub(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { sys::gsl_complex_sub(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the product of the complex numbers a and b, z=ab.
    #[doc(alias = "gsl_complex_mul")]
    pub fn mul(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { sys::gsl_complex_mul(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the quotient of the complex numbers a and b, z=a/b.
    #[doc(alias = "gsl_complex_div")]
    pub fn div(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { sys::gsl_complex_div(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the sum of the complex number a and the real number x, z=a+x.
    #[doc(alias = "gsl_complex_add_real")]
    pub fn add_real(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_add_real(self.unwrap(), x).wrap() }
    }

    /// This function returns the difference of the complex number a and the real number x, z=a-x.
    #[doc(alias = "gsl_complex_sub_real")]
    pub fn sub_real(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_sub_real(self.unwrap(), x).wrap() }
    }

    /// This function returns the product of the complex number a and the real number x, z=ax.
    #[doc(alias = "gsl_complex_mul_real")]
    pub fn mul_real(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_mul_real(self.unwrap(), x).wrap() }
    }

    /// This function returns the quotient of the complex number a and the real number x, z=a/x.
    #[doc(alias = "gsl_complex_div_real")]
    pub fn div_real(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_div_real(self.unwrap(), x).wrap() }
    }

    /// This function returns the sum of the complex number a and the imaginary number iy, z=a+iy.
    #[doc(alias = "gsl_complex_add_imag")]
    pub fn add_imag(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_add_imag(self.unwrap(), x).wrap() }
    }

    /// This function returns the difference of the complex number a and the imaginary number iy,
    /// z=a-iy.
    #[doc(alias = "gsl_complex_sub_imag")]
    pub fn sub_imag(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_sub_imag(self.unwrap(), x).wrap() }
    }

    /// This function returns the product of the complex number a and the imaginary number iy,
    /// z=a*(iy).
    #[doc(alias = "gsl_complex_mul_imag")]
    pub fn mul_imag(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_mul_imag(self.unwrap(), x).wrap() }
    }

    /// This function returns the quotient of the complex number a and the imaginary number iy,
    /// z=a/(iy).
    #[doc(alias = "gsl_complex_div_imag")]
    pub fn div_imag(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_div_imag(self.unwrap(), x).wrap() }
    }

    /// This function returns the complex conjugate of the complex number z, z^* = x - i y.
    #[doc(alias = "gsl_complex_conjugate")]
    pub fn conjugate(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_conjugate(self.unwrap()).wrap() }
    }

    /// This function returns the inverse, or reciprocal, of the complex number z, 1/z = (x - i y)/
    /// (x^2 + y^2).
    #[doc(alias = "gsl_complex_inverse")]
    pub fn inverse(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_inverse(self.unwrap()).wrap() }
    }

    /// This function returns the negative of the complex number z, -z = (-x) + i(-y).
    #[doc(alias = "gsl_complex_negative")]
    pub fn negative(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_negative(self.unwrap()).wrap() }
    }

    /// This function returns the square root of the complex number z, \sqrt z.
    ///
    /// The branch cut is the negative real axis. The result always lies in the right half of the
    /// omplex plane.
    #[doc(alias = "gsl_complex_sqrt")]
    pub fn sqrt(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_sqrt(self.unwrap()).wrap() }
    }

    /// This function returns the complex square root of the real number x, where x may be negative.
    #[doc(alias = "gsl_complex_sqrt_real")]
    pub fn sqrt_real(x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_sqrt_real(x).wrap() }
    }

    /// The function returns the complex number z raised to the complex power a, z^a.
    /// This is computed as \exp(\log(z)*a) using complex logarithms and complex exponentials.
    #[doc(alias = "gsl_complex_pow")]
    pub fn pow(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { sys::gsl_complex_pow(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the complex number z raised to the real power x, z^x.
    #[doc(alias = "gsl_complex_pow_real")]
    pub fn pow_real(&self, x: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_pow_real(self.unwrap(), x).wrap() }
    }

    /// This function returns the complex exponential of the complex number z, \exp(z).
    #[doc(alias = "gsl_complex_exp")]
    pub fn exp(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_exp(self.unwrap()).wrap() }
    }

    /// This function returns the complex natural logarithm (base e) of the complex number z,
    /// \log(z).
    ///
    /// The branch cut is the negative real axis.
    #[doc(alias = "gsl_complex_log")]
    pub fn log(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_log(self.unwrap()).wrap() }
    }

    /// This function returns the complex base-10 logarithm of the complex number z, \log_10 (z).
    #[doc(alias = "gsl_complex_log10")]
    pub fn log10(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_log10(self.unwrap()).wrap() }
    }

    /// This function returns the complex base-b logarithm of the complex number z, \log_b(z).
    /// This quantity is computed as the ratio \log(z)/\log(b).
    #[doc(alias = "gsl_complex_log_b")]
    pub fn log_b(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { sys::gsl_complex_log_b(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the complex sine of the complex number z, \sin(z) = (\exp(iz) -
    /// \exp(-iz))/(2i).
    #[doc(alias = "gsl_complex_sin")]
    pub fn sin(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_sin(self.unwrap()).wrap() }
    }

    /// This function returns the complex cosine of the complex number z, \cos(z) = (\exp(iz) +
    /// \exp(-iz))/2.
    #[doc(alias = "gsl_complex_cos")]
    pub fn cos(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_cos(self.unwrap()).wrap() }
    }

    /// This function returns the complex tangent of the complex number z, \tan(z) =
    /// \sin(z)/\cos(z).
    #[doc(alias = "gsl_complex_tan")]
    pub fn tan(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_tan(self.unwrap()).wrap() }
    }

    /// This function returns the complex secant of the complex number z, \sec(z) = 1/\cos(z).
    #[doc(alias = "gsl_complex_sec")]
    pub fn sec(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_sec(self.unwrap()).wrap() }
    }

    /// This function returns the complex cosecant of the complex number z, \csc(z) = 1/\sin(z).
    #[doc(alias = "gsl_complex_csc")]
    pub fn csc(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_csc(self.unwrap()).wrap() }
    }

    /// This function returns the complex cotangent of the complex number z, \cot(z) = 1/\tan(z).
    #[doc(alias = "gsl_complex_cot")]
    pub fn cot(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_cot(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsine of the complex number z, \arcsin(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arcsin")]
    pub fn arcsin(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arcsin(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsine of the real number z, \arcsin(z).
    ///
    /// * For z between -1 and 1, the function returns a real value in the range [-\pi/2,\pi/2].
    /// * For z less than -1 the result has a real part of -\pi/2 and a positive imaginary part.
    /// * For z greater than 1 the result has a real part of \pi/2 and a negative imaginary part.
    #[doc(alias = "gsl_complex_arcsin_real")]
    pub fn arcsin_real(z: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_arcsin_real(z).wrap() }
    }

    /// This function returns the complex arccosine of the complex number z, \arccos(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arccos")]
    pub fn arccos(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccos(self.unwrap()).wrap() }
    }

    /// This function returns the complex arccosine of the real number z, \arccos(z).
    ///
    /// * For z between -1 and 1, the function returns a real value in the range [0,\pi].
    /// * For z less than -1 the result has a real part of \pi and a negative imaginary part.
    /// * For z greater than 1 the result is purely imaginary and positive.
    #[doc(alias = "gsl_complex_arccos_real")]
    pub fn arccos_real(z: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccos_real(z).wrap() }
    }

    /// This function returns the complex arctangent of the complex number z, \arctan(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    #[doc(alias = "gsl_complex_arctan")]
    pub fn arctan(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arctan(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsecant of the complex number z, \arcsec(z) =
    /// \arccos(1/z).
    #[doc(alias = "gsl_complex_arcsec")]
    pub fn arcsec(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arcsec(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsecant of the real number z, \arcsec(z) = \arccos(1/z).
    #[doc(alias = "gsl_complex_arcsec_real")]
    pub fn arcsec_real(z: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_arcsec_real(z).wrap() }
    }

    /// This function returns the complex arccosecant of the complex number z, \arccsc(z) =
    /// \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsc")]
    pub fn arccsc(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccsc(self.unwrap()).wrap() }
    }

    /// This function returns the complex arccosecant of the real number z, \arccsc(z) =
    /// \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsc_real")]
    pub fn arccsc_real(z: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccsc_real(z).wrap() }
    }

    /// This function returns the complex arccotangent of the complex number z, \arccot(z) =
    /// \arctan(1/z).
    #[doc(alias = "gsl_complex_arccot")]
    pub fn arccot(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccot(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic sine of the complex number z, \sinh(z) =
    /// (\exp(z) - \exp(-z))/2.
    #[doc(alias = "gsl_complex_sinh")]
    pub fn sinh(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_sinh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic cosine of the complex number z, \cosh(z) =
    /// (\exp(z) + \exp(-z))/2.
    #[doc(alias = "gsl_complex_cosh")]
    pub fn cosh(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_cosh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic tangent of the complex number z, \tanh(z) =
    /// \sinh(z)/\cosh(z).
    #[doc(alias = "gsl_complex_tanh")]
    pub fn tanh(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_tanh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic secant of the complex number z, \sech(z) =
    /// 1/\cosh(z).
    #[doc(alias = "gsl_complex_sech")]
    pub fn sech(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_sech(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic cosecant of the complex number z, \csch(z) =
    /// 1/\sinh(z).
    #[doc(alias = "gsl_complex_csch")]
    pub fn csch(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_csch(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic cotangent of the complex number z, \coth(z) =
    /// 1/\tanh(z).
    #[doc(alias = "gsl_complex_coth")]
    pub fn coth(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_coth(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arcsine of the complex number z, \arcsinh(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    #[doc(alias = "gsl_complex_arcsinh")]
    pub fn arcsinh(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arcsinh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccosine of the complex number z, \arccosh(z).
    /// The branch cut is on the real axis, less than 1.
    /// Note that in this case we use the negative square root in formula 4.6.21 of Abramowitz &
    /// Stegun giving \arccosh(z)=\log(z-\sqrt{z^2-1}).
    #[doc(alias = "gsl_complex_arccosh")]
    pub fn arccosh(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccosh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccosine of the real number z, \arccosh(z).
    #[doc(alias = "gsl_complex_arccosh_real")]
    pub fn arccosh_real(z: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccosh_real(z).wrap() }
    }

    /// This function returns the complex hyperbolic arctangent of the complex number z,
    /// \arctanh(z).
    ///
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arctanh")]
    pub fn arctanh(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arctanh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arctangent of the real number z, \arctanh(z).
    #[doc(alias = "gsl_complex_arctanh_real")]
    pub fn arctanh_real(z: f64) -> ComplexF64 {
        unsafe { sys::gsl_complex_arctanh_real(z).wrap() }
    }

    /// This function returns the complex hyperbolic arcsecant of the complex number z, \arcsech(z)
    /// = \arccosh(1/z).
    #[doc(alias = "gsl_complex_arcsech")]
    pub fn arcsech(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arcsech(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccosecant of the complex number z,
    /// \arccsch(z) = \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsch")]
    pub fn arccsch(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccsch(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccotangent of the complex number z,
    /// \arccoth(z) = \arctanh(1/z).
    #[doc(alias = "gsl_complex_arccoth")]
    pub fn arccoth(&self) -> ComplexF64 {
        unsafe { sys::gsl_complex_arccoth(self.unwrap()).wrap() }
    }

    pub fn real(&self) -> f64 {
        self.dat[0]
    }

    pub fn imaginary(&self) -> f64 {
        self.dat[1]
    }
}

impl Debug for ComplexF64 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "[{}, {}]", self.dat[0], self.dat[1])
    }
}

impl Default for ComplexF64 {
    fn default() -> ComplexF64 {
        ComplexF64 { dat: [0f64, 0f64] }
    }
}

impl CFFI<sys::gsl_complex> for ComplexF64 {
    fn wrap(t: sys::gsl_complex) -> ComplexF64 {
        unsafe { std::mem::transmute(t) }
    }

    fn unwrap(self) -> sys::gsl_complex {
        unsafe { std::mem::transmute(self) }
    }
}

impl CFFI<sys::gsl_complex_float> for ComplexF64 {
    fn wrap(t: sys::gsl_complex_float) -> ComplexF64 {
        ComplexF64 {
            dat: [t.dat[0] as f64, t.dat[1] as f64],
        }
    }

    fn unwrap(self) -> sys::gsl_complex_float {
        sys::gsl_complex_float {
            dat: [self.dat[0] as f32, self.dat[1] as f32],
        }
    }
}

impl FFFI<ComplexF32> for sys::gsl_complex {
    fn wrap(self) -> ComplexF32 {
        ComplexF32 {
            dat: [self.dat[0] as f32, self.dat[1] as f32],
        }
    }

    fn unwrap(t: ComplexF32) -> sys::gsl_complex {
        sys::gsl_complex {
            dat: [t.dat[0] as f64, t.dat[1] as f64],
        }
    }
}

impl FFFI<ComplexF64> for sys::gsl_complex {
    fn wrap(self) -> ComplexF64 {
        unsafe { std::mem::transmute(self) }
    }

    fn unwrap(t: ComplexF64) -> sys::gsl_complex {
        unsafe { std::mem::transmute(t) }
    }
}

#[repr(C)]
#[derive(Clone, Copy, PartialEq)]
pub struct ComplexF32 {
    pub dat: [f32; 2],
}

impl ComplexF32 {
    /// This function uses the rectangular Cartesian components (x,y) to return the complex number
    /// z = x + i y.
    #[doc(alias = "gsl_complex_rect")]
    pub fn rect(x: f32, y: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_rect(x as f64, y as f64).wrap() }
    }

    /// This function returns the complex number z = r \exp(i \theta) = r (\cos(\theta) + i
    /// \sin(\theta)) from the polar representation (r,theta).
    #[doc(alias = "gsl_complex_polar")]
    pub fn polar(r: f32, theta: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_polar(r as f64, theta as f64).wrap() }
    }

    /// This function returns the argument of the complex number z, \arg(z), where -\pi < \arg(z)
    /// <= \pi.
    #[doc(alias = "gsl_complex_arg")]
    pub fn arg(&self) -> f32 {
        unsafe { sys::gsl_complex_arg(self.unwrap()) as f32 }
    }

    /// This function returns the magnitude of the complex number z, |z|.
    #[doc(alias = "gsl_complex_abs")]
    pub fn abs(&self) -> f32 {
        unsafe { sys::gsl_complex_abs(self.unwrap()) as f32 }
    }

    /// This function returns the squared magnitude of the complex number z, |z|^2.
    #[doc(alias = "gsl_complex_abs2")]
    pub fn abs2(&self) -> f32 {
        unsafe { sys::gsl_complex_abs2(self.unwrap()) as f32 }
    }

    /// This function returns the natural logarithm of the magnitude of the complex number z,
    /// \log|z|.
    ///
    /// It allows an accurate evaluation of \log|z| when |z| is close to one.
    /// The direct evaluation of log(gsl_complex_abs(z)) would lead to a loss of precision in
    /// this case.
    #[doc(alias = "gsl_complex_logabs")]
    pub fn logabs(&self) -> f32 {
        unsafe { sys::gsl_complex_logabs(self.unwrap()) as f32 }
    }

    /// This function returns the sum of the complex numbers a and b, z=a+b.
    #[doc(alias = "gsl_complex_add")]
    pub fn add(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { sys::gsl_complex_add(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the difference of the complex numbers a and b, z=a-b.
    #[doc(alias = "gsl_complex_sub")]
    pub fn sub(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { sys::gsl_complex_sub(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the product of the complex numbers a and b, z=ab.
    #[doc(alias = "gsl_complex_mul")]
    pub fn mul(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { sys::gsl_complex_mul(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the quotient of the complex numbers a and b, z=a/b.
    #[doc(alias = "gsl_complex_div")]
    pub fn div(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { sys::gsl_complex_div(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the sum of the complex number a and the real number x, z=a+x.
    #[doc(alias = "gsl_complex_add_real")]
    pub fn add_real(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_add_real(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the difference of the complex number a and the real number x, z=a-x.
    #[doc(alias = "gsl_complex_sub_real")]
    pub fn sub_real(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_sub_real(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the product of the complex number a and the real number x, z=ax.
    #[doc(alias = "gsl_complex_mul_real")]
    pub fn mul_real(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_mul_real(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the quotient of the complex number a and the real number x, z=a/x.
    #[doc(alias = "gsl_complex_div_real")]
    pub fn div_real(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_div_real(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the sum of the complex number a and the imaginary number iy, z=a+iy.
    #[doc(alias = "gsl_complex_add_imag")]
    pub fn add_imag(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_add_imag(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the difference of the complex number a and the imaginary number iy, z=a-iy.
    #[doc(alias = "gsl_complex_sub_imag")]
    pub fn sub_imag(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_sub_imag(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the product of the complex number a and the imaginary number iy, z=a*(iy).
    #[doc(alias = "gsl_complex_mul_imag")]
    pub fn mul_imag(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_mul_imag(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the quotient of the complex number a and the imaginary number iy, z=a/(iy).
    #[doc(alias = "gsl_complex_div_imag")]
    pub fn div_imag(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_div_imag(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the complex conjugate of the complex number z, z^* = x - i y.
    #[doc(alias = "gsl_complex_conjugate")]
    pub fn conjugate(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_conjugate(self.unwrap()).wrap() }
    }

    /// This function returns the inverse, or reciprocal, of the complex number z, 1/z = (x - i y)/
    /// (x^2 + y^2).
    #[doc(alias = "gsl_complex_inverse")]
    pub fn inverse(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_inverse(self.unwrap()).wrap() }
    }

    /// This function returns the negative of the complex number z, -z = (-x) + i(-y).
    #[doc(alias = "gsl_complex_negative")]
    pub fn negative(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_negative(self.unwrap()).wrap() }
    }

    /// This function returns the square root of the complex number z, \sqrt z.
    ///
    /// The branch cut is the negative real axis. The result always lies in the right half of the
    /// complex plane.
    #[doc(alias = "gsl_complex_sqrt")]
    pub fn sqrt(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_sqrt(self.unwrap()).wrap() }
    }

    /// This function returns the complex square root of the real number x, where x may be negative.
    #[doc(alias = "gsl_complex_sqrt_real")]
    pub fn sqrt_real(x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_sqrt_real(x as f64).wrap() }
    }

    /// The function returns the complex number z raised to the complex power a, z^a.
    ///
    /// This is computed as \exp(\log(z)*a) using complex logarithms and complex exponentials.
    #[doc(alias = "gsl_complex_pow")]
    pub fn pow(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { sys::gsl_complex_pow(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the complex number z raised to the real power x, z^x.
    #[doc(alias = "gsl_complex_pow_real")]
    pub fn pow_real(&self, x: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_pow_real(self.unwrap(), x as f64).wrap() }
    }

    /// This function returns the complex exponential of the complex number z, \exp(z).
    #[doc(alias = "gsl_complex_exp")]
    pub fn exp(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_exp(self.unwrap()).wrap() }
    }

    /// This function returns the complex natural logarithm (base e) of the complex number z, \log(z).
    /// The branch cut is the negative real axis.
    #[doc(alias = "gsl_complex_log")]
    pub fn log(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_log(self.unwrap()).wrap() }
    }

    /// This function returns the complex base-10 logarithm of the complex number z, \log_10 (z).
    #[doc(alias = "gsl_complex_log10")]
    pub fn log10(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_log10(self.unwrap()).wrap() }
    }

    /// This function returns the complex base-b logarithm of the complex number z, \log_b(z).
    /// This quantity is computed as the ratio \log(z)/\log(b).
    #[doc(alias = "gsl_complex_log_b")]
    pub fn log_b(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { sys::gsl_complex_log_b(self.unwrap(), other.unwrap()).wrap() }
    }

    /// This function returns the complex sine of the complex number z, \sin(z) = (\exp(iz) -
    /// \exp(-iz))/(2i).
    #[doc(alias = "gsl_complex_sin")]
    pub fn sin(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_sin(self.unwrap()).wrap() }
    }

    /// This function returns the complex cosine of the complex number z, \cos(z) = (\exp(iz) +
    /// \exp(-iz))/2.
    #[doc(alias = "gsl_complex_cos")]
    pub fn cos(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_cos(self.unwrap()).wrap() }
    }

    /// This function returns the complex tangent of the complex number z, \tan(z) =
    /// \sin(z)/\cos(z).
    #[doc(alias = "gsl_complex_tan")]
    pub fn tan(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_tan(self.unwrap()).wrap() }
    }

    /// This function returns the complex secant of the complex number z, \sec(z) = 1/\cos(z).
    #[doc(alias = "gsl_complex_sec")]
    pub fn sec(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_sec(self.unwrap()).wrap() }
    }

    /// This function returns the complex cosecant of the complex number z, \csc(z) = 1/\sin(z).
    #[doc(alias = "gsl_complex_csc")]
    pub fn csc(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_csc(self.unwrap()).wrap() }
    }

    /// This function returns the complex cotangent of the complex number z, \cot(z) = 1/\tan(z).
    #[doc(alias = "gsl_complex_cot")]
    pub fn cot(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_cot(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsine of the complex number z, \arcsin(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arcsin")]
    pub fn arcsin(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arcsin(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsine of the real number z, \arcsin(z).
    ///
    /// * For z between -1 and 1, the function returns a real value in the range [-\pi/2,\pi/2].
    /// * For z less than -1 the result has a real part of -\pi/2 and a positive imaginary part.
    /// * For z greater than 1 the result has a real part of \pi/2 and a negative imaginary part.
    #[doc(alias = "gsl_complex_arcsin_real")]
    pub fn arcsin_real(z: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_arcsin_real(z as f64).wrap() }
    }

    /// This function returns the complex arccosine of the complex number z, \arccos(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arccos")]
    pub fn arccos(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccos(self.unwrap()).wrap() }
    }

    /// This function returns the complex arccosine of the real number z, \arccos(z).
    ///
    /// * For z between -1 and 1, the function returns a real value in the range [0,\pi].
    /// * For z less than -1 the result has a real part of \pi and a negative imaginary part.
    /// * For z greater than 1 the result is purely imaginary and positive.
    #[doc(alias = "gsl_complex_arccos_real")]
    pub fn arccos_real(z: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccos_real(z as f64).wrap() }
    }

    /// This function returns the complex arctangent of the complex number z, \arctan(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    #[doc(alias = "gsl_complex_arctan")]
    pub fn arctan(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arctan(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsecant of the complex number z, \arcsec(z) =
    /// \arccos(1/z).
    #[doc(alias = "gsl_complex_arcsec")]
    pub fn arcsec(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arcsec(self.unwrap()).wrap() }
    }

    /// This function returns the complex arcsecant of the real number z, \arcsec(z) = \arccos(1/z).
    #[doc(alias = "gsl_complex_arcsec_real")]
    pub fn arcsec_real(z: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_arcsec_real(z as f64).wrap() }
    }

    /// This function returns the complex arccosecant of the complex number z, \arccsc(z) =
    /// \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsc")]
    pub fn arccsc(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccsc(self.unwrap()).wrap() }
    }

    /// This function returns the complex arccosecant of the real number z, \arccsc(z) =
    /// \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsc_real")]
    pub fn arccsc_real(z: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccsc_real(z as f64).wrap() }
    }

    /// This function returns the complex arccotangent of the complex number z, \arccot(z) =
    /// \arctan(1/z).
    #[doc(alias = "gsl_complex_arccot")]
    pub fn arccot(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccot(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic sine of the complex number z, \sinh(z) =
    /// (\exp(z) - \exp(-z))/2.
    #[doc(alias = "gsl_complex_sinh")]
    pub fn sinh(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_sinh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic cosine of the complex number z, \cosh(z) =
    /// (\exp(z) + \exp(-z))/2.
    #[doc(alias = "gsl_complex_cosh")]
    pub fn cosh(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_cosh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic tangent of the complex number z, \tanh(z) =
    /// \sinh(z)/\cosh(z).
    #[doc(alias = "gsl_complex_tanh")]
    pub fn tanh(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_tanh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic secant of the complex number z, \sech(z) =
    /// 1/\cosh(z).
    #[doc(alias = "gsl_complex_sech")]
    pub fn sech(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_sech(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic cosecant of the complex number z, \csch(z) =
    /// 1/\sinh(z).
    #[doc(alias = "gsl_complex_csch")]
    pub fn csch(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_csch(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic cotangent of the complex number z, \coth(z) =
    /// 1/\tanh(z).
    #[doc(alias = "gsl_complex_coth")]
    pub fn coth(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_coth(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arcsine of the complex number z, \arcsinh(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    #[doc(alias = "gsl_complex_arcsinh")]
    pub fn arcsinh(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arcsinh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccosine of the complex number z, \arccosh(z).
    ///
    /// The branch cut is on the real axis, less than 1.
    ///
    /// Note that in this case we use the negative square root in formula 4.6.21 of Abramowitz &
    /// Stegun giving \arccosh(z)=\log(z-\sqrt{z^2-1}).
    #[doc(alias = "gsl_complex_arccosh")]
    pub fn arccosh(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccosh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccosine of the real number z, \arccosh(z).
    #[doc(alias = "gsl_complex_arccosh_real")]
    pub fn arccosh_real(z: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccosh_real(z as f64).wrap() }
    }

    /// This function returns the complex hyperbolic arctangent of the complex number z,
    /// arctanh(z).
    ///
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    #[doc(alias = "gsl_complex_arctanh")]
    pub fn arctanh(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arctanh(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arctangent of the real number z, \arctanh(z).
    #[doc(alias = "gsl_complex_arctanh_real")]
    pub fn arctanh_real(z: f32) -> ComplexF32 {
        unsafe { sys::gsl_complex_arctanh_real(z as f64).wrap() }
    }

    /// This function returns the complex hyperbolic arcsecant of the complex number z, \arcsech(z)
    /// = \arccosh(1/z).
    #[doc(alias = "gsl_complex_arcsech")]
    pub fn arcsech(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arcsech(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccosecant of the complex number z,
    /// \arccsch(z) = \arcsin(1/z).
    #[doc(alias = "gsl_complex_arccsch")]
    pub fn arccsch(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccsch(self.unwrap()).wrap() }
    }

    /// This function returns the complex hyperbolic arccotangent of the complex number z,
    /// \arccoth(z) = \arctanh(1/z).
    #[doc(alias = "gsl_complex_arccoth")]
    pub fn arccoth(&self) -> ComplexF32 {
        unsafe { sys::gsl_complex_arccoth(self.unwrap()).wrap() }
    }

    pub fn real(&self) -> f32 {
        self.dat[0]
    }

    pub fn imaginary(&self) -> f32 {
        self.dat[1]
    }
}

impl Debug for ComplexF32 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "[{}, {}]", self.dat[0], self.dat[1])
    }
}

impl Default for ComplexF32 {
    fn default() -> ComplexF32 {
        ComplexF32 { dat: [0f32, 0f32] }
    }
}

impl CFFI<sys::gsl_complex> for ComplexF32 {
    fn wrap(s: sys::gsl_complex) -> ComplexF32 {
        ComplexF32 {
            dat: [s.dat[0] as f32, s.dat[1] as f32],
        }
    }

    fn unwrap(self) -> sys::gsl_complex {
        sys::gsl_complex {
            dat: [self.dat[0] as f64, self.dat[1] as f64],
        }
    }
}

impl CFFI<sys::gsl_complex_float> for ComplexF32 {
    fn wrap(s: sys::gsl_complex_float) -> ComplexF32 {
        unsafe { std::mem::transmute(s) }
    }

    fn unwrap(self) -> sys::gsl_complex_float {
        unsafe { std::mem::transmute(self) }
    }
}

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
fn complex_f64() {
    let v = ComplexF64::rect(10., 10.);
    assert_eq!(v, ComplexF64 { dat: [10., 10.] });
    let v2 = ComplexF64::rect(1., -1.);
    assert_eq!(v2, ComplexF64 { dat: [1., -1.] });
    let v = ComplexF64::polar(5., 7.);
    assert_eq!(
        format!("{:.4} {:.4}", v.dat[0], v.dat[1]),
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
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "4.7695 2.2849".to_owned()
    );
    let v3 = v.sub(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.7695 4.2849".to_owned()
    );
    let v3 = v.mul(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "7.0544 -0.4846".to_owned()
    );
    let v3 = v.div(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.2423 3.5272".to_owned()
    );
    let v3 = v.add_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "8.7695 3.2849".to_owned()
    );
    let v3 = v.sub_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-1.2305 3.2849".to_owned()
    );
    let v3 = v.mul_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "18.8476 16.4247".to_owned()
    );
    let v3 = v.div_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.7539 0.6570".to_owned()
    );

    let v3 = v.add_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "3.7695 8.2849".to_owned()
    );
    let v3 = v.sub_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "3.7695 -1.7151".to_owned()
    );
    let v3 = v.mul_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-16.4247 18.8476".to_owned()
    );
    let v3 = v.div_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.6570 -0.7539".to_owned()
    );

    let v3 = v.conjugate();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "3.7695 -3.2849".to_owned()
    );
    let v3 = v.inverse();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1508 -0.1314".to_owned()
    );
    let v3 = v.negative();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-3.7695 -3.2849".to_owned()
    );
    let v3 = v.sqrt();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.0940 0.7844".to_owned()
    );
    let v3 = ComplexF64::sqrt_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.2361 0.0000".to_owned()
    );

    let v3 = v.pow(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "6.4240 -7.9737".to_owned()
    );
    let v3 = v.pow_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-2824.0381 -1338.0708".to_owned()
    );
    let v3 = v.exp();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-42.9142 -6.1938".to_owned()
    );
    let v3 = v.log();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.6094 0.7168".to_owned()
    );
    let v3 = v.log10();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.6990 0.3113".to_owned()
    );
    let v3 = v.log_b(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0071 2.0523".to_owned()
    );
    let v3 = v.sin();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-7.8557 -10.7913".to_owned()
    );
    let v3 = v.cos();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-10.8216 7.8337".to_owned()
    );
    let v3 = v.tan();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.0027 0.9991".to_owned()
    );

    let v3 = v.sec();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0606 -0.0439".to_owned()
    );
    let v3 = v.csc();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0441 0.0606".to_owned()
    );
    let v3 = v.cot();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.0027 -1.0009".to_owned()
    );
    let v3 = v.arcsin();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.8440 2.3014".to_owned()
    );
    let v3 = ComplexF64::arcsin_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.5708 -2.2924".to_owned()
    );
    let v3 = v.arccos();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.7268 -2.3014".to_owned()
    );
    let v3 = ComplexF64::arccos_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.0000 2.2924".to_owned()
    );
    let v3 = v.arctan();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.4186 0.1291".to_owned()
    );
    let v3 = v.arcsec();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.4208 0.1325".to_owned()
    );
    let v3 = ComplexF64::arcsec_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.3694 0.0000".to_owned()
    );
    let v3 = v.arccsc();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1500 -0.1325".to_owned()
    );
    let v3 = ComplexF64::arccsc_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.2014 0.0000".to_owned()
    );
    let v3 = v.arccot();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1522 -0.1291".to_owned()
    );
    let v3 = v.sinh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-21.4457 -3.0986".to_owned()
    );
    let v3 = v.cosh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-21.4685 -3.0953".to_owned()
    );
    let v3 = v.tanh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.9990 0.0003".to_owned()
    );
    let v3 = v.sech();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0456 0.0066".to_owned()
    );
    let v3 = v.csch();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0457 0.0066".to_owned()
    );
    let v3 = v.coth();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.0010 -0.0003".to_owned()
    );
    let v3 = v.arcsinh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.3041 0.7070".to_owned()
    );
    let v3 = v.arccosh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.3014 0.7268".to_owned()
    );
    let v3 = ComplexF64::arccosh_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.2924 0.0000".to_owned()
    );
    let v3 = v.arctanh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1493 1.4372".to_owned()
    );
    let v3 = ComplexF64::arctanh_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.2027 -1.5708".to_owned()
    );
    let v3 = v.arcsech();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1325 -1.4208".to_owned()
    );
    let v3 = v.arccsch();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1515 -0.1303".to_owned()
    );
    let v3 = v.arccoth();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1493 -0.1336".to_owned()
    );
}

#[test]
fn complex_f32() {
    let v = ComplexF32::rect(10., 10.);
    assert_eq!(v, ComplexF32 { dat: [10., 10.] });
    let v2 = ComplexF32::rect(1., -1.);
    assert_eq!(v2, ComplexF32 { dat: [1., -1.] });
    let v = ComplexF32::polar(5., 7.);
    assert_eq!(
        format!("{:.4} {:.4}", v.dat[0], v.dat[1]),
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
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "4.7695 2.2849".to_owned()
    );
    let v3 = v.sub(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.7695 4.2849".to_owned()
    );
    let v3 = v.mul(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "7.0544 -0.4846".to_owned()
    );
    let v3 = v.div(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.2423 3.5272".to_owned()
    );
    let v3 = v.add_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "8.7695 3.2849".to_owned()
    );
    let v3 = v.sub_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-1.2305 3.2849".to_owned()
    );
    let v3 = v.mul_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "18.8476 16.4247".to_owned()
    );
    let v3 = v.div_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.7539 0.6570".to_owned()
    );

    let v3 = v.add_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "3.7695 8.2849".to_owned()
    );
    let v3 = v.sub_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "3.7695 -1.7151".to_owned()
    );
    let v3 = v.mul_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-16.4247 18.8476".to_owned()
    );
    let v3 = v.div_imag(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.6570 -0.7539".to_owned()
    );

    let v3 = v.conjugate();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "3.7695 -3.2849".to_owned()
    );
    let v3 = v.inverse();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1508 -0.1314".to_owned()
    );
    let v3 = v.negative();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-3.7695 -3.2849".to_owned()
    );
    let v3 = v.sqrt();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.0940 0.7844".to_owned()
    );
    let v3 = ComplexF32::sqrt_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.2361 0.0000".to_owned()
    );

    let v3 = v.pow(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "6.4240 -7.9737".to_owned()
    );
    let v3 = v.pow_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-2824.0381 -1338.0712".to_owned()
    );
    let v3 = v.exp();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-42.9142 -6.1938".to_owned()
    );
    let v3 = v.log();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.6094 0.7168".to_owned()
    );
    let v3 = v.log10();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.6990 0.3113".to_owned()
    );
    let v3 = v.log_b(&v2);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0071 2.0523".to_owned()
    );
    let v3 = v.sin();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-7.8557 -10.7913".to_owned()
    );
    let v3 = v.cos();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-10.8216 7.8337".to_owned()
    );
    let v3 = v.tan();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.0027 0.9991".to_owned()
    );

    let v3 = v.sec();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0606 -0.0439".to_owned()
    );
    let v3 = v.csc();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0441 0.0606".to_owned()
    );
    let v3 = v.cot();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.0027 -1.0009".to_owned()
    );
    let v3 = v.arcsin();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.8440 2.3014".to_owned()
    );
    let v3 = ComplexF32::arcsin_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.5708 -2.2924".to_owned()
    );
    let v3 = v.arccos();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.7268 -2.3014".to_owned()
    );
    let v3 = ComplexF32::arccos_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.0000 2.2924".to_owned()
    );
    let v3 = v.arctan();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.4186 0.1291".to_owned()
    );
    let v3 = v.arcsec();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.4208 0.1325".to_owned()
    );
    let v3 = ComplexF32::arcsec_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.3694 0.0000".to_owned()
    );
    let v3 = v.arccsc();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1500 -0.1325".to_owned()
    );
    let v3 = ComplexF32::arccsc_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.2014 0.0000".to_owned()
    );
    let v3 = v.arccot();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1522 -0.1291".to_owned()
    );
    let v3 = v.sinh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-21.4457 -3.0986".to_owned()
    );
    let v3 = v.cosh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-21.4685 -3.0953".to_owned()
    );
    let v3 = v.tanh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.9990 0.0003".to_owned()
    );
    let v3 = v.sech();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0456 0.0066".to_owned()
    );
    let v3 = v.csch();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "-0.0457 0.0066".to_owned()
    );
    let v3 = v.coth();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "1.0010 -0.0003".to_owned()
    );
    let v3 = v.arcsinh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.3041 0.7070".to_owned()
    );
    let v3 = v.arccosh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.3014 0.7268".to_owned()
    );
    let v3 = ComplexF32::arccosh_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "2.2924 0.0000".to_owned()
    );
    let v3 = v.arctanh();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1493 1.4372".to_owned()
    );
    let v3 = ComplexF32::arctanh_real(5.);
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.2027 -1.5708".to_owned()
    );
    let v3 = v.arcsech();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1325 -1.4208".to_owned()
    );
    let v3 = v.arccsch();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1515 -0.1303".to_owned()
    );
    let v3 = v.arccoth();
    assert_eq!(
        format!("{:.4} {:.4}", v3.dat[0], v3.dat[1]),
        "0.1493 -0.1336".to_owned()
    );
}
