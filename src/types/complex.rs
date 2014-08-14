//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// TODO : port to Rust type : http://doc.rust-lang.org/num/complex/struct.Complex.html

use std::fmt::{Formatter,Show};
use std::fmt;
use std::default::Default;

pub struct ComplexF64 {
    pub data: [f64, ..2]
}

impl ComplexF64 {
    /// This function uses the rectangular Cartesian components (x,y) to return the complex number z = x + i y.
    pub fn rect(x: f64, y: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_rect(x, y)) }
    }

    /// This function returns the complex number z = r \exp(i \theta) = r (\cos(\theta) + i \sin(\theta)) from the polar representation (r,theta).
    pub fn polar(r: f64, theta: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_polar(r, theta)) }
    }

    /// This function returns the argument of the complex number z, \arg(z), where -\pi < \arg(z) <= \pi.
    pub fn arg(&self) -> f64 {
        unsafe { ::ffi::gsl_complex_arg(::std::mem::transmute(self.data)) }
    }

    /// This function returns the magnitude of the complex number z, |z|.
    pub fn abs(&self) -> f64 {
        unsafe { ::ffi::gsl_complex_abs(::std::mem::transmute(self.data)) }
    }

    /// This function returns the squared magnitude of the complex number z, |z|^2.
    pub fn abs2(&self) -> f64 {
        unsafe { ::ffi::gsl_complex_abs2(::std::mem::transmute(self.data)) }
    }

    /// This function returns the natural logarithm of the magnitude of the complex number z, \log|z|.
    /// It allows an accurate evaluation of \log|z| when |z| is close to one.
    /// The direct evaluation of log(gsl_complex_abs(z)) would lead to a loss of precision in this case.
    pub fn logabs(&self) -> f64 {
        unsafe { ::ffi::gsl_complex_logabs(::std::mem::transmute(self.data)) }
    }

    /// This function returns the sum of the complex numbers a and b, z=a+b.
    pub fn add(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_add(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the difference of the complex numbers a and b, z=a-b.
    pub fn sub(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sub(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the product of the complex numbers a and b, z=ab.
    pub fn mul(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_mul(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the quotient of the complex numbers a and b, z=a/b.
    pub fn div(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_div(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the sum of the complex number a and the real number x, z=a+x.
    pub fn add_real(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_add_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the difference of the complex number a and the real number x, z=a-x.
    pub fn sub_real(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sub_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the product of the complex number a and the real number x, z=ax.
    pub fn mul_real(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_mul_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the quotient of the complex number a and the real number x, z=a/x.
    pub fn div_real(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_div_real(::std::mem::transmute(self.data), x)) }
    }
    
    /// This function returns the sum of the complex number a and the imaginary number iy, z=a+iy.
    pub fn add_imag(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_add_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the difference of the complex number a and the imaginary number iy, z=a-iy.
    pub fn sub_imag(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sub_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the product of the complex number a and the imaginary number iy, z=a*(iy).
    pub fn mul_imag(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_mul_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the quotient of the complex number a and the imaginary number iy, z=a/(iy).
    pub fn div_imag(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_div_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the complex conjugate of the complex number z, z^* = x - i y.
    pub fn conjugate(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_conjugate(::std::mem::transmute(self.data))) }
    }

    /// This function returns the inverse, or reciprocal, of the complex number z, 1/z = (x - i y)/(x^2 + y^2).
    pub fn inverse(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_inverse(::std::mem::transmute(self.data))) }
    }

    /// This function returns the negative of the complex number z, -z = (-x) + i(-y).
    pub fn negative(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_negative(::std::mem::transmute(self.data))) }
    }

    /// This function returns the square root of the complex number z, \sqrt z.
    /// The branch cut is the negative real axis. The result always lies in the right half of the complex plane.
    pub fn sqrt(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sqrt(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex square root of the real number x, where x may be negative.
    pub fn sqrt_real(x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sqrt_real(x)) }
    }

    /// The function returns the complex number z raised to the complex power a, z^a.
    /// This is computed as \exp(\log(z)*a) using complex logarithms and complex exponentials.
    pub fn pow(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_pow(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the complex number z raised to the real power x, z^x.
    pub fn pow_real(&self, x: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_pow_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the complex exponential of the complex number z, \exp(z).
    pub fn exp(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_exp(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex natural logarithm (base e) of the complex number z, \log(z).
    /// The branch cut is the negative real axis.
    pub fn log(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_log(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex base-10 logarithm of the complex number z, \log_10 (z).
    pub fn log10(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_log10(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex base-b logarithm of the complex number z, \log_b(z).
    /// This quantity is computed as the ratio \log(z)/\log(b).
    pub fn log_b(&self, other: &ComplexF64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_log_b(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the complex sine of the complex number z, \sin(z) = (\exp(iz) - \exp(-iz))/(2i).
    pub fn sin(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sin(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex cosine of the complex number z, \cos(z) = (\exp(iz) + \exp(-iz))/2.
    pub fn cos(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_cos(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex tangent of the complex number z, \tan(z) = \sin(z)/\cos(z).
    pub fn tan(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_tan(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex secant of the complex number z, \sec(z) = 1/\cos(z).
    pub fn sec(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sec(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex cosecant of the complex number z, \csc(z) = 1/\sin(z).
    pub fn csc(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_csc(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex cotangent of the complex number z, \cot(z) = 1/\tan(z).
    pub fn cot(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_cot(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsine of the complex number z, \arcsin(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    pub fn arcsin(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arcsin(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsine of the real number z, \arcsin(z).
    /// 
    /// * For z between -1 and 1, the function returns a real value in the range [-\pi/2,\pi/2].
    /// * For z less than -1 the result has a real part of -\pi/2 and a positive imaginary part.
    /// * For z greater than 1 the result has a real part of \pi/2 and a negative imaginary part.
    pub fn arcsin_real(z: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arcsin_real(z)) }
    }

    /// This function returns the complex arccosine of the complex number z, \arccos(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    pub fn arccos(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccos(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arccosine of the real number z, \arccos(z).
    /// 
    /// * For z between -1 and 1, the function returns a real value in the range [0,\pi].
    /// * For z less than -1 the result has a real part of \pi and a negative imaginary part.
    /// * For z greater than 1 the result is purely imaginary and positive.
    pub fn arccos_real(z: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccos_real(z)) }
    }

    /// This function returns the complex arctangent of the complex number z, \arctan(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    pub fn arctan(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arctan(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsecant of the complex number z, \arcsec(z) = \arccos(1/z).
    pub fn arcsec(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arcsec(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsecant of the real number z, \arcsec(z) = \arccos(1/z).
    pub fn arcsec_real(z: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arcsec_real(z)) }
    }

    /// This function returns the complex arccosecant of the complex number z, \arccsc(z) = \arcsin(1/z).
    pub fn arccsc(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccsc(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arccosecant of the real number z, \arccsc(z) = \arcsin(1/z).
    pub fn arccsc_real(z: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccsc_real(z)) }
    }

    /// This function returns the complex arccotangent of the complex number z, \arccot(z) = \arctan(1/z).
    pub fn arccot(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccot(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic sine of the complex number z, \sinh(z) = (\exp(z) - \exp(-z))/2.
    pub fn sinh(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sinh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic cosine of the complex number z, \cosh(z) = (\exp(z) + \exp(-z))/2.
    pub fn cosh(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_cosh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic tangent of the complex number z, \tanh(z) = \sinh(z)/\cosh(z).
    pub fn tanh(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_tanh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic secant of the complex number z, \sech(z) = 1/\cosh(z).
    pub fn sech(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_sech(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic cosecant of the complex number z, \csch(z) = 1/\sinh(z).
    pub fn csch(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_csch(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic cotangent of the complex number z, \coth(z) = 1/\tanh(z).
    pub fn coth(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_coth(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arcsine of the complex number z, \arcsinh(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    pub fn arcsinh(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arcsinh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccosine of the complex number z, \arccosh(z).
    /// The branch cut is on the real axis, less than 1.
    /// Note that in this case we use the negative square root in formula 4.6.21 of Abramowitz & Stegun giving \arccosh(z)=\log(z-\sqrt{z^2-1}).
    pub fn arccosh(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccosh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccosine of the real number z, \arccosh(z).
    pub fn arccosh_real(z: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccosh_real(z)) }
    }

    /// This function returns the complex hyperbolic arctangent of the complex number z, \arctanh(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    pub fn arctanh(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arctanh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arctangent of the real number z, \arctanh(z).
    pub fn arctanh_real(z: f64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arctanh_real(z)) }
    }

    /// This function returns the complex hyperbolic arcsecant of the complex number z, \arcsech(z) = \arccosh(1/z).
    pub fn arcsech(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arcsech(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccosecant of the complex number z, \arccsch(z) = \arcsin(1/z).
    pub fn arccsch(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccsch(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccotangent of the complex number z, \arccoth(z) = \arctanh(1/z).
    pub fn arccoth(&self) -> ComplexF64 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_arccoth(::std::mem::transmute(self.data))) }
    }
}

impl Show for ComplexF64 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "[{}, {}]", self.data[0], self.data[1])
    }
}

impl Clone for ComplexF64 {
    fn clone(&self) -> ComplexF64 {
        ComplexF64 {
            data: [self.data[0], self.data[1]]
        }
    }
}

impl Default for ComplexF64 {
    fn default() -> ComplexF64 {
        ComplexF64 {
            data: [0f64, 0f64]
        }
    }
}

pub struct ComplexF32 {
    pub data: [f32, ..2]
}

// I'll implement it in Rust directly
/*impl ComplexF32 {
    /// This function returns the argument of the complex number z, \arg(z), where -\pi < \arg(z) <= \pi.
    pub fn arg(&self) -> f32 {
        unsafe { ::ffi::gsl_complex_float_arg(::std::mem::transmute(self.data)) }
    }

    /// This function returns the magnitude of the complex number z, |z|.
    pub fn abs(&self) -> f32 {
        unsafe { ::ffi::gsl_complex_float_abs(::std::mem::transmute(self.data)) }
    }

    /// This function returns the squared magnitude of the complex number z, |z|^2.
    pub fn abs2(&self) -> f32 {
        unsafe { ::ffi::gsl_complex_float_abs2(::std::mem::transmute(self.data)) }
    }

    /// This function returns the natural logarithm of the magnitude of the complex number z, \log|z|.
    /// It allows an accurate evaluation of \log|z| when |z| is close to one.
    /// The direct evaluation of log(gsl_complex_float_abs(z)) would lead to a loss of precision in this case.
    pub fn logabs(&self) -> f32 {
        unsafe { ::ffi::gsl_complex_float_logabs(::std::mem::transmute(self.data)) }
    }

    /// This function returns the sum of the complex numbers a and b, z=a+b.
    pub fn add(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_add(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the difference of the complex numbers a and b, z=a-b.
    pub fn sub(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sub(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the product of the complex numbers a and b, z=ab.
    pub fn mul(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_mul(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the quotient of the complex numbers a and b, z=a/b.
    pub fn div(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_div(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the sum of the complex number a and the real number x, z=a+x.
    pub fn add_real(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_add_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the difference of the complex number a and the real number x, z=a-x.
    pub fn sub_real(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sub_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the product of the complex number a and the real number x, z=ax.
    pub fn mul_real(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_mul_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the quotient of the complex number a and the real number x, z=a/x.
    pub fn div_real(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_div_real(::std::mem::transmute(self.data), x)) }
    }
    
    /// This function returns the sum of the complex number a and the imaginary number iy, z=a+iy.
    pub fn add_imag(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_add_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the difference of the complex number a and the imaginary number iy, z=a-iy.
    pub fn sub_imag(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sub_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the product of the complex number a and the imaginary number iy, z=a*(iy).
    pub fn mul_imag(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_mul_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the quotient of the complex number a and the imaginary number iy, z=a/(iy).
    pub fn div_imag(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_div_imag(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the complex conjugate of the complex number z, z^* = x - i y.
    pub fn conjugate(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_conjugate(::std::mem::transmute(self.data))) }
    }

    /// This function returns the inverse, or reciprocal, of the complex number z, 1/z = (x - i y)/(x^2 + y^2).
    pub fn inverse(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_inverse(::std::mem::transmute(self.data))) }
    }

    /// This function returns the negative of the complex number z, -z = (-x) + i(-y).
    pub fn negative(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_negative(::std::mem::transmute(self.data))) }
    }

    /// This function returns the square root of the complex number z, \sqrt z.
    /// The branch cut is the negative real axis. The result always lies in the right half of the complex plane.
    pub fn sqrt(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sqrt(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex square root of the real number x, where x may be negative.
    pub fn sqrt_real(x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sqrt_real(x)) }
    }

    /// The function returns the complex number z raised to the complex power a, z^a.
    /// This is computed as \exp(\log(z)*a) using complex logarithms and complex exponentials.
    pub fn pow(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_pow(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the complex number z raised to the real power x, z^x.
    pub fn pow_real(&self, x: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_pow_real(::std::mem::transmute(self.data), x)) }
    }

    /// This function returns the complex exponential of the complex number z, \exp(z).
    pub fn exp(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_exp(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex natural logarithm (base e) of the complex number z, \log(z).
    /// The branch cut is the negative real axis.
    pub fn log(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_log(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex base-10 logarithm of the complex number z, \log_10 (z).
    pub fn log10(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_log10(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex base-b logarithm of the complex number z, \log_b(z).
    /// This quantity is computed as the ratio \log(z)/\log(b).
    pub fn log_b(&self, other: &ComplexF32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_log_b(::std::mem::transmute(self.data), ::std::mem::transmute(other.data))) }
    }

    /// This function returns the complex sine of the complex number z, \sin(z) = (\exp(iz) - \exp(-iz))/(2i).
    pub fn sin(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sin(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex cosine of the complex number z, \cos(z) = (\exp(iz) + \exp(-iz))/2.
    pub fn cos(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_cos(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex tangent of the complex number z, \tan(z) = \sin(z)/\cos(z).
    pub fn tan(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_tan(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex secant of the complex number z, \sec(z) = 1/\cos(z).
    pub fn sec(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sec(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex cosecant of the complex number z, \csc(z) = 1/\sin(z).
    pub fn csc(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_csc(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex cotangent of the complex number z, \cot(z) = 1/\tan(z).
    pub fn cot(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_cot(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsine of the complex number z, \arcsin(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    pub fn arcsin(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arcsin(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsine of the real number z, \arcsin(z).
    /// 
    /// * For z between -1 and 1, the function returns a real value in the range [-\pi/2,\pi/2].
    /// * For z less than -1 the result has a real part of -\pi/2 and a positive imaginary part.
    /// * For z greater than 1 the result has a real part of \pi/2 and a negative imaginary part.
    pub fn arcsin_real(z: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arcsin_real(z)) }
    }

    /// This function returns the complex arccosine of the complex number z, \arccos(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    pub fn arccos(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccos(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arccosine of the real number z, \arccos(z).
    /// 
    /// * For z between -1 and 1, the function returns a real value in the range [0,\pi].
    /// * For z less than -1 the result has a real part of \pi and a negative imaginary part.
    /// * For z greater than 1 the result is purely imaginary and positive.
    pub fn arccos_real(z: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccos_real(z)) }
    }

    /// This function returns the complex arctangent of the complex number z, \arctan(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    pub fn arctan(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arctan(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsecant of the complex number z, \arcsec(z) = \arccos(1/z).
    pub fn arcsec(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arcsec(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arcsecant of the real number z, \arcsec(z) = \arccos(1/z).
    pub fn arcsec_real(z: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arcsec_real(z)) }
    }

    /// This function returns the complex arccosecant of the complex number z, \arccsc(z) = \arcsin(1/z).
    pub fn arccsc(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccsc(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex arccosecant of the real number z, \arccsc(z) = \arcsin(1/z).
    pub fn arccsc_real(z: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccsc_real(z)) }
    }

    /// This function returns the complex arccotangent of the complex number z, \arccot(z) = \arctan(1/z).
    pub fn arccot(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccot(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic sine of the complex number z, \sinh(z) = (\exp(z) - \exp(-z))/2.
    pub fn sinh(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sinh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic cosine of the complex number z, \cosh(z) = (\exp(z) + \exp(-z))/2.
    pub fn cosh(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_cosh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic tangent of the complex number z, \tanh(z) = \sinh(z)/\cosh(z).
    pub fn tanh(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_tanh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic secant of the complex number z, \sech(z) = 1/\cosh(z).
    pub fn sech(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_sech(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic cosecant of the complex number z, \csch(z) = 1/\sinh(z).
    pub fn csch(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_csch(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic cotangent of the complex number z, \coth(z) = 1/\tanh(z).
    pub fn coth(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_coth(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arcsine of the complex number z, \arcsinh(z).
    /// The branch cuts are on the imaginary axis, below -i and above i.
    pub fn arcsinh(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arcsinh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccosine of the complex number z, \arccosh(z).
    /// The branch cut is on the real axis, less than 1.
    /// Note that in this case we use the negative square root in formula 4.6.21 of Abramowitz & Stegun giving \arccosh(z)=\log(z-\sqrt{z^2-1}).
    pub fn arccosh(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccosh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccosine of the real number z, \arccosh(z).
    pub fn arccosh_real(z: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccosh_real(z)) }
    }

    /// This function returns the complex hyperbolic arctangent of the complex number z, \arctanh(z).
    /// The branch cuts are on the real axis, less than -1 and greater than 1.
    pub fn arctanh(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arctanh(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arctangent of the real number z, \arctanh(z).
    pub fn arctanh_real(z: f32) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arctanh_real(z)) }
    }

    /// This function returns the complex hyperbolic arcsecant of the complex number z, \arcsech(z) = \arccosh(1/z).
    pub fn arcsech(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arcsech(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccosecant of the complex number z, \arccsch(z) = \arcsin(1/z).
    pub fn arccsch(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccsch(::std::mem::transmute(self.data))) }
    }

    /// This function returns the complex hyperbolic arccotangent of the complex number z, \arccoth(z) = \arctanh(1/z).
    pub fn arccoth(&self) -> ComplexF32 {
        unsafe { ::std::mem::transmute(::ffi::gsl_complex_float_arccoth(::std::mem::transmute(self.data))) }
    }
}*/

impl Show for ComplexF32 {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "[{}, {}]", self.data[0], self.data[1])
    }
}

impl Clone for ComplexF32 {
    fn clone(&self) -> ComplexF32 {
        ComplexF32 {
            data: [self.data[0], self.data[1]]
        }
    }
}

impl Default for ComplexF32 {
    fn default() -> ComplexF32 {
        ComplexF32 {
            data: [0f32, 0f32]
        }
    }
}