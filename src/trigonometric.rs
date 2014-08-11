//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::mem::zeroed;
use enums;

pub trait Trigonometric {
    /// This routine computes the sine function \sin(x).
    fn sin(&self) -> Self;
    /// This routine computes the sine function \sin(x).
    fn sin_e(&self) -> (enums::Value, ::types::Result);
    /// This routine computes the cosine function \sin(x).
    fn cos(&self) -> Self;
    /// This routine computes the cosine function \sin(x).
    fn cos_e(&self) -> (enums::Value, ::types::Result);
    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    fn sf_hypot(&self) -> Self;
    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    fn sf_hypot_e(&self) -> (enums::Value, ::types::Result);
    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    fn sinc(&self) -> Self;
    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    fn sinc_e(&self) -> (enums::Value, ::types::Result);
    /// This function computes the complex sine, \sin(z_r + i z_i) storing the real and imaginary parts in szr, szi.
    fn complex_sin_e(&self, zi: f64) -> (enums::Value, ::types::Result, ::types::Result);
    /// This function computes the complex cosine, \cos(z_r + i z_i) storing the real and imaginary parts in czr, czi.
    fn complex_cos_e(&self, zi: f64) -> (enums::Value, ::types::Result, ::types::Result);
    /// This function computes the logarithm of the complex sine, \log(\sin(z_r + i z_i)) storing the real and imaginary parts in lszr, lszi.
    fn complex_logsin_e(&self, zi: f64) -> (enums::Value, ::types::Result, ::types::Result);
    /// This routine computes \log(\sinh(x)) for x > 0.
    fn lnsinh(&self) -> Self;
    /// This routine computes \log(\sinh(x)) for x > 0.
    fn lnsinh_e(&self) -> (enums::Value, ::types::Result);
    /// This routine computes \log(\cosh(x)) for x > 0.
    fn lncosh(&self) -> Self;
    /// This routine computes \log(\cosh(x)) for x > 0.
    fn lncosh_e(&self) -> (enums::Value, ::types::Result);
    /// This function converts the polar coordinates (r,theta) to rectilinear coordinates (x,y), x = r\cos(\theta), y = r\sin(\theta).
    fn polar_to_rect(&self, theta: f64) -> (enums::Value, ::types::Result, ::types::Result);
    /// This function converts the rectilinear coordinates (x,y) to polar coordinates (r,theta), such that x = r\cos(\theta), y = r\sin(\theta).
    /// The argument theta lies in the range [-\pi, \pi].
    fn rect_to_polar(&self, y: f64) -> (enums::Value, ::types::Result, ::types::Result);
    /// This routine forces the angle theta to lie in the range (-\pi,\pi].
    /// 
    /// Note that the mathematical value of \pi is slightly greater than M_PI, so the machine numbers M_PI and -M_PI are included in the range.
    fn angle_restrict_symm(&self) -> Self;
    /// This routine forces the angle theta to lie in the range (-\pi,\pi].
    /// 
    /// Note that the mathematical value of \pi is slightly greater than M_PI, so the machine numbers M_PI and -M_PI are included in the range.
    fn angle_restrict_symm_e(&mut self) -> enums::Value;
    /// This routine forces the angle theta to lie in the range [0, 2\pi).
    /// 
    /// Note that the mathematical value of 2\pi is slightly greater than 2*M_PI, so the machine number 2*M_PI is included in the range.
    fn angle_restrict_pos(&self) -> Self;
    /// This routine forces the angle theta to lie in the range [0, 2\pi).
    /// 
    /// Note that the mathematical value of 2\pi is slightly greater than 2*M_PI, so the machine number 2*M_PI is included in the range.
    fn angle_restrict_pos_e(&mut self) -> enums::Value;
    /// This routine computes the sine of an angle x with an associated absolute error dx, \sin(x \pm dx).
    /// 
    /// Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.
    fn sin_err_e(&self, dx: f64) -> (enums::Value, ::types::Result);
    /// This routine computes the cosine of an angle x with an associated absolute error dx, \cos(x \pm dx).
    /// 
    /// Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.
    fn cos_err_e(&self, dx: f64) -> (enums::Value, ::types::Result);
}

impl Trigonometric for f64 {
    fn sin(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_sin(*self) }
    }

    fn sin_e(&self) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_sin_e(*self, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    fn cos(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_cos(*self) }
    }

    fn cos_e(&self) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_cos_e(*self, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    fn sf_hypot(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_hypot(*self) }
    }

    fn sf_hypot_e(&self) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_hypot_e(*self, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    fn sinc(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_sinc(*self) }
    }

    fn sinc_e(&self) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_sinc_e(*self, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    fn complex_sin_e(&self, zi: f64) -> (enums::Value, ::types::Result, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_complex_sin_e(*self, zi, &mut result, &mut result2) };

        (ret, ::types::Result{val: result.val, err: result.err}, ::types::Result{val: result2.val, err: result2.err})
    }

    fn complex_cos_e(&self, zi: f64) -> (enums::Value, ::types::Result, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_complex_cos_e(*self, zi, &mut result, &mut result2) };

        (ret, ::types::Result{val: result.val, err: result.err}, ::types::Result{val: result2.val, err: result2.err})
    }

    fn complex_logsin_e(&self, zi: f64) -> (enums::Value, ::types::Result, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_complex_logsin_e(*self, zi, &mut result, &mut result2) };

        (ret, ::types::Result{val: result.val, err: result.err}, ::types::Result{val: result2.val, err: result2.err})
    }

    fn lnsinh(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_lnsinh(*self) }
    }

    fn lnsinh_e(&self) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_lnsinh_e(*self, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    fn lncosh(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_lncosh(*self) }
    }

    fn lncosh_e(&self) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_lncosh_e(*self, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    fn polar_to_rect(&self, theta: f64) -> (enums::Value, ::types::Result, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_polar_to_rect(*self, theta, &mut result, &mut result2) };

        (ret, ::types::Result{val: result.val, err: result.err}, ::types::Result{val: result2.val, err: result2.err})
    }

    fn rect_to_polar(&self, y: f64) -> (enums::Value, ::types::Result, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_rect_to_polar(*self, y, &mut result, &mut result2) };

        (ret, ::types::Result{val: result.val, err: result.err}, ::types::Result{val: result2.val, err: result2.err})
    }

    fn angle_restrict_symm(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_angle_restrict_symm(*self) }
    }

    fn angle_restrict_symm_e(&mut self) -> enums::Value {
        unsafe { ::ffi::gsl_sf_angle_restrict_symm_e(self) }
    }

    fn angle_restrict_pos(&self) -> f64 {
        unsafe { ::ffi::gsl_sf_angle_restrict_pos(*self) }
    }

    fn angle_restrict_pos_e(&mut self) -> enums::Value {
        unsafe { ::ffi::gsl_sf_angle_restrict_pos_e(self) }
    }

    fn sin_err_e(&self, dx: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_sin_err_e(*self, dx, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    fn cos_err_e(&self, dx: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_cos_err_e(*self, dx, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }
}