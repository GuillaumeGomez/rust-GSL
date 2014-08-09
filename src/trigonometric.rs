//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod Trigonometric {
    use Gsl;
    use std::mem::zeroed;
    use enums;

    /// This routine computes the sine function \sin(x).
    pub fn sin(x: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_sin(x) }
    }

    /// This routine computes the sine function \sin(x).
    pub fn sin_e(x: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_sin_e(x, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This routine computes the cosine function \sin(x).
    pub fn cos(x: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_cos(x) }
    }

    /// This routine computes the cosine function \sin(x).
    pub fn cos_e(x: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_cos_e(x, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    pub fn hypot(x: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_hypot(x) }
    }

    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    pub fn hypot_e(x: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_hypot_e(x, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    pub fn sinc(x: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_sinc(x) }
    }

    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    pub fn sinc_e(x: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_sinc_e(x, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This function computes the complex sine, \sin(z_r + i z_i) storing the real and imaginary parts in szr, szi.
    pub fn complex_sin_e(zr: f64, zi: f64) -> (enums::GslValue, Gsl::Result, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_complex_sin_e(zr, zi, &mut result, &mut result2) };

        (ret, Gsl::Result{val: result.val, err: result.err}, Gsl::Result{val: result2.val, err: result2.err})
    }

    /// This function computes the complex cosine, \cos(z_r + i z_i) storing the real and imaginary parts in czr, czi.
    pub fn complex_cos_e(zr: f64, zi: f64) -> (enums::GslValue, Gsl::Result, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_complex_cos_e(zr, zi, &mut result, &mut result2) };

        (ret, Gsl::Result{val: result.val, err: result.err}, Gsl::Result{val: result2.val, err: result2.err})
    }

    /// This function computes the logarithm of the complex sine, \log(\sin(z_r + i z_i)) storing the real and imaginary parts in lszr, lszi.
    pub fn complex_logsin_e(zr: f64, zi: f64) -> (enums::GslValue, Gsl::Result, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_complex_logsin_e(zr, zi, &mut result, &mut result2) };

        (ret, Gsl::Result{val: result.val, err: result.err}, Gsl::Result{val: result2.val, err: result2.err})
    }

    /// This routine computes \log(\sinh(x)) for x > 0.
    pub fn lnsinh(x: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_lnsinh(x) }
    }

    /// This routine computes \log(\sinh(x)) for x > 0.
    pub fn lnsinh_e(x: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_lnsinh_e(x, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This routine computes \log(\cosh(x)) for x > 0.
    pub fn lncosh(x: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_lncosh(x) }
    }

    /// This routine computes \log(\cosh(x)) for x > 0.
    pub fn lncosh_e(x: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_lncosh_e(x, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This function converts the polar coordinates (r,theta) to rectilinear coordinates (x,y), x = r\cos(\theta), y = r\sin(\theta).
    pub fn polar_to_rect(r: f64, theta: f64) -> (enums::GslValue, Gsl::Result, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_polar_to_rect(r, theta, &mut result, &mut result2) };

        (ret, Gsl::Result{val: result.val, err: result.err}, Gsl::Result{val: result2.val, err: result2.err})
    }

    /// This function converts the rectilinear coordinates (x,y) to polar coordinates (r,theta), such that x = r\cos(\theta), y = r\sin(\theta).
    /// The argument theta lies in the range [-\pi, \pi].
    pub fn rect_to_polar(x: f64, y: f64) -> (enums::GslValue, Gsl::Result, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_rect_to_polar(x, y, &mut result, &mut result2) };

        (ret, Gsl::Result{val: result.val, err: result.err}, Gsl::Result{val: result2.val, err: result2.err})
    }

    /// This routine forces the angle theta to lie in the range (-\pi,\pi].
    /// 
    /// Note that the mathematical value of \pi is slightly greater than M_PI, so the machine numbers M_PI and -M_PI are included in the range.
    pub fn angle_restrict_symm(theta: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_angle_restrict_symm(theta) }
    }

    /// This routine forces the angle theta to lie in the range (-\pi,\pi].
    /// 
    /// Note that the mathematical value of \pi is slightly greater than M_PI, so the machine numbers M_PI and -M_PI are included in the range.
    pub fn angle_restrict_symm_e(theta: &mut f64) -> i32 {
        unsafe { ::ffi::gsl_sf_angle_restrict_symm_e(theta) }
    }

    /// This routine forces the angle theta to lie in the range [0, 2\pi).
    /// 
    /// Note that the mathematical value of 2\pi is slightly greater than 2*M_PI, so the machine number 2*M_PI is included in the range.
    pub fn angle_restrict_pos(theta: f64) -> f64 {
        unsafe { ::ffi::gsl_sf_angle_restrict_pos(theta) }
    }

    /// This routine forces the angle theta to lie in the range [0, 2\pi).
    /// 
    /// Note that the mathematical value of 2\pi is slightly greater than 2*M_PI, so the machine number 2*M_PI is included in the range.
    pub fn angle_restrict_pos_e(theta: &mut f64) -> i32 {
        unsafe { ::ffi::gsl_sf_angle_restrict_pos_e(theta) }
    }

    /// This routine computes the sine of an angle x with an associated absolute error dx, \sin(x \pm dx).
    /// 
    /// Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.
    pub fn sin_err_e(x: f64, dx: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_sin_err_e(x, dx, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This routine computes the cosine of an angle x with an associated absolute error dx, \cos(x \pm dx).
    /// 
    /// Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.
    pub fn cos_err_e(x: f64, dx: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_cos_err_e(x, dx, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }
}