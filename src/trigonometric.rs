//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use enums;
use std::mem::zeroed;
use types;

pub trait Trigonometric {
    /// This routine computes the sine function \sin(x).
    fn sin(&self) -> Self;
    /// This routine computes the sine function \sin(x).
    fn sin_e(&self) -> Result<types::Result, enums::Value>;
    /// This routine computes the cosine function \sin(x).
    fn cos(&self) -> Self;
    /// This routine computes the cosine function \sin(x).
    fn cos_e(&self) -> Result<types::Result, enums::Value>;
    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    fn sf_hypot(&self, y: f64) -> Self;
    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    fn sf_hypot_e(&self, y: f64) -> Result<types::Result, enums::Value>;
    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    fn sinc(&self) -> Self;
    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    fn sinc_e(&self) -> Result<types::Result, enums::Value>;
    /// This function computes the complex sine, \sin(z_r + i z_i) storing the real and imaginary parts in szr, szi.
    fn complex_sin_e(&self, zi: f64) -> Result<(types::Result, types::Result), enums::Value>;
    /// This function computes the complex cosine, \cos(z_r + i z_i) storing the real and imaginary parts in czr, czi.
    fn complex_cos_e(&self, zi: f64) -> Result<(types::Result, types::Result), enums::Value>;
    /// This function computes the logarithm of the complex sine, \log(\sin(z_r + i z_i)) storing the real and imaginary parts in lszr, lszi.
    fn complex_logsin_e(&self, zi: f64) -> Result<(types::Result, types::Result), enums::Value>;
    /// This routine computes \log(\sinh(x)) for x > 0.
    fn lnsinh(&self) -> Self;
    /// This routine computes \log(\sinh(x)) for x > 0.
    fn lnsinh_e(&self) -> Result<types::Result, enums::Value>;
    /// This routine computes \log(\cosh(x)) for x > 0.
    fn lncosh(&self) -> Self;
    /// This routine computes \log(\cosh(x)) for x > 0.
    fn lncosh_e(&self) -> Result<types::Result, enums::Value>;
    /// This function converts the polar coordinates (r,theta) to rectilinear coordinates (x,y), x = r\cos(\theta), y = r\sin(\theta).
    fn polar_to_rect(&self, theta: f64) -> Result<(types::Result, types::Result), enums::Value>;
    /// This function converts the rectilinear coordinates (x,y) to polar coordinates (r,theta), such that x = r\cos(\theta), y = r\sin(\theta).
    /// The argument theta lies in the range [-\pi, \pi].
    fn rect_to_polar(&self, y: f64) -> Result<(types::Result, types::Result), enums::Value>;
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
    fn sin_err_e(&self, dx: f64) -> Result<types::Result, enums::Value>;
    /// This routine computes the cosine of an angle x with an associated absolute error dx, \cos(x \pm dx).
    ///
    /// Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.
    fn cos_err_e(&self, dx: f64) -> Result<types::Result, enums::Value>;
}

impl Trigonometric for f64 {
    fn sin(&self) -> f64 {
        unsafe { ::sys::gsl_sf_sin(*self) }
    }

    fn sin_e(&self) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_sin_e(*self, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }

    fn cos(&self) -> f64 {
        unsafe { ::sys::gsl_sf_cos(*self) }
    }

    fn cos_e(&self) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_cos_e(*self, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }

    fn sf_hypot(&self, y: f64) -> f64 {
        unsafe { ::sys::gsl_sf_hypot(*self, y) }
    }

    fn sf_hypot_e(&self, y: f64) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_hypot_e(*self, y, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }

    fn sinc(&self) -> f64 {
        unsafe { ::sys::gsl_sf_sinc(*self) }
    }

    fn sinc_e(&self) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_sinc_e(*self, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }

    fn complex_sin_e(&self, zi: f64) -> Result<(types::Result, types::Result), enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_complex_sin_e(*self, zi, &mut result, &mut result2) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok((result.into(), result2.into()))
        } else {
            Err(ret)
        }
    }

    fn complex_cos_e(&self, zi: f64) -> Result<(types::Result, types::Result), enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_complex_cos_e(*self, zi, &mut result, &mut result2) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok((result.into(), result2.into()))
        } else {
            Err(ret)
        }
    }

    fn complex_logsin_e(&self, zi: f64) -> Result<(types::Result, types::Result), enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_complex_logsin_e(*self, zi, &mut result, &mut result2) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok((result.into(), result2.into()))
        } else {
            Err(ret)
        }
    }

    fn lnsinh(&self) -> f64 {
        unsafe { ::sys::gsl_sf_lnsinh(*self) }
    }

    fn lnsinh_e(&self) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_lnsinh_e(*self, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }

    fn lncosh(&self) -> f64 {
        unsafe { ::sys::gsl_sf_lncosh(*self) }
    }

    fn lncosh_e(&self) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_lncosh_e(*self, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }

    fn polar_to_rect(&self, theta: f64) -> Result<(types::Result, types::Result), enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_polar_to_rect(*self, theta, &mut result, &mut result2) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok((result.into(), result2.into()))
        } else {
            Err(ret)
        }
    }

    fn rect_to_polar(&self, y: f64) -> Result<(types::Result, types::Result), enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let mut result2 = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_rect_to_polar(*self, y, &mut result, &mut result2) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok((result.into(), result2.into()))
        } else {
            Err(ret)
        }
    }

    fn angle_restrict_symm(&self) -> f64 {
        unsafe { ::sys::gsl_sf_angle_restrict_symm(*self) }
    }

    fn angle_restrict_symm_e(&mut self) -> enums::Value {
        enums::Value::from(unsafe { ::sys::gsl_sf_angle_restrict_symm_e(self) })
    }

    fn angle_restrict_pos(&self) -> f64 {
        unsafe { ::sys::gsl_sf_angle_restrict_pos(*self) }
    }

    fn angle_restrict_pos_e(&mut self) -> enums::Value {
        enums::Value::from(unsafe { ::sys::gsl_sf_angle_restrict_pos_e(self) })
    }

    fn sin_err_e(&self, dx: f64) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_sin_err_e(*self, dx, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }

    fn cos_err_e(&self, dx: f64) -> Result<types::Result, enums::Value> {
        let mut result = unsafe { zeroed::<::sys::gsl_sf_result>() };
        let ret = unsafe { ::sys::gsl_sf_cos_err_e(*self, dx, &mut result) };

        let ret = enums::Value::from(ret);
        if ret.is_success() {
            Ok(result.into())
        } else {
            Err(ret)
        }
    }
}
