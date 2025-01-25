//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::types;
use crate::Error;
use std::mem::MaybeUninit;

pub trait Trigonometric {
    /// This routine computes the sine function \sin(x).
    fn sin(&self) -> Self;
    /// This routine computes the sine function \sin(x).
    fn sin_e(&self) -> Result<types::Result, Error>;
    /// This routine computes the cosine function \sin(x).
    fn cos(&self) -> Self;
    /// This routine computes the cosine function \sin(x).
    fn cos_e(&self) -> Result<types::Result, Error>;
    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    fn sf_hypot(&self, y: f64) -> Self;
    /// This routine computes the hypotenuse function \sqrt{x^2 + y^2} avoiding overflow and underflow.
    fn sf_hypot_e(&self, y: f64) -> Result<types::Result, Error>;
    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    fn sinc(&self) -> Self;
    /// This routine computes \sinc(x) = \sin(\pi x) / (\pi x) for any value of x.
    fn sinc_e(&self) -> Result<types::Result, Error>;
    /// This function computes the complex sine, \sin(z_r + i z_i) storing the real and imaginary parts in szr, szi.
    fn complex_sin_e(&self, zi: f64) -> Result<(types::Result, types::Result), Error>;
    /// This function computes the complex cosine, \cos(z_r + i z_i) storing the real and imaginary parts in czr, czi.
    fn complex_cos_e(&self, zi: f64) -> Result<(types::Result, types::Result), Error>;
    /// This function computes the logarithm of the complex sine, \log(\sin(z_r + i z_i)) storing the real and imaginary parts in lszr, lszi.
    fn complex_logsin_e(&self, zi: f64) -> Result<(types::Result, types::Result), Error>;
    /// This routine computes \log(\sinh(x)) for x > 0.
    fn lnsinh(&self) -> Self;
    /// This routine computes \log(\sinh(x)) for x > 0.
    fn lnsinh_e(&self) -> Result<types::Result, Error>;
    /// This routine computes \log(\cosh(x)) for x > 0.
    fn lncosh(&self) -> Self;
    /// This routine computes \log(\cosh(x)) for x > 0.
    fn lncosh_e(&self) -> Result<types::Result, Error>;
    /// This function converts the polar coordinates (r,theta) to rectilinear coordinates (x,y), x = r\cos(\theta), y = r\sin(\theta).
    fn polar_to_rect(&self, theta: f64) -> Result<(types::Result, types::Result), Error>;
    /// This function converts the rectilinear coordinates (x,y) to polar coordinates (r,theta), such that x = r\cos(\theta), y = r\sin(\theta).
    /// The argument theta lies in the range [-\pi, \pi].
    fn rect_to_polar(&self, y: f64) -> Result<(types::Result, types::Result), Error>;
    /// This routine forces the angle theta to lie in the range (-\pi,\pi].
    ///
    /// Note that the mathematical value of \pi is slightly greater than M_PI, so the machine numbers M_PI and -M_PI are included in the range.
    fn angle_restrict_symm(&self) -> Self;
    /// This routine forces the angle theta to lie in the range (-\pi,\pi].
    ///
    /// Note that the mathematical value of \pi is slightly greater than M_PI, so the machine numbers M_PI and -M_PI are included in the range.
    fn angle_restrict_symm_e(&mut self) -> Result<(), Error>;
    /// This routine forces the angle theta to lie in the range [0, 2\pi).
    ///
    /// Note that the mathematical value of 2\pi is slightly greater than 2*M_PI, so the machine number 2*M_PI is included in the range.
    fn angle_restrict_pos(&self) -> Self;
    /// This routine forces the angle theta to lie in the range [0, 2\pi).
    ///
    /// Note that the mathematical value of 2\pi is slightly greater than 2*M_PI, so the machine number 2*M_PI is included in the range.
    fn angle_restrict_pos_e(&mut self) -> Result<(), Error>;
    /// This routine computes the sine of an angle x with an associated absolute error dx, \sin(x \pm dx).
    ///
    /// Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.
    fn sin_err_e(&self, dx: f64) -> Result<types::Result, Error>;
    /// This routine computes the cosine of an angle x with an associated absolute error dx, \cos(x \pm dx).
    ///
    /// Note that this function is provided in the error-handling form only since its purpose is to compute the propagated error.
    fn cos_err_e(&self, dx: f64) -> Result<types::Result, Error>;
}

impl Trigonometric for f64 {
    #[doc(alias = "gsl_sf_sin")]
    fn sin(&self) -> f64 {
        unsafe { sys::gsl_sf_sin(*self) }
    }

    #[doc(alias = "gsl_sf_sin_e")]
    fn sin_e(&self) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_sin_e(*self, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[doc(alias = "gsl_sf_cos")]
    fn cos(&self) -> f64 {
        unsafe { sys::gsl_sf_cos(*self) }
    }

    #[doc(alias = "gsl_sf_cos_e")]
    fn cos_e(&self) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_cos_e(*self, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[doc(alias = "gsl_sf_hypot")]
    fn sf_hypot(&self, y: f64) -> f64 {
        unsafe { sys::gsl_sf_hypot(*self, y) }
    }

    #[doc(alias = "gsl_sf_hypot_e")]
    fn sf_hypot_e(&self, y: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_hypot_e(*self, y, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[doc(alias = "gsl_sf_sinc")]
    fn sinc(&self) -> f64 {
        unsafe { sys::gsl_sf_sinc(*self) }
    }

    #[doc(alias = "gsl_sf_sinc_e")]
    fn sinc_e(&self) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_sinc_e(*self, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[doc(alias = "gsl_sf_complex_sin_e")]
    fn complex_sin_e(&self, zi: f64) -> Result<(types::Result, types::Result), Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let mut result2 = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe {
            sys::gsl_sf_complex_sin_e(*self, zi, result.as_mut_ptr(), result2.as_mut_ptr())
        };

        Error::handle(
            ret,
            (
                unsafe { result.assume_init() }.into(),
                unsafe { result2.assume_init() }.into(),
            ),
        )
    }

    #[doc(alias = "gsl_sf_complex_cos_e")]
    fn complex_cos_e(&self, zi: f64) -> Result<(types::Result, types::Result), Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let mut result2 = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe {
            sys::gsl_sf_complex_cos_e(*self, zi, result.as_mut_ptr(), result2.as_mut_ptr())
        };

        Error::handle(
            ret,
            (
                unsafe { result.assume_init() }.into(),
                unsafe { result2.assume_init() }.into(),
            ),
        )
    }

    #[doc(alias = "gsl_sf_complex_logsin_e")]
    fn complex_logsin_e(&self, zi: f64) -> Result<(types::Result, types::Result), Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let mut result2 = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe {
            sys::gsl_sf_complex_logsin_e(*self, zi, result.as_mut_ptr(), result2.as_mut_ptr())
        };

        Error::handle(
            ret,
            (
                unsafe { result.assume_init() }.into(),
                unsafe { result2.assume_init() }.into(),
            ),
        )
    }

    #[doc(alias = "gsl_sf_lnsinh")]
    fn lnsinh(&self) -> f64 {
        unsafe { sys::gsl_sf_lnsinh(*self) }
    }

    #[doc(alias = "gsl_sf_lnsinh_e")]
    fn lnsinh_e(&self) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lnsinh_e(*self, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[doc(alias = "gsl_sf_lncosh")]
    fn lncosh(&self) -> f64 {
        unsafe { sys::gsl_sf_lncosh(*self) }
    }

    #[doc(alias = "gsl_sf_lncosh_e")]
    fn lncosh_e(&self) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lncosh_e(*self, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[doc(alias = "gsl_sf_polar_to_rect")]
    fn polar_to_rect(&self, theta: f64) -> Result<(types::Result, types::Result), Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let mut result2 = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe {
            sys::gsl_sf_polar_to_rect(*self, theta, result.as_mut_ptr(), result2.as_mut_ptr())
        };

        Error::handle(
            ret,
            (
                unsafe { result.assume_init() }.into(),
                unsafe { result2.assume_init() }.into(),
            ),
        )
    }

    #[doc(alias = "gsl_sf_rect_to_polar")]
    fn rect_to_polar(&self, y: f64) -> Result<(types::Result, types::Result), Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let mut result2 = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe {
            sys::gsl_sf_rect_to_polar(*self, y, result.as_mut_ptr(), result2.as_mut_ptr())
        };

        Error::handle(
            ret,
            (
                unsafe { result.assume_init() }.into(),
                unsafe { result2.assume_init() }.into(),
            ),
        )
    }

    #[doc(alias = "gsl_sf_angle_restrict_symm")]
    fn angle_restrict_symm(&self) -> f64 {
        unsafe { sys::gsl_sf_angle_restrict_symm(*self) }
    }

    #[doc(alias = "gsl_sf_angle_restrict_symm_e")]
    fn angle_restrict_symm_e(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_sf_angle_restrict_symm_e(self) };
        Error::handle(ret, ())
    }

    #[doc(alias = "gsl_sf_angle_restrict_pos")]
    fn angle_restrict_pos(&self) -> f64 {
        unsafe { sys::gsl_sf_angle_restrict_pos(*self) }
    }

    #[doc(alias = "gsl_sf_angle_restrict_pos_e")]
    fn angle_restrict_pos_e(&mut self) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_sf_angle_restrict_pos_e(self) };
        Error::handle(ret, ())
    }

    #[doc(alias = "gsl_sf_sin_err_e")]
    fn sin_err_e(&self, dx: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_sin_err_e(*self, dx, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    #[doc(alias = "gsl_sf_cos_err_e")]
    fn cos_err_e(&self, dx: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_cos_err_e(*self, dx, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }
}
