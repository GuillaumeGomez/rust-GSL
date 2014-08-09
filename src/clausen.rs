//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The Clausen function is defined by the following integral,

Cl_2(x) = - \int_0^x dt \log(2 \sin(t/2))

It is related to the dilogarithm by Cl_2(\theta) = \Im Li_2(\exp(i\theta)). 
!*/

pub mod Clausen {
    use ffi;
    use Gsl;
    use enums;
    use std::mem::zeroed;

    /// This routine computes the Clausen integral Cl_2(x).
    pub fn clausen(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_clausen(x) }
    }

    /// This routine computes the Clausen integral Cl_2(x).
    pub fn clausen_e(x: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_clausen_e(x, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }
}