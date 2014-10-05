//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The polygamma functions of order n are defined by

\psi^{(n)}(x) = (d/dx)^n \psi(x) = (d/dx)^{n+1} \log(\Gamma(x))

where \psi(x) = \Gamma'(x)/\Gamma(x) is known as the digamma function.
!*/

pub mod diagamma {
    use std::mem::zeroed;
    use enums;
    use ffi;

    /// This routine computes the digamma function \psi(n) for positive integer n. The digamma function is also called the Psi function.
    pub fn psi_int(n: i32) -> f64 {
        unsafe { ffi::gsl_sf_psi_int(n) }
    }

    /// This routine computes the digamma function \psi(n) for positive integer n. The digamma function is also called the Psi function.
    pub fn psi_int_e(n: i32) -> (enums::value::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_psi_int_e(n, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the digamma function \psi(x) for general x, x \ne 0.
    pub fn psi(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_psi(x) }
    }

    /// This routine computes the digamma function \psi(x) for general x, x \ne 0.
    pub fn psi_e(x: f64) -> (enums::value::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_psi_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the real part of the digamma function on the line 1+i y, \Re[\psi(1 + i y)].
    pub fn psi_1piy(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_psi_1piy(x) }
    }

    /// This routine computes the real part of the digamma function on the line 1+i y, \Re[\psi(1 + i y)].
    pub fn psi_1piy_e(x: f64) -> (enums::value::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_psi_1piy_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }
}

pub mod trigamma {
    use std::mem::zeroed;
    use enums;
    use ffi;

    /// This routine computes the Trigamma function \psi'(n) for positive integer n.
    pub fn psi_1_int(n: i32) -> f64 {
        unsafe { ffi::gsl_sf_psi_1_int(n) }
    }

    /// This routine computes the Trigamma function \psi'(n) for positive integer n.
    pub fn psi_1_int_e(n: i32) -> (enums::value::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_psi_1_int_e(n, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the Trigamma function \psi'(x) for general x.
    pub fn psi_1(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_psi_1(x) }
    }

    /// This routine computes the Trigamma function \psi'(x) for general x.
    pub fn psi_1_e(x: f64) -> (enums::value::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_psi_1_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }
}

pub mod polygamma {
    use std::mem::zeroed;
    use enums;
    use ffi;

    /// This routine computes the polygamma function \psi^{(n)}(x) for n >= 0, x > 0.
    pub fn psi_n(n: i32, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_psi_n(n, x) }
    }

    /// This routine computes the polygamma function \psi^{(n)}(x) for n >= 0, x > 0.
    pub fn psi_n_e(n: i32, x: f64) -> (enums::value::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_psi_n_e(n, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }
}