//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The Riemann zeta function is defined in Abramowitz & Stegun, Section 23.2.

/// The Riemann zeta function is defined by the infinite sum \zeta(s) = \sum_{k=1}^\infty k^{-s}.
pub mod riemann {
    use crate::Value;
    use std::mem::MaybeUninit;

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    #[doc(alias = "gsl_sf_zeta_int")]
    pub fn zeta_int(n: i32) -> f64 {
        unsafe { sys::gsl_sf_zeta_int(n) }
    }

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    #[doc(alias = "gsl_sf_zeta_int_e")]
    pub fn zeta_int_e(n: i32) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_zeta_int_e(n, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the Riemann zeta function \zeta(s) for arbitrary s, s \ne 1.
    #[doc(alias = "gsl_sf_zeta")]
    pub fn zeta(x: f64) -> f64 {
        unsafe { sys::gsl_sf_zeta(x) }
    }

    /// This routine computes the Riemann zeta function \zeta(s) for arbitrary s, s \ne 1.
    #[doc(alias = "gsl_sf_zeta_e")]
    pub fn zeta_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_zeta_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }
}

/// For large positive argument, the Riemann zeta function approaches one.
/// In this region the fractional part is interesting, and therefore we need a function to evaluate it explicitly.
pub mod riemann_mins_one {
    use crate::Value;
    use std::mem::MaybeUninit;

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    #[doc(alias = "gsl_sf_zetam1_int")]
    pub fn zetam1_int(n: i32) -> f64 {
        unsafe { sys::gsl_sf_zetam1_int(n) }
    }

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    #[doc(alias = "gsl_sf_zetam1_int_e")]
    pub fn zetam1_int_e(n: i32) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_zetam1_int_e(n, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes \zeta(s) - 1 for arbitrary s, s \ne 1.
    #[doc(alias = "gsl_sf_zetam1")]
    pub fn zetam1(x: f64) -> f64 {
        unsafe { sys::gsl_sf_zetam1(x) }
    }

    /// This routine computes \zeta(s) - 1 for arbitrary s, s \ne 1.
    #[doc(alias = "gsl_sf_zetam1_e")]
    pub fn zetam1_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_zetam1_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }
}

/// The Hurwitz zeta function is defined by \zeta(s,q) = \sum_0^\infty (k+q)^{-s}.
pub mod hurwitz {
    use crate::Value;
    use std::mem::MaybeUninit;

    /// This routine computes the Hurwitz zeta function \zeta(s,q) for s > 1, q > 0.
    #[doc(alias = "gsl_sf_hzeta")]
    pub fn hzeta(s: f64, q: f64) -> f64 {
        unsafe { sys::gsl_sf_hzeta(s, q) }
    }

    /// This routine computes the Hurwitz zeta function \zeta(s,q) for s > 1, q > 0.
    #[doc(alias = "gsl_sf_hzeta_e")]
    pub fn hzeta_e(s: f64, q: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_hzeta_e(s, q, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }
}

/// The eta function is defined by \eta(s) = (1-2^{1-s}) \zeta(s).
pub mod eta {
    use crate::Value;
    use std::mem::MaybeUninit;

    /// This routine computes the eta function \eta(n) for integer n.
    #[doc(alias = "gsl_sf_eta_int")]
    pub fn eta_int(n: i32) -> f64 {
        unsafe { sys::gsl_sf_eta_int(n) }
    }

    /// This routine computes the eta function \eta(n) for integer n.
    #[doc(alias = "gsl_sf_eta_int_e")]
    pub fn eta_int_e(n: i32) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_eta_int_e(n, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the eta function \eta(s) for arbitrary s.
    #[doc(alias = "gsl_sf_eta")]
    pub fn eta(s: f64) -> f64 {
        unsafe { sys::gsl_sf_eta(s) }
    }

    /// This routine computes the eta function \eta(s) for arbitrary s.
    #[doc(alias = "gsl_sf_eta_e")]
    pub fn eta_e(s: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_eta_e(s, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }
}
