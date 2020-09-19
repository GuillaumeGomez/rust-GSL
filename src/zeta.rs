//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The Riemann zeta function is defined in Abramowitz & Stegun, Section 23.2.

/// The Riemann zeta function is defined by the infinite sum \zeta(s) = \sum_{k=1}^\infty k^{-s}.
pub mod riemann {
    use enums;
    use std::mem::zeroed;

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    pub fn zeta_int(n: i32) -> f64 {
        unsafe { sys::gsl_sf_zeta_int(n) }
    }

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    pub fn zeta_int_e(n: i32) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
        let ret = unsafe { sys::gsl_sf_zeta_int_e(n, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }

    /// This routine computes the Riemann zeta function \zeta(s) for arbitrary s, s \ne 1.
    pub fn zeta(x: f64) -> f64 {
        unsafe { sys::gsl_sf_zeta(x) }
    }

    /// This routine computes the Riemann zeta function \zeta(s) for arbitrary s, s \ne 1.
    pub fn zeta_e(x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
        let ret = unsafe { sys::gsl_sf_zeta_e(x, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }
}

/// For large positive argument, the Riemann zeta function approaches one.
/// In this region the fractional part is interesting, and therefore we need a function to evaluate it explicitly.
pub mod riemann_mins_one {
    use enums;
    use std::mem::zeroed;

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    pub fn zetam1_int(n: i32) -> f64 {
        unsafe { sys::gsl_sf_zetam1_int(n) }
    }

    /// This routine computes the Riemann zeta function \zeta(n) for integer n, n \ne 1.
    pub fn zetam1_int_e(n: i32) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
        let ret = unsafe { sys::gsl_sf_zetam1_int_e(n, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }

    /// This routine computes \zeta(s) - 1 for arbitrary s, s \ne 1.
    pub fn zetam1(x: f64) -> f64 {
        unsafe { sys::gsl_sf_zetam1(x) }
    }

    /// This routine computes \zeta(s) - 1 for arbitrary s, s \ne 1.
    pub fn zetam1_e(x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
        let ret = unsafe { sys::gsl_sf_zetam1_e(x, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }
}

/// The Hurwitz zeta function is defined by \zeta(s,q) = \sum_0^\infty (k+q)^{-s}.
pub mod hurwitz {
    use enums;
    use std::mem::zeroed;

    /// This routine computes the Hurwitz zeta function \zeta(s,q) for s > 1, q > 0.
    pub fn hzeta(s: f64, q: f64) -> f64 {
        unsafe { sys::gsl_sf_hzeta(s, q) }
    }

    /// This routine computes the Hurwitz zeta function \zeta(s,q) for s > 1, q > 0.
    pub fn hzeta_e(s: f64, q: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
        let ret = unsafe { sys::gsl_sf_hzeta_e(s, q, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }
}

/// The eta function is defined by \eta(s) = (1-2^{1-s}) \zeta(s).
pub mod eta {
    use enums;
    use std::mem::zeroed;

    /// This routine computes the eta function \eta(n) for integer n.
    pub fn eta_int(n: i32) -> f64 {
        unsafe { sys::gsl_sf_eta_int(n) }
    }

    /// This routine computes the eta function \eta(n) for integer n.
    pub fn eta_int_e(n: i32) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
        let ret = unsafe { sys::gsl_sf_eta_int_e(n, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }

    /// This routine computes the eta function \eta(s) for arbitrary s.
    pub fn eta(s: f64) -> f64 {
        unsafe { sys::gsl_sf_eta(s) }
    }

    /// This routine computes the eta function \eta(s) for arbitrary s.
    pub fn eta_e(s: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<sys::gsl_sf_result>() };
        let ret = unsafe { sys::gsl_sf_eta_e(s, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }
}
