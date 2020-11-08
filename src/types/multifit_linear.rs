//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{MatrixF64, Value, VectorF64};
use ffi::FFI;

ffi_wrapper!(
    MultifitLinearWorkspace,
    *mut sys::gsl_multifit_linear_workspace,
    gsl_multifit_linear_free
);

impl MultifitLinearWorkspace {
    pub fn new(n: usize, p: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_multifit_linear_alloc(n, p) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    /// Returns `(Value, chisq)`.
    pub fn linear(
        &mut self,
        x: &MatrixF64,
        y: &VectorF64,
        c: &mut VectorF64,
        cov: &mut MatrixF64,
    ) -> (Value, f64) {
        let mut chisq = 0.;
        let ret = unsafe {
            sys::gsl_multifit_linear(
                x.unwrap_shared(),
                y.unwrap_shared(),
                c.unwrap_unique(),
                cov.unwrap_unique(),
                &mut chisq,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), chisq)
    }

    /// Returns `(Value, chisq, rank)`.
    #[cfg(feature = "v2_3")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_3")))]
    pub fn linear_tsvd(
        &mut self,
        x: &MatrixF64,
        y: &VectorF64,
        tol: f64,
        c: &mut VectorF64,
        cov: &mut MatrixF64,
    ) -> (Value, f64, usize) {
        let mut chisq = 0.;
        let mut rank = 0;
        let ret = unsafe {
            sys::gsl_multifit_linear_tsvd(
                x.unwrap_shared(),
                y.unwrap_shared(),
                tol,
                c.unwrap_unique(),
                cov.unwrap_unique(),
                &mut chisq,
                &mut rank,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), chisq, rank)
    }

    pub fn linear_svd(&mut self, x: &MatrixF64) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_svd(
                x.unwrap_shared(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_linear_bsvd(&mut self, x: &MatrixF64) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_bsvd(
                x.unwrap_shared(),
                self.unwrap_unique(),
            ))
        }
    }

    #[cfg(feature = "v2_3")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_3")))]
    pub fn linear_rank(&self, tol: f64) -> usize {
        unsafe { sys::gsl_multifit_linear_rank(tol, self.unwrap_shared()) }
    }

    /// Returns `(Value, rnorm, snorm)`.
    pub fn linear_solve(
        &mut self,
        lambda: f64,
        x: &MatrixF64,
        y: &VectorF64,
        c: &mut VectorF64,
    ) -> (Value, f64, f64) {
        let mut rnorm = 0.;
        let mut snorm = 0.;
        let ret = unsafe {
            sys::gsl_multifit_linear_solve(
                lambda,
                x.unwrap_shared(),
                y.unwrap_shared(),
                c.unwrap_unique(),
                &mut rnorm,
                &mut snorm,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), rnorm, snorm)
    }

    pub fn linear_stdform1(
        &mut self,
        l: &VectorF64,
        x: &MatrixF64,
        y: &VectorF64,
        xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_stdform1(
                l.unwrap_shared(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_wstdform1(
        &mut self,
        l: &VectorF64,
        x: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_wstdform1(
                l.unwrap_shared(),
                x.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_stdform2(
        &mut self,
        lqr: &MatrixF64,
        ltau: &VectorF64,
        x: &MatrixF64,
        y: &VectorF64,
        xs: &mut MatrixF64,
        ys: &mut VectorF64,
        m: &mut MatrixF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_stdform2(
                lqr.unwrap_shared(),
                ltau.unwrap_shared(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                xs.unwrap_unique(),
                ys.unwrap_unique(),
                m.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_wstdform2(
        &mut self,
        lqr: &MatrixF64,
        ltau: &VectorF64,
        x: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        xs: &mut MatrixF64,
        ys: &mut VectorF64,
        m: &mut MatrixF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_wstdform2(
                lqr.unwrap_shared(),
                ltau.unwrap_shared(),
                x.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                xs.unwrap_unique(),
                ys.unwrap_unique(),
                m.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_genform1(&mut self, l: &VectorF64, cs: &VectorF64, c: &mut VectorF64) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_genform1(
                l.unwrap_shared(),
                cs.unwrap_shared(),
                c.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_genform2(
        &mut self,
        lqr: &MatrixF64,
        ltau: &VectorF64,
        x: &MatrixF64,
        y: &VectorF64,
        cs: &VectorF64,
        m: &MatrixF64,
        c: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_genform2(
                lqr.unwrap_shared(),
                ltau.unwrap_shared(),
                x.unwrap_shared(),
                y.unwrap_shared(),
                cs.unwrap_shared(),
                m.unwrap_shared(),
                c.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_wgenform2(
        &mut self,
        lqr: &MatrixF64,
        ltau: &VectorF64,
        x: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        cs: &VectorF64,
        m: &MatrixF64,
        c: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_wgenform2(
                lqr.unwrap_shared(),
                ltau.unwrap_shared(),
                x.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                cs.unwrap_shared(),
                m.unwrap_shared(),
                c.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_lcurve(
        &mut self,
        y: &VectorF64,
        reg_param: &mut VectorF64,
        rho: &mut VectorF64,
        eta: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_lcurve(
                y.unwrap_shared(),
                reg_param.unwrap_unique(),
                rho.unwrap_unique(),
                eta.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn linear_Lsobolev(
        &mut self,
        p: usize,
        kmax: usize,
        alpha: &VectorF64,
        l: &mut MatrixF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_Lsobolev(
                p,
                kmax,
                alpha.unwrap_shared(),
                l.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    /// Returns `(Value, chisq)`.
    pub fn wlinear(
        &mut self,
        x: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        c: &mut VectorF64,
        cov: &mut MatrixF64,
    ) -> (Value, f64) {
        let mut chisq = 0.;
        let ret = unsafe {
            sys::gsl_multifit_wlinear(
                x.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                c.unwrap_unique(),
                cov.unwrap_unique(),
                &mut chisq,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), chisq)
    }

    /// Returns `(Value, chisq, rank)`.
    #[cfg(feature = "v2_3")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_3")))]
    pub fn wlinear_tsvd(
        &mut self,
        x: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        tol: f64,
        c: &mut VectorF64,
        cov: &mut MatrixF64,
    ) -> (Value, f64, usize) {
        let mut chisq = 0.;
        let mut rank = 0;
        let ret = unsafe {
            sys::gsl_multifit_wlinear_tsvd(
                x.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                tol,
                c.unwrap_unique(),
                cov.unwrap_unique(),
                &mut chisq,
                &mut rank,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), chisq, rank)
    }

    /// Returns `(Value, rank, chisq)`.
    pub fn wlinear_svd(
        &mut self,
        x: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        tol: f64,
        c: &mut VectorF64,
        cov: &mut MatrixF64,
    ) -> (Value, usize, f64) {
        let mut rank = 0;
        let mut chisq = 0.;
        let ret = unsafe {
            sys::gsl_multifit_wlinear_svd(
                x.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                tol,
                &mut rank,
                c.unwrap_unique(),
                cov.unwrap_unique(),
                &mut chisq,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), rank, chisq)
    }

    /// Returns `(Value, rank, chisq)`.
    pub fn wlinear_usvd(
        &mut self,
        x: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        tol: f64,
        c: &mut VectorF64,
        cov: &mut MatrixF64,
    ) -> (Value, usize, f64) {
        let mut rank = 0;
        let mut chisq = 0.;
        let ret = unsafe {
            sys::gsl_multifit_wlinear_usvd(
                x.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                tol,
                &mut rank,
                c.unwrap_unique(),
                cov.unwrap_unique(),
                &mut chisq,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), rank, chisq)
    }

    #[cfg(feature = "v2_1")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_1")))]
    pub fn linear_rcond(&mut self) -> f64 {
        unsafe { sys::gsl_multifit_linear_rcond(self.unwrap_unique()) }
    }

    /// Returns `(Value, delta0)`.
    pub fn linear_gcv_init(
        &mut self,
        y: &VectorF64,
        reg_param: &mut VectorF64,
        UTy: &mut VectorF64,
    ) -> (Value, f64) {
        let mut delta0 = 0.;
        let ret = unsafe {
            sys::gsl_multifit_linear_gcv_init(
                y.unwrap_shared(),
                reg_param.unwrap_unique(),
                UTy.unwrap_unique(),
                &mut delta0,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), delta0)
    }

    pub fn linear_gcv_curve(
        &mut self,
        reg_param: &VectorF64,
        UTy: &VectorF64,
        delta0: f64,
        g: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multifit_linear_gcv_curve(
                reg_param.unwrap_shared(),
                UTy.unwrap_shared(),
                delta0,
                g.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    /// Returns `(Value, lambda)`.
    pub fn linear_gcv_min(
        &mut self,
        reg_param: &VectorF64,
        UTy: &VectorF64,
        g: &VectorF64,
        delta0: f64,
    ) -> (Value, f64) {
        let mut lambda = 0.;
        let ret = unsafe {
            sys::gsl_multifit_linear_gcv_min(
                reg_param.unwrap_shared(),
                UTy.unwrap_shared(),
                g.unwrap_shared(),
                delta0,
                &mut lambda,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), lambda)
    }

    pub fn linear_gcv_calc(&mut self, lambda: f64, UTy: &VectorF64, delta0: f64) -> f64 {
        unsafe {
            sys::gsl_multifit_linear_gcv_calc(
                lambda,
                UTy.unwrap_shared(),
                delta0,
                self.unwrap_unique(),
            )
        }
    }

    /// Returns `(Value, lambda, g_lambda)`.
    pub fn linear_gcv(
        &mut self,
        y: &VectorF64,
        reg_param: &mut VectorF64,
        g: &mut VectorF64,
    ) -> (Value, f64, f64) {
        let mut lambda = 0.;
        let mut g_lambda = 0.;
        let ret = unsafe {
            sys::gsl_multifit_linear_gcv(
                y.unwrap_shared(),
                reg_param.unwrap_unique(),
                g.unwrap_unique(),
                &mut lambda,
                &mut g_lambda,
                self.unwrap_unique(),
            )
        };
        (Value::from(ret), lambda, g_lambda)
    }
}
