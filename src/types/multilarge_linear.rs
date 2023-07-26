//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{MatrixF64, Value, VectorF64};
use ffi::FFI;

ffi_wrapper!(MultilargeLinearType, *const sys::gsl_multilarge_linear_type);

impl MultilargeLinearType {
    #[doc(alias = "gsl_multilarge_linear_normal")]
    pub fn normal() -> MultilargeLinearType {
        ffi_wrap!(gsl_multilarge_linear_normal)
    }

    #[doc(alias = "gsl_multilarge_linear_tsqr")]
    pub fn tsqr() -> MultilargeLinearType {
        ffi_wrap!(gsl_multilarge_linear_tsqr)
    }
}

ffi_wrapper!(
    MultilargeLinearWorkspace,
    *mut sys::gsl_multilarge_linear_workspace,
    gsl_multilarge_linear_free
);

impl MultilargeLinearWorkspace {
    #[doc(alias = "gsl_multilarge_linear_alloc")]
    pub fn new(t: MultilargeLinearType, p: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_multilarge_linear_alloc(t.unwrap_shared(), p) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

    #[doc(alias = "gsl_multilarge_linear_name")]
    pub fn name(&self) -> Option<String> {
        let n = unsafe { sys::gsl_multilarge_linear_name(self.unwrap_shared()) };
        if n.is_null() {
            return None;
        }
        let mut len = 0;
        loop {
            if unsafe { *n.offset(len) } == 0 {
                break;
            }
            len += 1;
        }
        let slice = unsafe { std::slice::from_raw_parts(n as _, len as _) };
        std::str::from_utf8(slice).ok().map(|x| x.to_owned())
    }

    #[doc(alias = "gsl_multilarge_linear_reset")]
    pub fn reset(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_multilarge_linear_reset(self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_multilarge_linear_accumulate")]
    pub fn accumulate(&mut self, x: &mut MatrixF64, y: &mut VectorF64) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_accumulate(
                x.unwrap_unique(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// Returns `(rnorm, snorm)`.
    #[doc(alias = "gsl_multilarge_linear_solve")]
    pub fn solve(&mut self, lambda: f64, c: &mut VectorF64) -> Result<(f64, f64), Value> {
        let mut rnorm = 0.;
        let mut snorm = 0.;
        let ret = unsafe {
            sys::gsl_multilarge_linear_solve(
                lambda,
                c.unwrap_unique(),
                &mut rnorm,
                &mut snorm,
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, (rnorm, snorm))
    }

    /// Returns `rcond`.
    #[doc(alias = "gsl_multilarge_linear_rcond")]
    pub fn rcond(&mut self) -> Result<f64, Value> {
        let mut rcond = 0.;
        let ret = unsafe { sys::gsl_multilarge_linear_rcond(&mut rcond, self.unwrap_unique()) };
        result_handler!(ret, rcond)
    }

    #[cfg(feature = "v2_2")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
    #[doc(alias = "gsl_multilarge_linear_lcurve")]
    pub fn lcurve(
        &mut self,
        reg_param: &mut VectorF64,
        rho: &mut VectorF64,
        eta: &mut VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_lcurve(
                reg_param.unwrap_unique(),
                rho.unwrap_unique(),
                eta.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_multilarge_linear_wstdform1")]
    pub fn wstdform1(
        &mut self,
        L: &VectorF64,
        X: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_wstdform1(
                L.unwrap_shared(),
                X.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_multilarge_linear_stdform1")]
    pub fn stdform1(
        &mut self,
        L: &VectorF64,
        X: &MatrixF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_stdform1(
                L.unwrap_shared(),
                X.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_multilarge_linear_wstdform2")]
    pub fn wstdform2(
        &mut self,
        LQR: &MatrixF64,
        Ltau: &VectorF64,
        X: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_wstdform2(
                LQR.unwrap_shared(),
                Ltau.unwrap_shared(),
                X.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_multilarge_linear_stdform2")]
    pub fn stdform2(
        &mut self,
        LQR: &MatrixF64,
        Ltau: &VectorF64,
        X: &MatrixF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_stdform2(
                LQR.unwrap_shared(),
                Ltau.unwrap_shared(),
                X.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_multilarge_linear_genform1")]
    pub fn genform1(
        &mut self,
        L: &VectorF64,
        cs: &VectorF64,
        c: &mut VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_genform1(
                L.unwrap_shared(),
                cs.unwrap_shared(),
                c.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = "gsl_multilarge_linear_genform2")]
    pub fn genform2(
        &mut self,
        LQR: &MatrixF64,
        Ltau: &VectorF64,
        cs: &VectorF64,
        c: &mut VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_multilarge_linear_genform2(
                LQR.unwrap_shared(),
                Ltau.unwrap_shared(),
                cs.unwrap_shared(),
                c.unwrap_unique(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[cfg(feature = "v2_7")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_7")))]
    #[doc(alias = "gsl_multilarge_linear_matrix_ptr")]
    pub fn matrix<F: FnOnce(&MatrixF64)>(&self, f: F) {
        f(&MatrixF64::soft_wrap(unsafe {
            sys::gsl_multilarge_linear_matrix_ptr(self.unwrap_shared()) as _
        }))
    }

    #[cfg(feature = "v2_7")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_7")))]
    #[doc(alias = "gsl_multilarge_linear_rhs_ptr")]
    pub fn rhs<F: FnOnce(&VectorF64)>(&self, f: F) {
        f(&VectorF64::soft_wrap(unsafe {
            sys::gsl_multilarge_linear_rhs_ptr(self.unwrap_shared()) as _
        }))
    }
}
