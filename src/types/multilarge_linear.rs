//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::{MatrixF64, Value, VectorF64};
use ffi::FFI;

ffi_wrapper!(MultilargeLinearType, *const sys::gsl_multilarge_linear_type);

impl MultilargeLinearType {
    pub fn normal() -> MultilargeLinearType {
        ffi_wrap!(gsl_multilarge_linear_normal)
    }

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
    pub fn new(t: MultilargeLinearType, p: usize) -> Option<Self> {
        let s = unsafe { sys::gsl_multilarge_linear_alloc(t.unwrap_shared(), p) };
        if s.is_null() {
            None
        } else {
            Some(Self::wrap(s))
        }
    }

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
        let slice = unsafe { ::std::slice::from_raw_parts(n as _, len as _) };
        ::std::str::from_utf8(slice).ok().map(|x| x.to_owned())
    }

    pub fn reset(&mut self) -> Value {
        unsafe { Value::from(sys::gsl_multilarge_linear_reset(self.unwrap_unique())) }
    }

    pub fn accumulate(&mut self, x: &mut MatrixF64, y: &mut VectorF64) -> Value {
        unsafe {
            Value::from(sys::gsl_multilarge_linear_accumulate(
                x.unwrap_unique(),
                y.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    /// Returns `(Value, rnorm, snorm)`.
    pub fn solve(&mut self, lambda: f64, c: &mut VectorF64) -> (Value, f64, f64) {
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
        (Value::from(ret), rnorm, snorm)
    }

    /// Returns `(Value, rcond)`.
    pub fn rcond(&mut self) -> (Value, f64) {
        let mut rcond = 0.;
        let ret = unsafe { sys::gsl_multilarge_linear_rcond(&mut rcond, self.unwrap_unique()) };
        (Value::from(ret), rcond)
    }

    #[cfg(feature = "v2_2")]
    pub fn lcurve(
        &mut self,
        reg_param: &mut VectorF64,
        rho: &mut VectorF64,
        eta: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multilarge_linear_lcurve(
                reg_param.unwrap_unique(),
                rho.unwrap_unique(),
                eta.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn wstdform1(
        &mut self,
        L: &VectorF64,
        X: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multilarge_linear_wstdform1(
                L.unwrap_shared(),
                X.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn stdform1(
        &mut self,
        L: &VectorF64,
        X: &MatrixF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multilarge_linear_stdform1(
                L.unwrap_shared(),
                X.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn wstdform2(
        &mut self,
        LQR: &MatrixF64,
        Ltau: &VectorF64,
        X: &MatrixF64,
        w: &VectorF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multilarge_linear_wstdform2(
                LQR.unwrap_shared(),
                Ltau.unwrap_shared(),
                X.unwrap_shared(),
                w.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn stdform2(
        &mut self,
        LQR: &MatrixF64,
        Ltau: &VectorF64,
        X: &MatrixF64,
        y: &VectorF64,
        Xs: &mut MatrixF64,
        ys: &mut VectorF64,
    ) -> Value {
        unsafe {
            Value::from(sys::gsl_multilarge_linear_stdform2(
                LQR.unwrap_shared(),
                Ltau.unwrap_shared(),
                X.unwrap_shared(),
                y.unwrap_shared(),
                Xs.unwrap_unique(),
                ys.unwrap_unique(),
                self.unwrap_unique(),
            ))
        }
    }

    pub fn genform1(&mut self, L: &VectorF64, cs: &VectorF64, c: &mut VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_multilarge_linear_genform1(
                L.unwrap_shared(),
                cs.unwrap_shared(),
                c.unwrap_unique(),
                self.unwrap_unique(),
            )
        })
    }

    pub fn genform2(
        &mut self,
        LQR: &MatrixF64,
        Ltau: &VectorF64,
        cs: &VectorF64,
        c: &mut VectorF64,
    ) -> Value {
        Value::from(unsafe {
            sys::gsl_multilarge_linear_genform2(
                LQR.unwrap_shared(),
                Ltau.unwrap_shared(),
                cs.unwrap_shared(),
                c.unwrap_unique(),
                self.unwrap_unique(),
            )
        })
    }

    #[cfg(feature = "v2_7")]
    pub fn matrix<F: FnOnce(&MatrixF64)>(&self, f: F) {
        f(&MatrixF64::soft_wrap(unsafe {
            sys::gsl_multilarge_linear_matrix_ptr(self.unwrap_shared()) as _
        }))
    }

    #[cfg(feature = "v2_7")]
    pub fn rhs<F: FnOnce(&VectorF64)>(&self, f: F) {
        f(&VectorF64::soft_wrap(unsafe {
            sys::gsl_multilarge_linear_rhs_ptr(self.unwrap_shared()) as _
        }))
    }
}
