/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

pub mod Cblas {
    pub mod Level1 {
        pub fn sdsdot(N: i32, alpha: f32, x: &[f32], incx: i32, y: &[f32], incy: i32) -> f32 {
            unsafe { ::ffi::cblas_sdsdot(N, alpha, x.as_ptr(), incx, y.as_ptr(), incy) }
        }

        pub fn dsdot(N: i32, x: &[f32], incx: i32, y: &[f32], incy: i32) -> f64 {
            unsafe { ::ffi::cblas_dsdot(N, x.as_ptr(), incx, y.as_ptr(), incy) }
        }

        pub fn sdot(N: i32, x: &[f32], incx: i32, y: &[f32], incy: i32) -> f32 {
            unsafe { ::ffi::cblas_sdot(N, x.as_ptr(), incx, y.as_ptr(), incy) }
        }

        pub fn ddot(N: i32, x: &[f32], incx: i32, y: &[f32], incy: i32) -> f64 {
            unsafe { ::ffi::cblas_ddot(N, x.as_ptr(), incx, y.as_ptr(), incy) }
        }

        pub fn cdotu_sub<T>(N: i32, x: &[T], incx: i32, y: &[T], incy: i32, dotu: &mut [T]) {
            unsafe { ::ffi::cblas_cdotu_sub(N, x.as_ptr() as *const ::libc::c_void, incx, y.as_ptr() as *const ::libc::c_void,
                incy, dotu.as_mut_ptr() as *mut ::libc::c_void) }
        }

        pub fn cdotc_sub<T>(N: i32, x: &[T], incx: i32, y: &[T], incy: i32, dotc: &mut [T]) {
            unsafe { ::ffi::cblas_cdotc_sub(N, x.as_ptr() as *const ::libc::c_void, incx, y.as_ptr() as *const ::libc::c_void,
                incy, dotc.as_mut_ptr() as *mut ::libc::c_void) }
        }

        pub fn zdotu_sub<T>(N: i32, x: &[T], incx: i32, y: &[T], incy: i32, dotu: &mut [T]) {
            unsafe { ::ffi::cblas_zdotu_sub(N, x.as_ptr() as *const ::libc::c_void, incx, y.as_ptr() as *const ::libc::c_void,
                incy, dotu.as_mut_ptr() as *mut ::libc::c_void) }
        }

        pub fn zdotc_sub<T>(N: i32, x: &[T], incx: i32, y: &[T], incy: i32, dotc: &mut [T]) {
            unsafe { ::ffi::cblas_zdotc_sub(N, x.as_ptr() as *const ::libc::c_void, incx, y.as_ptr() as *const ::libc::c_void,
                incy, dotc.as_mut_ptr() as *mut ::libc::c_void) }
        }

        pub fn snrm2(N: i32, x: &[f32], incx: i32) -> f32 {
            unsafe { ::ffi::cblas_snrm2(N, x.as_ptr(), incx) }
        }

        pub fn sasum(N: i32, x: &[f32], incx: i32) -> f32 {
            unsafe { ::ffi::cblas_sasum(N, x.as_ptr(), incx) }
        }

        pub fn dnrm2(N: i32, x: &[f64], incx: i32) -> f64 {
            unsafe { ::ffi::cblas_dnrm2(N, x.as_ptr(), incx) }
        }

        pub fn dasum(N: i32, x: &[f64], incx: i32) -> f64 {
            unsafe { ::ffi::cblas_dasum(N, x.as_ptr(), incx) }
        }

        pub fn scnrm2<T>(N: i32, x: &[T], incx: i32) -> f32 {
            unsafe { ::ffi::cblas_scnrm2(N, x.as_ptr() as *const ::libc::c_void, incx) }
        }

        pub fn scasum<T>(N: i32, x: &[T], incx: i32) -> f32 {
            unsafe { ::ffi::cblas_scasum(N, x.as_ptr() as *const ::libc::c_void, incx) }
        }

        pub fn dznrm2<T>(N: i32, x: &[T], incx: i32) -> f64 {
            unsafe { ::ffi::cblas_dznrm2(N, x.as_ptr() as *const ::libc::c_void, incx) }
        }

        pub fn dzasum<T>(N: i32, x: &[T], incx: i32) -> f64 {
            unsafe { ::ffi::cblas_dzasum(N, x.as_ptr() as *const ::libc::c_void, incx) }
        }

        pub fn isamax(N: i32, x: &[f32], incx: i32) -> ::types::CblasIndex {
            ::types::CblasIndex(unsafe { ::ffi::cblas_isamax(N, x.as_ptr(), incx) })
        }

        pub fn idamax(N: i32, x: &[f64], incx: i32) -> ::types::CblasIndex {
            ::types::CblasIndex(unsafe { ::ffi::cblas_idamax(N, x.as_ptr(), incx) })
        }

        pub fn icamax<T>(N: i32, x: &[T], incx: i32) -> ::types::CblasIndex {
            ::types::CblasIndex(unsafe { ::ffi::cblas_icamax(N, x.as_ptr() as *const ::libc::c_void, incx) })
        }

        pub fn izamax<T>(N: i32, x: &[T], incx: i32) -> ::types::CblasIndex {
            ::types::CblasIndex(unsafe { ::ffi::cblas_izamax(N, x.as_ptr() as *const ::libc::c_void, incx) })
        }

        pub fn sswap(N: i32, x: &mut [f32], incx: i32, y: &mut [f32], incy: i32) {
            unsafe { ::ffi::cblas_sswap(N, x.as_mut_ptr(), incx, y.as_mut_ptr(), incy) }
        }

        pub fn scopy(N: i32, x: &[f32], incx: i32, y: &mut [f32], incy: i32) {
            unsafe { ::ffi::cblas_scopy(N, x.as_ptr(), incx, y.as_mut_ptr(), incy) }
        }

        pub fn saxpy(N: i32, alpha: f32, x: &[f32], incx: i32, y: &mut [f32], incy: i32) {
            unsafe { ::ffi::cblas_saxpy(N, alpha, x.as_ptr(), incx, y.as_mut_ptr(), incy) }
        }

        pub fn dswap(N: i32, x: &mut [f64], incx: i32, y: &mut [f64], incy: i32) {
            unsafe { ::ffi::cblas_dswap(N, x.as_mut_ptr(), incx, y.as_mut_ptr(), incy) }
        }

        pub fn dcopy(N: i32, x: &[f64], incx: i32, y: &mut [f64], incy: i32) {
            unsafe { ::ffi::cblas_dcopy(N, x.as_ptr(), incx, y.as_mut_ptr(), incy) }
        }

        pub fn daxpy(N: i32, alpha: f64, x: &[f64], incx: i32, y: &mut [f64], incy: i32) {
            unsafe { ::ffi::cblas_daxpy(N, alpha, x.as_ptr(), incx, y.as_mut_ptr(), incy) }
        }

        pub fn cswap<T>(N: i32, x: &mut [T], incx: i32, y: &mut [T], incy: i32) {
            unsafe { ::ffi::cblas_cswap(N, x.as_mut_ptr() as *mut ::libc::c_void, incx, y.as_mut_ptr() as *mut ::libc::c_void, incy) }
        }

        pub fn ccopy<T>(N: i32, x: &[T], incx: i32, y: &mut [T], incy: i32) {
            unsafe { ::ffi::cblas_ccopy(N, x.as_ptr() as *const ::libc::c_void, incx, y.as_mut_ptr() as *mut ::libc::c_void, incy) }
        }

        pub fn caxpy<T>(N: i32, alpha: &[T], x: &[T], incx: i32, y: &mut [T], incy: i32) {
            unsafe { ::ffi::cblas_caxpy(N, alpha.as_ptr() as *const ::libc::c_void, x.as_ptr() as *const ::libc::c_void,
                incx, y.as_mut_ptr() as *mut ::libc::c_void, incy) }
        }

        pub fn zswap<T>(N: i32, x: &mut [T], incx: i32, y: &mut [T], incy: i32) {
            unsafe { ::ffi::cblas_zswap(N, x.as_mut_ptr() as *mut ::libc::c_void, incx, y.as_mut_ptr() as *mut ::libc::c_void, incy) }
        }

        pub fn zcopy<T>(N: i32, x: &[T], incx: i32, y: &mut [T], incy: i32) {
            unsafe { ::ffi::cblas_zcopy(N, x.as_ptr() as *const ::libc::c_void, incx, y.as_mut_ptr() as *mut ::libc::c_void, incy) }
        }

        pub fn zaxpy<T>(N: i32, alpha: &[T], x: &[T], incx: i32, y: &mut [T], incy: i32) {
            unsafe { ::ffi::cblas_zaxpy(N, alpha.as_ptr() as *const ::libc::c_void, x.as_ptr() as *const ::libc::c_void,
                incx, y.as_mut_ptr() as *mut ::libc::c_void, incy) }
        }

        pub fn srotg(a: &mut [f32], b: &mut [f32], c: &mut [f32], s: &mut [f32]) {
            unsafe { ::ffi::cblas_srotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), s.as_mut_ptr()) }
        }

        pub fn srotmg(d1: &mut [f32], d2: &mut [f32], b1: &mut [f32], b2: &[f32], P: &mut [f32]) {
            unsafe { ::ffi::cblas_srotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2.as_ptr(), P.as_mut_ptr()) }
        }

        pub fn srot(N: i32, x: &mut [f32], incx: i32, y: &mut [f32], incy: i32, c: f32, s: f32) {
            unsafe { ::ffi::cblas_srot(N, x.as_mut_ptr(), incx, y.as_mut_ptr(), incy, c, s) }
        }

        pub fn srotm(N: i32, x: &mut [f32], incx: i32, y: &mut [f32], incy: i32, p: &[f32]) {
            unsafe { ::ffi::cblas_srotm(N, x.as_mut_ptr(), incx, y.as_mut_ptr(), incy, p.as_ptr()) }
        }

        pub fn drotg(a: &mut [f64], b: &mut [f64], c: &mut [f64], s: &mut [f64]) {
            unsafe { ::ffi::cblas_drotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), s.as_mut_ptr()) }
        }

        pub fn drotmg(d1: &mut [f64], d2: &mut [f64], b1: &mut [f64], b2: &[f64], P: &mut [f64]) {
            unsafe { ::ffi::cblas_drotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2.as_ptr(), P.as_mut_ptr()) }
        }

        pub fn drot(N: i32, x: &mut [f64], incx: i32, y: &mut [f64], incy: i32, c: f64, s: f64) {
            unsafe { ::ffi::cblas_drot(N, x.as_mut_ptr(), incx, y.as_mut_ptr(), incy, c, s) }
        }

        pub fn drotm(N: i32, x: &mut [f64], incx: i32, y: &mut [f64], incy: i32, p: &[f64]) {
            unsafe { ::ffi::cblas_drotm(N, x.as_mut_ptr(), incx, y.as_mut_ptr(), incy, p.as_ptr()) }
        }

        pub fn sscal(N: i32, alpha: f32, x: &mut [f32], incx: i32) {
            unsafe { ::ffi::cblas_sscal(N, alpha, x.as_mut_ptr(), incx) }
        }

        pub fn dscal(N: i32, alpha: f64, x: &mut [f64], incx: i32) {
            unsafe { ::ffi::cblas_dscal(N, alpha, x.as_mut_ptr(), incx) }
        }

        pub fn cscal<T>(N: i32, alpha: &[T], x: &mut [T], incx: i32) {
            unsafe { ::ffi::cblas_cscal(N, alpha.as_ptr() as *const ::libc::c_void, x.as_mut_ptr() as *mut ::libc::c_void, incx) }
        }

        pub fn zscal<T>(N: i32, alpha: &[T], x: &mut [T], incx: i32) {
            unsafe { ::ffi::cblas_zscal(N, alpha.as_ptr() as *const ::libc::c_void, x.as_mut_ptr() as *mut ::libc::c_void, incx) }
        }

        pub fn csscal<T>(N: i32, alpha: f32, x: &mut [T], incx: i32) {
            unsafe { ::ffi::cblas_csscal(N, alpha, x.as_mut_ptr() as *mut ::libc::c_void, incx) }
        }

        pub fn zdscal<T>(N: i32, alpha: f64, x: &mut [T], incx: i32) {
            unsafe { ::ffi::cblas_zdscal(N, alpha, x.as_mut_ptr() as *mut ::libc::c_void, incx) }
        }
    }
}