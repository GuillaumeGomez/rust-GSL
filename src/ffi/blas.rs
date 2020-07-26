//! BLAS and CBLAS support

use cblas;

pub type CBLAS_INDEX = c_uint;
pub type CBLAS_INDEX_t = CBLAS_INDEX;
pub type CBLAS_TRANSPOSE_t = cblas::Transpose;
pub type CBLAS_UPLO_t = cblas::Uplo;
pub type CBLAS_DIAG_t = cblas::Diag;
pub type CBLAS_SIDE_t = cblas::Side;

use libc::{c_double, c_float, c_int, c_uint, c_void};

use super::{
    gsl_complex, gsl_complex_float, gsl_matrix, gsl_matrix_complex, gsl_matrix_complex_float,
    gsl_matrix_float, gsl_vector, gsl_vector_complex, gsl_vector_complex_float, gsl_vector_float,
};

extern "C" {
    // Level 1 CBLAS functions
    pub fn cblas_sdsdot(
        N: c_int,
        alpha: c_float,
        x: *const c_float,
        incx: c_int,
        y: *const c_float,
        incy: c_int,
    ) -> c_float;
    pub fn cblas_dsdot(
        N: c_int,
        x: *const c_float,
        incx: c_int,
        y: *const c_float,
        incy: c_int,
    ) -> c_double;
    pub fn cblas_sdot(
        N: c_int,
        x: *const c_float,
        incx: c_int,
        y: *const c_float,
        incy: c_int,
    ) -> c_float;
    pub fn cblas_ddot(
        N: c_int,
        x: *const c_float,
        incx: c_int,
        y: *const c_float,
        incy: c_int,
    ) -> c_double;
    pub fn cblas_cdotu_sub(
        N: c_int,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        dotu: *mut c_void,
    );
    pub fn cblas_cdotc_sub(
        N: c_int,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        dotc: *mut c_void,
    );
    pub fn cblas_zdotu_sub(
        N: c_int,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        dotu: *mut c_void,
    );
    pub fn cblas_zdotc_sub(
        N: c_int,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        dotc: *mut c_void,
    );
    pub fn cblas_snrm2(N: c_int, x: *const c_float, incx: c_int) -> c_float;
    pub fn cblas_sasum(N: c_int, x: *const c_float, incx: c_int) -> c_float;
    pub fn cblas_dnrm2(N: c_int, x: *const c_double, incx: c_int) -> c_double;
    pub fn cblas_dasum(N: c_int, x: *const c_double, incx: c_int) -> c_double;
    pub fn cblas_scnrm2(N: c_int, x: *const c_void, incx: c_int) -> c_float;
    pub fn cblas_scasum(N: c_int, x: *const c_void, incx: c_int) -> c_float;
    pub fn cblas_dznrm2(N: c_int, x: *const c_void, incx: c_int) -> c_double;
    pub fn cblas_dzasum(N: c_int, x: *const c_void, incx: c_int) -> c_double;
    pub fn cblas_isamax(N: c_int, x: *const c_float, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_idamax(N: c_int, x: *const c_double, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_icamax(N: c_int, x: *const c_void, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_izamax(N: c_int, x: *const c_void, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_sswap(N: c_int, x: *mut c_float, incx: c_int, y: *mut c_float, incy: c_int);
    pub fn cblas_scopy(N: c_int, x: *const c_float, incx: c_int, y: *mut c_float, incy: c_int);
    pub fn cblas_saxpy(
        N: c_int,
        alpha: c_float,
        x: *const c_float,
        incx: c_int,
        y: *mut c_float,
        incy: c_int,
    );
    pub fn cblas_dswap(N: c_int, x: *mut c_double, incx: c_int, y: *mut c_double, incy: c_int);
    pub fn cblas_dcopy(N: c_int, x: *const c_double, incx: c_int, y: *mut c_double, incy: c_int);
    pub fn cblas_daxpy(
        N: c_int,
        alpha: c_double,
        x: *const c_double,
        incx: c_int,
        y: *mut c_double,
        incy: c_int,
    );
    pub fn cblas_cswap(N: c_int, x: *mut c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_ccopy(N: c_int, x: *const c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_caxpy(
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zswap(N: c_int, x: *mut c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_zcopy(N: c_int, x: *const c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_zaxpy(
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_srotg(a: *mut c_float, b: *mut c_float, c: *mut c_float, s: *mut c_float);
    pub fn cblas_srotmg(
        d1: *mut c_float,
        d2: *mut c_float,
        b1: *mut c_float,
        b2: *const c_float,
        P: *mut c_float,
    );
    pub fn cblas_srot(
        N: c_int,
        x: *mut c_float,
        incx: c_int,
        y: *mut c_float,
        incy: c_int,
        c: c_float,
        s: c_float,
    );
    pub fn cblas_srotm(
        N: c_int,
        x: *mut c_float,
        incx: c_int,
        y: *mut c_float,
        incy: c_int,
        P: *const c_float,
    );
    pub fn cblas_drotg(a: *mut c_double, b: *mut c_double, c: *mut c_double, s: *mut c_double);
    pub fn cblas_drotmg(
        d1: *mut c_double,
        d2: *mut c_double,
        b1: *mut c_double,
        b2: *const c_double,
        P: *mut c_double,
    );
    pub fn cblas_drot(
        N: c_int,
        x: *mut c_double,
        incx: c_int,
        y: *mut c_double,
        incy: c_int,
        c: c_double,
        s: c_double,
    );
    pub fn cblas_drotm(
        N: c_int,
        x: *mut c_double,
        incx: c_int,
        y: *mut c_double,
        incy: c_int,
        P: *const c_double,
    );
    pub fn cblas_sscal(N: c_int, alpha: c_float, x: *mut c_float, incx: c_int);
    pub fn cblas_dscal(N: c_int, alpha: c_double, x: *mut c_double, incx: c_int);
    pub fn cblas_cscal(N: c_int, alpha: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_zscal(N: c_int, alpha: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_csscal(N: c_int, alpha: c_float, x: *mut c_void, incx: c_int);
    pub fn cblas_zdscal(N: c_int, alpha: c_double, x: *mut c_void, incx: c_int);

    // Level 2 CBLAS functions
    pub fn cblas_sgemv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        x: *const c_float,
        incx: c_int,
        beta: c_float,
        y: *mut c_float,
        incy: c_int,
    );
    pub fn cblas_sgbmv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        KL: c_int,
        KU: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        x: *const c_float,
        incx: c_int,
        beta: c_float,
        y: *mut c_float,
        incy: c_int,
    );
    pub fn cblas_strmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_float,
        lda: c_int,
        x: *mut c_float,
        incx: c_int,
    );
    pub fn cblas_stbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_float,
        lda: c_int,
        x: *mut c_float,
        incx: c_int,
    );
    pub fn cblas_stpmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_float,
        x: *mut c_float,
        incx: c_int,
    );
    pub fn cblas_strsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_float,
        lda: c_int,
        x: *mut c_float,
        incx: c_int,
    );
    pub fn cblas_stbsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_float,
        lda: c_int,
        x: *mut c_float,
        incx: c_int,
    );
    pub fn cblas_stpsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_float,
        x: *mut c_float,
        incx: c_int,
    );
    pub fn cblas_dgemv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        x: *const c_double,
        incx: c_int,
        beta: c_double,
        y: *mut c_double,
        incy: c_int,
    );
    pub fn cblas_dgbmv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        KL: c_int,
        KU: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        x: *const c_double,
        incx: c_int,
        beta: c_double,
        y: *mut c_double,
        incy: c_int,
    );
    pub fn cblas_dtrmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_double,
        lda: c_int,
        x: *mut c_double,
        incx: c_int,
    );
    pub fn cblas_dtbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_double,
        lda: c_int,
        x: *mut c_double,
        incx: c_int,
    );
    pub fn cblas_dtpmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_double,
        x: *mut c_double,
        incx: c_int,
    );
    pub fn cblas_dtrsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_double,
        lda: c_int,
        x: *mut c_double,
        incx: c_int,
    );
    pub fn cblas_dtbsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_double,
        lda: c_int,
        x: *mut c_double,
        incx: c_int,
    );
    pub fn cblas_dtpsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_double,
        x: *mut c_double,
        incx: c_int,
    );
    pub fn cblas_cgemv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_cgbmv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        KL: c_int,
        KU: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_ctrmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ctbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ctpmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_void,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ctrsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ctbsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ctpsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_void,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_zgemv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zgbmv(
        order: cblas::Order,
        transA: cblas::Transpose,
        M: c_int,
        N: c_int,
        KL: c_int,
        KU: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_ztrmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ztbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ztpmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_void,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ztrsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ztbsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        K: c_int,
        A: *const c_void,
        lda: c_int,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ztpsv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        N: c_int,
        Ap: *const c_void,
        x: *mut c_void,
        incx: c_int,
    );
    pub fn cblas_ssymv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        x: *const c_float,
        incx: c_int,
        beta: c_float,
        y: *mut c_float,
        incy: c_int,
    );
    pub fn cblas_ssbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        K: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        x: *const c_float,
        incx: c_int,
        beta: c_float,
        y: *mut c_float,
        incy: c_int,
    );
    pub fn cblas_sspmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_float,
        Ap: *const c_float,
        x: *const c_float,
        incx: c_int,
        beta: c_float,
        y: *mut c_float,
        incy: c_int,
    );
    pub fn cblas_sger(
        order: cblas::Order,
        M: c_int,
        N: c_int,
        alpha: c_float,
        x: *const c_float,
        incx: c_int,
        y: *const c_float,
        incy: c_int,
        A: *mut c_float,
        lda: c_int,
    );
    pub fn cblas_ssyr(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_float,
        x: *const c_float,
        incx: c_int,
        A: *mut c_float,
        lda: c_int,
    );
    pub fn cblas_sspr(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_float,
        x: *const c_float,
        incx: c_int,
        Ap: *mut c_float,
    );
    pub fn cblas_ssyr2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_float,
        x: *const c_float,
        incx: c_int,
        y: *const c_float,
        incy: c_int,
        A: *mut c_float,
        lda: c_int,
    );
    pub fn cblas_sspr2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_float,
        x: *const c_float,
        incx: c_int,
        y: *const c_float,
        incy: c_int,
        A: *mut c_float,
    );
    pub fn cblas_dsymv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        x: *const c_double,
        incx: c_int,
        beta: c_double,
        y: *mut c_double,
        incy: c_int,
    );
    pub fn cblas_dsbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        K: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        x: *const c_double,
        incx: c_int,
        beta: c_double,
        y: *mut c_double,
        incy: c_int,
    );
    pub fn cblas_dspmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_double,
        Ap: *const c_double,
        x: *const c_double,
        incx: c_int,
        beta: c_double,
        y: *mut c_double,
        incy: c_int,
    );
    pub fn cblas_dger(
        order: cblas::Order,
        M: c_int,
        N: c_int,
        alpha: c_double,
        x: *const c_double,
        incx: c_int,
        y: *const c_double,
        incy: c_int,
        A: *mut c_double,
        lda: c_int,
    );
    pub fn cblas_dsyr(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_double,
        x: *const c_double,
        incx: c_int,
        A: *mut c_double,
        lda: c_int,
    );
    pub fn cblas_dspr(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_double,
        x: *const c_double,
        incx: c_int,
        Ap: *mut c_double,
    );
    pub fn cblas_dsyr2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_double,
        x: *const c_double,
        incx: c_int,
        y: *const c_double,
        incy: c_int,
        A: *mut c_double,
        lda: c_int,
    );
    pub fn cblas_dspr2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: c_double,
        x: *const c_double,
        incx: c_int,
        y: *const c_double,
        incy: c_int,
        A: *mut c_double,
    );
    pub fn cblas_chemv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_chbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_chpmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        Ap: *const c_void,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_csymv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_csbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_cspmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        Ap: *const c_void,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_cgeru(
        order: cblas::Order,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_cgerc(
        order: cblas::Order,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_cher(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_chpr(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        Ap: *mut c_void,
    );
    pub fn cblas_cher2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_chpr2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        Ap: *mut c_void,
    );
    pub fn cblas_zhemv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zhbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zhpmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        Ap: *const c_void,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zsymv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zsbmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zspmv(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        Ap: *const c_void,
        x: *const c_void,
        incx: c_int,
        beta: *const c_void,
        y: *mut c_void,
        incy: c_int,
    );
    pub fn cblas_zgeru(
        order: cblas::Order,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_zgerc(
        order: cblas::Order,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_zher(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_zhpr(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        Ap: *mut c_void,
    );
    pub fn cblas_zher2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        A: *mut c_void,
        lda: c_int,
    );
    pub fn cblas_zhpr2(
        order: cblas::Order,
        uplo: cblas::Uplo,
        N: c_int,
        alpha: *const c_void,
        x: *const c_void,
        incx: c_int,
        y: *const c_void,
        incy: c_int,
        Ap: *mut c_void,
    );

    // Level 3 CBLAS functions
    pub fn cblas_sgemm(
        order: cblas::Order,
        transA: cblas::Transpose,
        transB: cblas::Transpose,
        M: c_int,
        N: c_int,
        K: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        B: *const c_float,
        ldb: c_int,
        beta: c_float,
        C: *mut c_float,
        ldc: c_int,
    );
    pub fn cblas_ssymm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        M: c_int,
        N: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        B: *const c_float,
        ldb: c_int,
        beta: c_float,
        C: *mut c_float,
        ldc: c_int,
    );
    pub fn cblas_ssyrk(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        beta: c_float,
        C: *mut c_float,
        ldc: c_int,
    );
    pub fn cblas_ssyr2k(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        B: *const c_float,
        ldb: c_int,
        beta: c_float,
        C: *mut c_float,
        ldc: c_int,
    );
    pub fn cblas_strmm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        B: *mut c_float,
        ldb: c_int,
    );
    pub fn cblas_strsm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: c_float,
        A: *const c_float,
        lda: c_int,
        B: *mut c_float,
        ldb: c_int,
    );
    pub fn cblas_dgemm(
        order: cblas::Order,
        transA: cblas::Transpose,
        transB: cblas::Transpose,
        M: c_int,
        N: c_int,
        K: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        B: *const c_double,
        ldb: c_int,
        beta: c_double,
        C: *mut c_double,
        ldc: c_int,
    );
    pub fn cblas_dsymm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        M: c_int,
        N: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        B: *const c_double,
        ldb: c_int,
        beta: c_double,
        C: *mut c_double,
        ldc: c_int,
    );
    pub fn cblas_dsyrk(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        beta: c_double,
        C: *mut c_double,
        ldc: c_int,
    );
    pub fn cblas_dsyr2k(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        B: *const c_double,
        ldb: c_int,
        beta: c_double,
        C: *mut c_double,
        ldc: c_int,
    );
    pub fn cblas_dtrmm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        B: *mut c_double,
        ldb: c_int,
    );
    pub fn cblas_dtrsm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: c_double,
        A: *const c_double,
        lda: c_int,
        B: *mut c_double,
        ldb: c_int,
    );
    pub fn cblas_cgemm(
        order: cblas::Order,
        transA: cblas::Transpose,
        transB: cblas::Transpose,
        M: c_int,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_csymm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_csyrk(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_csyr2k(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_ctrmm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *mut c_void,
        ldb: c_int,
    );
    pub fn cblas_ctrsm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *mut c_void,
        ldb: c_int,
    );
    pub fn cblas_zgemm(
        order: cblas::Order,
        transA: cblas::Transpose,
        transB: cblas::Transpose,
        M: c_int,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_zsymm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_zsyrk(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_zsyr2k(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_ztrmm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *mut c_void,
        ldb: c_int,
    );
    pub fn cblas_ztrsm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        transA: cblas::Transpose,
        diag: cblas::Diag,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *mut c_void,
        ldb: c_int,
    );
    pub fn cblas_chemm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_cherk(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_cher2k(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_zhemm(
        order: cblas::Order,
        side: cblas::Side,
        uplo: cblas::Uplo,
        M: c_int,
        N: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_zherk(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: c_double,
        A: *const c_void,
        lda: c_int,
        beta: c_double,
        C: *mut c_void,
        ldc: c_int,
    );
    pub fn cblas_zher2k(
        order: cblas::Order,
        uplo: cblas::Uplo,
        trans: cblas::Transpose,
        N: c_int,
        K: c_int,
        alpha: *const c_void,
        A: *const c_void,
        lda: c_int,
        B: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        C: *mut c_void,
        ldc: c_int,
    );
    //to bind later
    //pub fn cblas_xerbla(p: c_int, rout: *const c_char, form: *const c_char, ...);

    // Level 1 BLAS functions
    pub fn gsl_blas_sdsdot(
        alpha: c_float,
        x: *const gsl_vector_float,
        y: *const gsl_vector_float,
        result: *mut c_float,
    ) -> c_int;
    pub fn gsl_blas_sdot(
        x: *const gsl_vector_float,
        y: *const gsl_vector_float,
        result: *mut c_float,
    ) -> c_int;
    pub fn gsl_blas_dsdot(
        x: *const gsl_vector_float,
        y: *const gsl_vector_float,
        result: *mut c_double,
    ) -> c_int;
    pub fn gsl_blas_ddot(
        x: *const gsl_vector,
        y: *const gsl_vector,
        result: *mut c_double,
    ) -> c_int;
    pub fn gsl_blas_cdotu(
        x: *const gsl_vector_complex_float,
        y: *const gsl_vector_complex_float,
        dotu: *mut gsl_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zdotu(
        x: *const gsl_vector_complex,
        y: *const gsl_vector_complex,
        dotu: *mut gsl_complex,
    ) -> c_int;
    pub fn gsl_blas_cdotc(
        x: *const gsl_vector_complex_float,
        y: *const gsl_vector_complex_float,
        dotc: *mut gsl_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zdotc(
        x: *const gsl_vector_complex,
        y: *const gsl_vector_complex,
        dotc: *mut gsl_complex,
    ) -> c_int;
    pub fn gsl_blas_snrm2(x: *const gsl_vector_float) -> c_float;
    pub fn gsl_blas_dnrm2(x: *const gsl_vector) -> c_double;
    pub fn gsl_blas_scnrm2(x: *const gsl_vector_complex_float) -> c_float;
    pub fn gsl_blas_dznrm2(x: *const gsl_vector_complex) -> c_double;
    pub fn gsl_blas_sasum(x: *const gsl_vector_float) -> c_float;
    pub fn gsl_blas_dasum(x: *const gsl_vector) -> c_double;
    pub fn gsl_blas_scasum(x: *const gsl_vector_complex_float) -> c_float;
    pub fn gsl_blas_dzasum(x: *const gsl_vector_complex) -> c_double;
    pub fn gsl_blas_isamax(x: *const gsl_vector_float) -> CBLAS_INDEX_t;
    pub fn gsl_blas_idamax(x: *const gsl_vector) -> CBLAS_INDEX_t;
    pub fn gsl_blas_icamax(x: *const gsl_vector_complex_float) -> CBLAS_INDEX_t;
    pub fn gsl_blas_izamax(x: *const gsl_vector_complex) -> CBLAS_INDEX_t;
    pub fn gsl_blas_sswap(x: *mut gsl_vector_float, y: *mut gsl_vector_float) -> c_int;
    pub fn gsl_blas_dswap(x: *mut gsl_vector, y: *mut gsl_vector) -> c_int;
    pub fn gsl_blas_cswap(
        x: *mut gsl_vector_complex_float,
        y: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zswap(x: *mut gsl_vector_complex, y: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_blas_scopy(x: *const gsl_vector_float, y: *mut gsl_vector_float) -> c_int;
    pub fn gsl_blas_dcopy(x: *const gsl_vector, y: *mut gsl_vector) -> c_int;
    pub fn gsl_blas_ccopy(
        x: *const gsl_vector_complex_float,
        y: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zcopy(x: *const gsl_vector_complex, y: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_blas_saxpy(
        alpha: c_float,
        x: *const gsl_vector_float,
        y: *mut gsl_vector_float,
    ) -> c_int;
    pub fn gsl_blas_daxpy(alpha: c_double, x: *const gsl_vector, y: *mut gsl_vector) -> c_int;
    pub fn gsl_blas_caxpy(
        alpha: gsl_complex_float,
        x: *const gsl_vector_complex_float,
        y: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zaxpy(
        alpha: gsl_complex,
        x: *const gsl_vector_complex,
        y: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_blas_sscal(alpha: c_float, x: *mut gsl_vector_float);
    pub fn gsl_blas_dscal(alpha: c_double, x: *mut gsl_vector);
    pub fn gsl_blas_cscal(alpha: gsl_complex_float, x: *mut gsl_vector_complex_float);
    pub fn gsl_blas_zscal(alpha: gsl_complex, x: *mut gsl_vector_complex);
    pub fn gsl_blas_csscal(alpha: c_float, x: *mut gsl_vector_complex_float);
    pub fn gsl_blas_zdscal(alpha: c_double, x: *mut gsl_vector_complex);
    pub fn gsl_blas_srotg(
        a: *mut c_float,
        b: *mut c_float,
        c: *mut c_float,
        d: *mut c_float,
    ) -> c_int;
    pub fn gsl_blas_drotg(
        a: *mut c_double,
        b: *mut c_double,
        c: *mut c_double,
        d: *mut c_double,
    ) -> c_int;
    pub fn gsl_blas_srot(
        a: *mut gsl_vector_float,
        b: *mut gsl_vector_float,
        c: c_float,
        d: c_float,
    ) -> c_int;
    pub fn gsl_blas_drot(a: *mut gsl_vector, b: *mut gsl_vector, c: c_double, d: c_double)
        -> c_int;
    pub fn gsl_blas_srotmg(
        d1: *mut c_float,
        d2: *mut c_float,
        b1: *mut c_float,
        b2: c_float,
        P: *mut c_float,
    ) -> c_int;
    pub fn gsl_blas_drotmg(
        d1: *mut c_double,
        d2: *mut c_double,
        b1: *mut c_double,
        b2: c_double,
        P: *mut c_double,
    ) -> c_int;
    pub fn gsl_blas_srotm(
        x: *mut gsl_vector_float,
        y: *mut gsl_vector_float,
        P: *mut c_float,
    ) -> c_int;
    pub fn gsl_blas_drotm(x: *mut gsl_vector, y: *mut gsl_vector, P: *mut c_double) -> c_int;

    // Level 2 BLAS functions
    pub fn gsl_blas_sgemv(
        transA: CBLAS_TRANSPOSE_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        x: *const gsl_vector_float,
        beta: c_float,
        y: *mut gsl_vector_float,
    ) -> c_int;
    pub fn gsl_blas_dgemv(
        transA: CBLAS_TRANSPOSE_t,
        alpha: c_double,
        A: *const gsl_matrix,
        x: *const gsl_vector,
        beta: c_double,
        y: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_blas_cgemv(
        transA: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        x: *const gsl_vector_complex_float,
        beta: gsl_complex_float,
        y: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zgemv(
        transA: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        x: *const gsl_vector_complex,
        beta: gsl_complex,
        y: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_blas_strmv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix_float,
        x: *mut gsl_vector_float,
    ) -> c_int;
    pub fn gsl_blas_dtrmv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_blas_ctrmv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix_complex_float,
        x: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_blas_ztrmv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix_complex,
        x: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_blas_strsv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix_float,
        x: *mut gsl_vector_float,
    ) -> c_int;
    pub fn gsl_blas_dtrsv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix,
        x: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_blas_ctrsv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix_complex_float,
        x: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_blas_ztrsv(
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        A: *const gsl_matrix_complex,
        x: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_blas_ssymv(
        uplo: CBLAS_UPLO_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        x: *const gsl_vector_float,
        beta: c_float,
        y: *mut gsl_vector_float,
    ) -> c_int;
    pub fn gsl_blas_dsymv(
        uplo: CBLAS_UPLO_t,
        alpha: c_double,
        A: *const gsl_matrix,
        x: *const gsl_vector,
        beta: c_double,
        y: *mut gsl_vector,
    ) -> c_int;
    pub fn gsl_blas_chemv(
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        x: *const gsl_vector_complex_float,
        beta: gsl_complex_float,
        y: *mut gsl_vector_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zhemv(
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        x: *const gsl_vector_complex,
        beta: gsl_complex,
        y: *mut gsl_vector_complex,
    ) -> c_int;
    pub fn gsl_blas_sger(
        alpha: c_float,
        x: *const gsl_vector_float,
        y: *const gsl_vector_float,
        A: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dger(
        alpha: c_double,
        x: *const gsl_vector,
        y: *const gsl_vector,
        A: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_cgeru(
        alpha: gsl_complex_float,
        x: *const gsl_vector_complex_float,
        y: *const gsl_vector_complex_float,
        A: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zgeru(
        alpha: gsl_complex,
        x: *const gsl_vector_complex,
        y: *const gsl_vector_complex,
        A: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_cgerc(
        alpha: gsl_complex_float,
        x: *const gsl_vector_complex_float,
        y: *const gsl_vector_complex_float,
        A: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zgerc(
        alpha: gsl_complex,
        x: *const gsl_vector_complex,
        y: *const gsl_vector_complex,
        A: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_ssyr(
        uplo: CBLAS_UPLO_t,
        alpha: c_float,
        x: *const gsl_vector_float,
        A: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dsyr(
        uplo: CBLAS_UPLO_t,
        alpha: c_double,
        x: *const gsl_vector,
        A: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_cher(
        uplo: CBLAS_UPLO_t,
        alpha: c_float,
        x: *const gsl_vector_complex_float,
        A: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zher(
        uplo: CBLAS_UPLO_t,
        alpha: c_double,
        x: *const gsl_vector_complex,
        A: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_ssyr2(
        uplo: CBLAS_UPLO_t,
        alpha: c_float,
        x: *const gsl_vector_float,
        y: *const gsl_vector_float,
        A: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dsyr2(
        uplo: CBLAS_UPLO_t,
        alpha: c_double,
        x: *const gsl_vector,
        y: *const gsl_vector,
        A: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_cher2(
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex_float,
        x: *const gsl_vector_complex_float,
        y: *const gsl_vector_complex_float,
        A: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zher2(
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex,
        x: *const gsl_vector_complex,
        y: *const gsl_vector_complex,
        A: *mut gsl_matrix_complex,
    ) -> c_int;

    // Level 3 BLAS functions
    pub fn gsl_blas_sgemm(
        transA: CBLAS_TRANSPOSE_t,
        transB: CBLAS_TRANSPOSE_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        B: *const gsl_matrix_float,
        beta: c_float,
        C: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dgemm(
        transA: CBLAS_TRANSPOSE_t,
        transB: CBLAS_TRANSPOSE_t,
        alpha: c_double,
        A: *const gsl_matrix,
        B: *const gsl_matrix,
        beta: c_double,
        C: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_cgemm(
        transA: CBLAS_TRANSPOSE_t,
        transB: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float,
        beta: gsl_complex_float,
        C: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zgemm(
        transA: CBLAS_TRANSPOSE_t,
        transB: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex,
        beta: gsl_complex,
        C: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_ssymm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        B: *const gsl_matrix_float,
        beta: c_float,
        C: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dsymm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        alpha: c_double,
        A: *const gsl_matrix,
        B: *const gsl_matrix,
        beta: c_double,
        C: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_csymm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float,
        beta: gsl_complex_float,
        C: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zsymm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex,
        beta: gsl_complex,
        C: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_chemm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float,
        beta: gsl_complex_float,
        C: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zhemm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex,
        beta: gsl_complex,
        C: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_strmm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        B: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dtrmm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: c_double,
        A: *const gsl_matrix,
        B: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_ctrmm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        B: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_ztrmm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        B: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_strsm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        B: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dtrsm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: c_double,
        A: *const gsl_matrix,
        B: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_ctrsm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        B: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_ztrsm(
        side: CBLAS_SIDE_t,
        uplo: CBLAS_UPLO_t,
        transA: CBLAS_TRANSPOSE_t,
        diag: CBLAS_DIAG_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        B: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_ssyrk(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        beta: c_float,
        C: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dsyrk(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: c_double,
        A: *const gsl_matrix,
        beta: c_double,
        C: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_csyrk(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        beta: gsl_complex_float,
        C: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zsyrk(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        beta: gsl_complex,
        C: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_cherk(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: c_float,
        A: *const gsl_matrix_complex_float,
        beta: c_float,
        C: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zherk(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: c_double,
        A: *const gsl_matrix_complex,
        beta: c_double,
        C: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_ssyr2k(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: c_float,
        A: *const gsl_matrix_float,
        B: *const gsl_matrix_float,
        beta: c_float,
        C: *mut gsl_matrix_float,
    ) -> c_int;
    pub fn gsl_blas_dsyr2k(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: c_double,
        A: *const gsl_matrix,
        B: *const gsl_matrix,
        beta: c_double,
        C: *mut gsl_matrix,
    ) -> c_int;
    pub fn gsl_blas_csyr2k(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float,
        beta: gsl_complex_float,
        C: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zsyr2k(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex,
        beta: gsl_complex,
        C: *mut gsl_matrix_complex,
    ) -> c_int;
    pub fn gsl_blas_cher2k(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float,
        beta: c_float,
        C: *mut gsl_matrix_complex_float,
    ) -> c_int;
    pub fn gsl_blas_zher2k(
        uplo: CBLAS_UPLO_t,
        trans: CBLAS_TRANSPOSE_t,
        alpha: gsl_complex,
        A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex,
        beta: c_double,
        C: *mut gsl_matrix_complex,
    ) -> c_int;
}
