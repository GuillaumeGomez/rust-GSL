//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub use sys::CBLAS_DIAG_t as Diag;
pub use sys::CBLAS_INDEX_t as Index;
pub use sys::CBLAS_ORDER_t as Order;
pub use sys::CBLAS_SIDE_t as Side;
pub use sys::CBLAS_TRANSPOSE_t as Transpose;
pub use sys::CBLAS_UPLO_t as Uplo;

pub mod level1 {
    use enums;
    use ffi;

    /// This function computes the sum \alpha + x^T y for the vectors x and y, returning the result in result.
    pub fn sdsdot(
        alpha: f32,
        x: &::types::VectorF32,
        y: &::types::VectorF32,
        result: &mut f32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sdsdot(
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                result,
            )
        })
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
    pub fn sdot(x: &::types::VectorF32, y: &::types::VectorF32, result: &mut f32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sdot(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                result,
            )
        })
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
    pub fn dsdot(x: &::types::VectorF32, y: &::types::VectorF32, result: &mut f64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsdot(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                result,
            )
        })
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
    pub fn ddot(x: &::types::VectorF64, y: &::types::VectorF64, result: &mut f64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ddot(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                result,
            )
        })
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning the result in dotu.
    pub fn cdotu(
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
        dotu: &mut ::types::ComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cdotu(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ::std::mem::transmute(dotu),
            )
        })
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning the result in dotu.
    pub fn zdotu(
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
        dotu: &mut ::types::ComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zdotu(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ::std::mem::transmute(dotu),
            )
        })
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y, returning the result in dotc.
    pub fn cdotc(
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
        dotc: &mut ::types::ComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cdotc(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ::std::mem::transmute(dotc),
            )
        })
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y, returning the result in dotc.
    pub fn zdotc(
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
        dotc: &mut ::types::ComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zdotc(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ::std::mem::transmute(dotc),
            )
        })
    }

    /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
    pub fn snrm2(x: &::types::VectorF32) -> f32 {
        unsafe { sys::gsl_blas_snrm2(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
    pub fn dnrm2(x: &::types::VectorF64) -> f64 {
        unsafe { sys::gsl_blas_dnrm2(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function computes the Euclidean norm of the complex vector x,
    ///
    /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
    pub fn scnrm2(x: &::types::VectorComplexF32) -> f32 {
        unsafe { sys::gsl_blas_scnrm2(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function computes the Euclidean norm of the complex vector x,
    ///
    /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
    pub fn dznrm2(x: &::types::VectorComplexF64) -> f64 {
        unsafe { sys::gsl_blas_dznrm2(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
    pub fn sasum(x: &::types::VectorF32) -> f32 {
        unsafe { sys::gsl_blas_sasum(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
    pub fn dasum(x: &::types::VectorF64) -> f64 {
        unsafe { sys::gsl_blas_dasum(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
    pub fn scasum(x: &::types::VectorComplexF32) -> f32 {
        unsafe { sys::gsl_blas_scasum(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
    pub fn dzasum(x: &::types::VectorComplexF64) -> f64 {
        unsafe { sys::gsl_blas_dzasum(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn isamax(x: &::types::VectorF32) -> u32 {
        unsafe { sys::gsl_blas_isamax(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn idamax(x: &::types::VectorF64) -> u32 {
        unsafe { sys::gsl_blas_idamax(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn icamax(x: &::types::VectorComplexF32) -> u32 {
        unsafe { sys::gsl_blas_icamax(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn izamax(x: &::types::VectorComplexF64) -> u32 {
        unsafe { sys::gsl_blas_izamax(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn sswap(x: &mut ::types::VectorF32, y: &mut ::types::VectorF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sswap(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn dswap(x: &mut ::types::VectorF64, y: &mut ::types::VectorF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dswap(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn cswap(
        x: &mut ::types::VectorComplexF32,
        y: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cswap(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn zswap(
        x: &mut ::types::VectorComplexF64,
        y: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zswap(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn scopy(x: &mut ::types::VectorF32, y: &mut ::types::VectorF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_scopy(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn dcopy(x: &mut ::types::VectorF64, y: &mut ::types::VectorF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dcopy(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn ccopy(
        x: &mut ::types::VectorComplexF32,
        y: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ccopy(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn zcopy(
        x: &mut ::types::VectorComplexF64,
        y: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zcopy(ffi::FFI::unwrap_unique(x), ffi::FFI::unwrap_unique(y))
        })
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn saxpy(alpha: f32, x: &::types::VectorF32, y: &mut ::types::VectorF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_saxpy(
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn daxpy(alpha: f64, x: &::types::VectorF64, y: &mut ::types::VectorF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_daxpy(
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn caxpy(
        alpha: &::types::ComplexF32,
        x: &::types::VectorComplexF32,
        y: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_caxpy(
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn zaxpy(
        alpha: &::types::ComplexF64,
        x: &::types::VectorComplexF64,
        y: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zaxpy(
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn sscal(alpha: f32, x: &mut ::types::VectorF32) {
        unsafe { sys::gsl_blas_sscal(alpha, ffi::FFI::unwrap_unique(x)) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn dscal(alpha: f64, x: &mut ::types::VectorF64) {
        unsafe { sys::gsl_blas_dscal(alpha, ffi::FFI::unwrap_unique(x)) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn cscal(alpha: &::types::ComplexF32, x: &mut ::types::VectorComplexF32) {
        unsafe {
            sys::gsl_blas_cscal(::std::mem::transmute(*alpha), ffi::FFI::unwrap_unique(x))
        }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn zscal(alpha: &::types::ComplexF64, x: &mut ::types::VectorComplexF64) {
        unsafe {
            sys::gsl_blas_zscal(::std::mem::transmute(*alpha), ffi::FFI::unwrap_unique(x))
        }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn csscal(alpha: f32, x: &mut ::types::VectorComplexF32) {
        unsafe { sys::gsl_blas_csscal(alpha, ffi::FFI::unwrap_unique(x)) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn zdscal(alpha: f64, x: &mut ::types::VectorComplexF64) {
        unsafe { sys::gsl_blas_zdscal(alpha, ffi::FFI::unwrap_unique(x)) }
    }

    /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
    ///
    /// [  c  s ] [ a ] = [ r ]
    ///
    /// [ -s  c ] [ b ]   [ 0 ]
    ///
    /// The variables a and b are overwritten by the routine.
    pub fn srotg(a: &mut [f32], b: &mut [f32], c: &mut [f32], d: &mut [f32]) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_srotg(
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                c.as_mut_ptr(),
                d.as_mut_ptr(),
            )
        })
    }

    /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
    ///
    /// [  c  s ] [ a ] = [ r ]
    ///
    /// [ -s  c ] [ b ]   [ 0 ]
    ///
    /// The variables a and b are overwritten by the routine.
    pub fn drotg(a: &mut [f64], b: &mut [f64], c: &mut [f64], d: &mut [f64]) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_drotg(
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                c.as_mut_ptr(),
                d.as_mut_ptr(),
            )
        })
    }

    /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
    pub fn srot(
        a: &mut ::types::VectorF32,
        b: &mut ::types::VectorF32,
        c: f32,
        d: f32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_srot(ffi::FFI::unwrap_unique(a), ffi::FFI::unwrap_unique(b), c, d)
        })
    }

    /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
    pub fn drot(
        a: &mut ::types::VectorF64,
        b: &mut ::types::VectorF64,
        c: f64,
        d: f64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_drot(ffi::FFI::unwrap_unique(a), ffi::FFI::unwrap_unique(b), c, d)
        })
    }

    /// This function computes a modified Givens transformation.
    /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
    pub fn srotmg(
        d1: &mut [f32],
        d2: &mut [f32],
        b1: &mut [f32],
        b2: f32,
        P: &mut [f32],
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_srotmg(
                d1.as_mut_ptr(),
                d2.as_mut_ptr(),
                b1.as_mut_ptr(),
                b2,
                P.as_mut_ptr(),
            )
        })
    }

    /// This function computes a modified Givens transformation.
    /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
    pub fn drotmg(
        d1: &mut [f64],
        d2: &mut [f64],
        b1: &mut [f64],
        b2: f64,
        P: &mut [f64],
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_drotmg(
                d1.as_mut_ptr(),
                d2.as_mut_ptr(),
                b1.as_mut_ptr(),
                b2,
                P.as_mut_ptr(),
            )
        })
    }

    /// This function applies a modified Givens transformation.
    pub fn srotm(
        x: &mut ::types::VectorF32,
        y: &mut ::types::VectorF32,
        P: &mut [f32],
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_srotm(
                ffi::FFI::unwrap_unique(x),
                ffi::FFI::unwrap_unique(y),
                P.as_mut_ptr(),
            )
        })
    }

    /// This function applies a modified Givens transformation.
    pub fn drotm(
        x: &mut ::types::VectorF64,
        y: &mut ::types::VectorF64,
        P: &mut [f64],
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_drotm(
                ffi::FFI::unwrap_unique(x),
                ffi::FFI::unwrap_unique(y),
                P.as_mut_ptr(),
            )
        })
    }
}

pub mod level2 {
    use enums;
    use ffi;

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn sgemv(
        transA: ::blas::Transpose,
        alpha: f32,
        A: &::types::MatrixF32,
        x: &::types::VectorF32,
        beta: f32,
        y: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sgemv(
                transA,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                beta,
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn dgemv(
        transA: ::blas::Transpose,
        alpha: f64,
        A: &::types::MatrixF64,
        x: &::types::VectorF64,
        beta: f64,
        y: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dgemv(
                transA,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                beta,
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn cgemv(
        transA: ::blas::Transpose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        x: &::types::VectorComplexF32,
        beta: &::types::ComplexF32,
        y: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cgemv(
                transA,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn zgemv(
        transA: ::blas::Transpose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        x: &::types::VectorComplexF64,
        beta: &::types::ComplexF64,
        y: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zgemv(
                transA,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strmv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixF32,
        x: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strmv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrmv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixF64,
        x: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrmv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrmv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixComplexF32,
        x: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrmv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrmv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixComplexF64,
        x: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrmv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strsv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixF32,
        x: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strsv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrsv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixF64,
        x: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrsv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrsv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixComplexF32,
        x: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrsv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrsv(
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        A: &::types::MatrixComplexF64,
        x: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrsv(
                uplo,
                transA,
                diag,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssymv(
        uplo: ::blas::Uplo,
        alpha: f32,
        A: &::types::MatrixF32,
        x: &::types::VectorF32,
        beta: f32,
        y: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssymv(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                beta,
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsymv(
        uplo: ::blas::Uplo,
        alpha: f64,
        A: &::types::MatrixF64,
        x: &::types::VectorF64,
        beta: f64,
        y: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsymv(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                beta,
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored. When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically assumed to be zero and are not referenced.
    pub fn chemv(
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        x: &::types::VectorComplexF32,
        beta: &::types::ComplexF32,
        y: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_chemv(
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored. When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically assumed to be zero and are not referenced.
    pub fn zhemv(
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        x: &::types::VectorComplexF64,
        beta: &::types::ComplexF64,
        y: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zhemv(
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(x),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(y),
            )
        })
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn sger(
        alpha: f32,
        x: &::types::VectorF32,
        y: &::types::VectorF32,
        A: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sger(
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn dger(
        alpha: f64,
        x: &::types::VectorF64,
        y: &::types::VectorF64,
        A: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dger(
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn cgeru(
        alpha: &::types::ComplexF32,
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
        A: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cgeru(
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn zgeru(
        alpha: &::types::ComplexF64,
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
        A: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zgeru(
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the conjugate rank-1 update A = \alpha x y^H + A of the matrix A.
    pub fn cgerc(
        alpha: &::types::ComplexF32,
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
        A: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cgerc(
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the conjugate rank-1 update A = \alpha x y^H + A of the matrix A.
    pub fn zgerc(
        alpha: &::types::ComplexF64,
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
        A: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zgerc(
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the symmetric rank-1 update A = \alpha x x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssyr(
        uplo: ::blas::Uplo,
        alpha: f32,
        x: &::types::VectorF32,
        A: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyr(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the symmetric rank-1 update A = \alpha x x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsyr(
        uplo: ::blas::Uplo,
        alpha: f64,
        x: &::types::VectorF64,
        A: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyr(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// These functions compute the hermitian rank-1 update A = \alpha x x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cher(
        uplo: ::blas::Uplo,
        alpha: f32,
        x: &::types::VectorComplexF32,
        A: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cher(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// These functions compute the hermitian rank-1 update A = \alpha x x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zher(
        uplo: ::blas::Uplo,
        alpha: f64,
        x: &::types::VectorComplexF64,
        A: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zher(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// These functions compute the symmetric rank-2 update A = \alpha x y^T + \alpha y x^T + A of the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssyr2(
        uplo: ::blas::Uplo,
        alpha: f32,
        x: &::types::VectorF32,
        y: &::types::VectorF32,
        A: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyr2(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// These functions compute the symmetric rank-2 update A = \alpha x y^T + \alpha y x^T + A of the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsyr2(
        uplo: ::blas::Uplo,
        alpha: f64,
        x: &::types::VectorF64,
        y: &::types::VectorF64,
        A: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyr2(
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// These functions compute the hermitian rank-2 update A = \alpha x y^H + \alpha^* y x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cher2(
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF32,
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
        A: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cher2(
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// These functions compute the hermitian rank-2 update A = \alpha x y^H + \alpha^* y x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zher2(
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF64,
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
        A: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zher2(
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }
}

pub mod level3 {
    use enums;
    use ffi;

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn sgemm(
        transA: ::blas::Transpose,
        transB: ::blas::Transpose,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sgemm(
                transA,
                transB,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn dgemm(
        transA: ::blas::Transpose,
        transB: ::blas::Transpose,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dgemm(
                transA,
                transB,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn cgemm(
        transA: ::blas::Transpose,
        transB: ::blas::Transpose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cgemm(
                transA,
                transB,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn zgemm(
        transA: ::blas::Transpose,
        transB: ::blas::Transpose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zgemm(
                transA,
                transB,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssymm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssymm(
                side,
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsymm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsymm(
                side,
                uplo,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn csymm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_csymm(
                side,
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn zsymm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zsymm(
                side,
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is Left and C = \alpha B A + \beta C for Side is Right, where the matrix A is hermitian.
    /// When Uplo is Upper then the upper triangle and diagonal of A are used, and when Uplo is Lower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn chemm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_chemm(
                side,
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is hermitian.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zhemm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zhemm(
                side,
                uplo,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strmm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strmm(
                side,
                uplo,
                transA,
                diag,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrmm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrmm(
                side,
                uplo,
                transA,
                diag,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrmm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrmm(
                side,
                uplo,
                transA,
                diag,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrmm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrmm(
                side,
                uplo,
                transA,
                diag,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strsm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strsm(
                side,
                uplo,
                transA,
                diag,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrsm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrsm(
                side,
                uplo,
                transA,
                diag,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrsm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrsm(
                side,
                uplo,
                transA,
                diag,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrsm(
        side: ::blas::Side,
        uplo: ::blas::Uplo,
        transA: ::blas::Transpose,
        diag: ::blas::Diag,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrsm(
                side,
                uplo,
                transA,
                diag,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(B),
            )
        })
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn ssyrk(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: f32,
        A: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyrk(
                uplo,
                trans,
                alpha,
                ffi::FFI::unwrap_shared(A),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn dsyrk(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: f64,
        A: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyrk(
                uplo,
                trans,
                alpha,
                ffi::FFI::unwrap_shared(A),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn csyrk(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_csyrk(
                uplo,
                trans,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn zsyrk(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zsyrk(
                uplo,
                trans,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// These functions compute a rank-k update of the hermitian matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and C = \alpha A^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cherk(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: f32,
        A: &::types::MatrixComplexF32,
        beta: f32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cherk(
                uplo,
                trans,
                alpha,
                ffi::FFI::unwrap_shared(A),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// These functions compute a rank-k update of the hermitian matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and C = \alpha A^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zherk(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: f64,
        A: &::types::MatrixComplexF64,
        beta: f64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zherk(
                uplo,
                trans,
                alpha,
                ffi::FFI::unwrap_shared(A),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn ssyr2k(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyr2k(
                uplo,
                trans,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn dsyr2k(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyr2k(
                uplo,
                trans,
                alpha,
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn csyr2k(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_csyr2k(
                uplo,
                trans,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn zsyr2k(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zsyr2k(
                uplo,
                trans,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                ::std::mem::transmute(*beta),
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-2k update of the hermitian matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cher2k(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: f32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cher2k(
                uplo,
                trans,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }

    /// This function computes a rank-2k update of the hermitian matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zher2k(
        uplo: ::blas::Uplo,
        trans: ::blas::Transpose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: f64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zher2k(
                uplo,
                trans,
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }
}
