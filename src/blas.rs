//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod level1 {
    use enums;
    use ffi;
    use types::complex::CFFI;

    /// This function computes the sum \alpha + x^T y for the vectors x and y, returning the result
    /// in result.
    ///
    /// Returns `result` if everything went fine.
    pub fn sdsdot(
        alpha: f32,
        x: &::types::VectorF32,
        y: &::types::VectorF32,
    ) -> (enums::Value, f32) {
        let mut result = 0.;
        let ret = unsafe {
            sys::gsl_blas_sdsdot(
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut result,
            )
        };
        (::Value::from(ret), result)
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the
    /// result in result.
    ///
    /// Returns `result` if everything went fine.
    pub fn sdot(x: &::types::VectorF32, y: &::types::VectorF32) -> (enums::Value, f32) {
        let mut result = 0.;
        let ret = unsafe {
            sys::gsl_blas_sdot(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut result,
            )
        };
        (::Value::from(ret), result)
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the
    /// result in result.
    ///
    /// Returns `result` if everything went fine.
    pub fn dsdot(x: &::types::VectorF32, y: &::types::VectorF32) -> (enums::Value, f64) {
        let mut result = 0.;
        let ret = unsafe {
            sys::gsl_blas_dsdot(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut result,
            )
        };
        (::Value::from(ret), result)
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the
    /// result in result.
    ///
    /// Returns `result` if everything went fine.
    pub fn ddot(x: &::types::VectorF64, y: &::types::VectorF64) -> (enums::Value, f64) {
        let mut result = 0.;
        let ret = unsafe {
            sys::gsl_blas_ddot(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut result,
            )
        };
        (::Value::from(ret), result)
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning
    /// the result in dotu.
    ///
    /// Returns `dotu` if everything went fine.
    pub fn cdotu(
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
    ) -> (enums::Value, ::types::ComplexF32) {
        let mut dotu = ::types::ComplexF32::default().unwrap();
        let ret = unsafe {
            sys::gsl_blas_cdotu(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut dotu,
            )
        };
        (::Value::from(ret), ::types::ComplexF32::wrap(dotu))
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning
    /// the result in dotu.
    ///
    /// Returns `dotu` if everything went fine.
    pub fn zdotu(
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
    ) -> (enums::Value, ::types::ComplexF64) {
        let mut dotu = ::types::ComplexF64::default().unwrap();
        let ret = unsafe {
            sys::gsl_blas_zdotu(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut dotu,
            )
        };
        (::Value::from(ret), ::types::ComplexF64::wrap(dotu))
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y,
    /// returning the result in dotc.
    ///
    /// Returns `dotc` if everything went fine.
    pub fn cdotc(
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
    ) -> (enums::Value, ::types::ComplexF32) {
        let mut dotc = ::types::ComplexF32::default().unwrap();
        let ret = unsafe {
            sys::gsl_blas_cdotc(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut dotc,
            )
        };
        (::Value::from(ret), ::types::ComplexF32::wrap(dotc))
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y,
    /// returning the result in dotc.
    ///
    /// Returns `dotc` if everything went fine.
    pub fn zdotc(
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
    ) -> (enums::Value, ::types::ComplexF64) {
        let mut dotc = ::types::ComplexF64::default().unwrap();
        let ret = unsafe {
            sys::gsl_blas_zdotc(
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_shared(y),
                &mut dotc,
            )
        };
        (::Value::from(ret), ::types::ComplexF64::wrap(dotc))
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
    pub fn isamax(x: &::types::VectorF32) -> u64 {
        unsafe { sys::gsl_blas_isamax(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn idamax(x: &::types::VectorF64) -> u64 {
        unsafe { sys::gsl_blas_idamax(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn icamax(x: &::types::VectorComplexF32) -> u64 {
        unsafe { sys::gsl_blas_icamax(ffi::FFI::unwrap_shared(x)) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn izamax(x: &::types::VectorComplexF64) -> u64 {
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
        unsafe { sys::gsl_blas_cscal(::std::mem::transmute(*alpha), ffi::FFI::unwrap_unique(x)) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn zscal(alpha: &::types::ComplexF64, x: &mut ::types::VectorComplexF64) {
        unsafe { sys::gsl_blas_zscal(::std::mem::transmute(*alpha), ffi::FFI::unwrap_unique(x)) }
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
        transA: enums::CblasTranspose,
        alpha: f32,
        A: &::types::MatrixF32,
        x: &::types::VectorF32,
        beta: f32,
        y: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sgemv(
                transA.into(),
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
        transA: enums::CblasTranspose,
        alpha: f64,
        A: &::types::MatrixF64,
        x: &::types::VectorF64,
        beta: f64,
        y: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dgemv(
                transA.into(),
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
        transA: enums::CblasTranspose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        x: &::types::VectorComplexF32,
        beta: &::types::ComplexF32,
        y: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cgemv(
                transA.into(),
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
        transA: enums::CblasTranspose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        x: &::types::VectorComplexF64,
        beta: &::types::ComplexF64,
        y: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zgemv(
                transA.into(),
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
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixF32,
        x: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrmv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixF64,
        x: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrmv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixComplexF32,
        x: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrmv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixComplexF64,
        x: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixF32,
        x: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixF64,
        x: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixComplexF32,
        x: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &::types::MatrixComplexF64,
        x: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_unique(x),
            )
        })
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssymv(
        uplo: enums::CblasUplo,
        alpha: f32,
        A: &::types::MatrixF32,
        x: &::types::VectorF32,
        beta: f32,
        y: &mut ::types::VectorF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssymv(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: f64,
        A: &::types::MatrixF64,
        x: &::types::VectorF64,
        beta: f64,
        y: &mut ::types::VectorF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsymv(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        x: &::types::VectorComplexF32,
        beta: &::types::ComplexF32,
        y: &mut ::types::VectorComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_chemv(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        x: &::types::VectorComplexF64,
        beta: &::types::ComplexF64,
        y: &mut ::types::VectorComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zhemv(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: f32,
        x: &::types::VectorF32,
        A: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyr(
                uplo.into(),
                alpha,
                ffi::FFI::unwrap_shared(x),
                ffi::FFI::unwrap_unique(A),
            )
        })
    }

    /// This function computes the symmetric rank-1 update A = \alpha x x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsyr(
        uplo: enums::CblasUplo,
        alpha: f64,
        x: &::types::VectorF64,
        A: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyr(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: f32,
        x: &::types::VectorComplexF32,
        A: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cher(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: f64,
        x: &::types::VectorComplexF64,
        A: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zher(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: f32,
        x: &::types::VectorF32,
        y: &::types::VectorF32,
        A: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyr2(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: f64,
        x: &::types::VectorF64,
        y: &::types::VectorF64,
        A: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyr2(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF32,
        x: &::types::VectorComplexF32,
        y: &::types::VectorComplexF32,
        A: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cher2(
                uplo.into(),
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
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF64,
        x: &::types::VectorComplexF64,
        y: &::types::VectorComplexF64,
        A: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zher2(
                uplo.into(),
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
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_sgemm(
                transA.into(),
                transB.into(),
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
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dgemm(
                transA.into(),
                transB.into(),
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
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cgemm(
                transA.into(),
                transB.into(),
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
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zgemm(
                transA.into(),
                transB.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssymm(
                side.into(),
                uplo.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsymm(
                side.into(),
                uplo.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_csymm(
                side.into(),
                uplo.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zsymm(
                side.into(),
                uplo.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_chemm(
                side.into(),
                uplo.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zhemm(
                side.into(),
                uplo.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_strsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dtrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ctrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ztrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f32,
        A: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyrk(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f64,
        A: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyrk(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_csyrk(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zsyrk(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f32,
        A: &::types::MatrixComplexF32,
        beta: f32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cherk(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f64,
        A: &::types::MatrixComplexF64,
        beta: f64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zherk(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f32,
        A: &::types::MatrixF32,
        B: &::types::MatrixF32,
        beta: f32,
        C: &mut ::types::MatrixF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_ssyr2k(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f64,
        A: &::types::MatrixF64,
        B: &::types::MatrixF64,
        beta: f64,
        C: &mut ::types::MatrixF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_dsyr2k(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: &::types::ComplexF32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_csyr2k(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: &::types::ComplexF64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zsyr2k(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &::types::ComplexF32,
        A: &::types::MatrixComplexF32,
        B: &::types::MatrixComplexF32,
        beta: f32,
        C: &mut ::types::MatrixComplexF32,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_cher2k(
                uplo.into(),
                trans.into(),
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
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &::types::ComplexF64,
        A: &::types::MatrixComplexF64,
        B: &::types::MatrixComplexF64,
        beta: f64,
        C: &mut ::types::MatrixComplexF64,
    ) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_blas_zher2k(
                uplo.into(),
                trans.into(),
                ::std::mem::transmute(*alpha),
                ffi::FFI::unwrap_shared(A),
                ffi::FFI::unwrap_shared(B),
                beta,
                ffi::FFI::unwrap_unique(C),
            )
        })
    }
}
