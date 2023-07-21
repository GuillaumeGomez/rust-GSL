//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod level1 {
    use crate::{types, Value};
    use ffi::FFI;
    use types::complex::CFFI;

    /// This function computes the sum \alpha + x^T y for the vectors x and y, returning the result
    /// in result.
    ///
    /// Returns `result`.
    #[doc(alias = "gsl_blas_sdsdot")]
    pub fn sdsdot(alpha: f32, x: &types::VectorF32, y: &types::VectorF32) -> Result<f32, Value> {
        let mut result = 0.;
        let ret = unsafe {
            sys::gsl_blas_sdsdot(alpha, x.unwrap_shared(), y.unwrap_shared(), &mut result)
        };
        result_handler!(ret, result)
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the
    /// result in result.
    ///
    /// Returns `result`.
    #[doc(alias = "gsl_blas_sdot")]
    pub fn sdot(x: &types::VectorF32, y: &types::VectorF32) -> Result<f32, Value> {
        let mut result = 0.;
        let ret = unsafe { sys::gsl_blas_sdot(x.unwrap_shared(), y.unwrap_shared(), &mut result) };
        result_handler!(ret, result)
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the
    /// result in result.
    ///
    /// Returns `result`.
    #[doc(alias = "gsl_blas_dsdot")]
    pub fn dsdot(x: &types::VectorF32, y: &types::VectorF32) -> Result<f64, Value> {
        let mut result = 0.;
        let ret = unsafe { sys::gsl_blas_dsdot(x.unwrap_shared(), y.unwrap_shared(), &mut result) };
        result_handler!(ret, result)
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the
    /// result in result.
    ///
    /// Returns `result`.
    #[doc(alias = "gsl_blas_ddot")]
    pub fn ddot(x: &types::VectorF64, y: &types::VectorF64) -> Result<f64, Value> {
        let mut result = 0.;
        let ret = unsafe { sys::gsl_blas_ddot(x.unwrap_shared(), y.unwrap_shared(), &mut result) };
        result_handler!(ret, result)
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning
    /// the result in dotu.
    ///
    /// Returns `dotu`.
    #[doc(alias = "gsl_blas_cdotu")]
    pub fn cdotu(
        x: &types::VectorComplexF32,
        y: &types::VectorComplexF32,
    ) -> Result<types::ComplexF32, Value> {
        let mut dotu = types::ComplexF32::default().unwrap();
        let ret = unsafe { sys::gsl_blas_cdotu(x.unwrap_shared(), y.unwrap_shared(), &mut dotu) };
        result_handler!(ret, types::ComplexF32::wrap(dotu))
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning
    /// the result in dotu.
    ///
    /// Returns `dotu`.
    #[doc(alias = "gsl_blas_zdotu")]
    pub fn zdotu(
        x: &types::VectorComplexF64,
        y: &types::VectorComplexF64,
    ) -> Result<types::ComplexF64, Value> {
        let mut dotu = types::ComplexF64::default().unwrap();
        let ret = unsafe { sys::gsl_blas_zdotu(x.unwrap_shared(), y.unwrap_shared(), &mut dotu) };
        result_handler!(ret, types::ComplexF64::wrap(dotu))
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y,
    /// returning the result in dotc.
    ///
    /// Returns `dotc`.
    #[doc(alias = "gsl_blas_cdotc")]
    pub fn cdotc(
        x: &types::VectorComplexF32,
        y: &types::VectorComplexF32,
    ) -> Result<types::ComplexF32, Value> {
        let mut dotc = types::ComplexF32::default().unwrap();
        let ret = unsafe { sys::gsl_blas_cdotc(x.unwrap_shared(), y.unwrap_shared(), &mut dotc) };
        result_handler!(ret, types::ComplexF32::wrap(dotc))
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y,
    /// returning the result in dotc.
    ///
    /// Returns `dotc`.
    #[doc(alias = "gsl_blas_zdotc")]
    pub fn zdotc(
        x: &types::VectorComplexF64,
        y: &types::VectorComplexF64,
    ) -> Result<types::ComplexF64, Value> {
        let mut dotc = types::ComplexF64::default().unwrap();
        let ret = unsafe { sys::gsl_blas_zdotc(x.unwrap_shared(), y.unwrap_shared(), &mut dotc) };
        result_handler!(ret, types::ComplexF64::wrap(dotc))
    }

    /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
    #[doc(alias = "gsl_blas_snrm2")]
    pub fn snrm2(x: &types::VectorF32) -> f32 {
        unsafe { sys::gsl_blas_snrm2(x.unwrap_shared()) }
    }

    /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
    #[doc(alias = "gsl_blas_dnrm2")]
    pub fn dnrm2(x: &types::VectorF64) -> f64 {
        unsafe { sys::gsl_blas_dnrm2(x.unwrap_shared()) }
    }

    /// This function computes the Euclidean norm of the complex vector x,
    ///
    /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
    #[doc(alias = "gsl_blas_scnrm2")]
    pub fn scnrm2(x: &types::VectorComplexF32) -> f32 {
        unsafe { sys::gsl_blas_scnrm2(x.unwrap_shared()) }
    }

    /// This function computes the Euclidean norm of the complex vector x,
    ///
    /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
    #[doc(alias = "gsl_blas_dznrm2")]
    pub fn dznrm2(x: &types::VectorComplexF64) -> f64 {
        unsafe { sys::gsl_blas_dznrm2(x.unwrap_shared()) }
    }

    /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
    #[doc(alias = "gsl_blas_sasum")]
    pub fn sasum(x: &types::VectorF32) -> f32 {
        unsafe { sys::gsl_blas_sasum(x.unwrap_shared()) }
    }

    /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
    #[doc(alias = "gsl_blas_dasum")]
    pub fn dasum(x: &types::VectorF64) -> f64 {
        unsafe { sys::gsl_blas_dasum(x.unwrap_shared()) }
    }

    /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
    #[doc(alias = "gsl_blas_scasum")]
    pub fn scasum(x: &types::VectorComplexF32) -> f32 {
        unsafe { sys::gsl_blas_scasum(x.unwrap_shared()) }
    }

    /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
    #[doc(alias = "gsl_blas_dzasum")]
    pub fn dzasum(x: &types::VectorComplexF64) -> f64 {
        unsafe { sys::gsl_blas_dzasum(x.unwrap_shared()) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    #[doc(alias = "gsl_blas_isamax")]
    pub fn isamax(x: &types::VectorF32) -> usize {
        unsafe { sys::gsl_blas_isamax(x.unwrap_shared()) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    #[doc(alias = "gsl_blas_idamax")]
    pub fn idamax(x: &types::VectorF64) -> usize {
        unsafe { sys::gsl_blas_idamax(x.unwrap_shared()) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    #[doc(alias = "gsl_blas_icamax")]
    pub fn icamax(x: &types::VectorComplexF32) -> usize {
        unsafe { sys::gsl_blas_icamax(x.unwrap_shared()) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    #[doc(alias = "gsl_blas_izamax")]
    pub fn izamax(x: &types::VectorComplexF64) -> usize {
        unsafe { sys::gsl_blas_izamax(x.unwrap_shared()) }
    }

    /// This function exchanges the elements of the vectors x and y.
    #[doc(alias = "gsl_blas_sswap")]
    pub fn sswap(x: &mut types::VectorF32, y: &mut types::VectorF32) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_sswap(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function exchanges the elements of the vectors x and y.
    #[doc(alias = "gsl_blas_dswap")]
    pub fn dswap(x: &mut types::VectorF64, y: &mut types::VectorF64) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_dswap(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function exchanges the elements of the vectors x and y.
    #[doc(alias = "gsl_blas_cswap")]
    pub fn cswap(
        x: &mut types::VectorComplexF32,
        y: &mut types::VectorComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_cswap(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function exchanges the elements of the vectors x and y.
    #[doc(alias = "gsl_blas_zswap")]
    pub fn zswap(
        x: &mut types::VectorComplexF64,
        y: &mut types::VectorComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_zswap(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function copy the elements of the vector x into the vector y.
    #[doc(alias = "gsl_blas_scopy")]
    pub fn scopy(x: &mut types::VectorF32, y: &mut types::VectorF32) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_scopy(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function copy the elements of the vector x into the vector y.
    #[doc(alias = "gsl_blas_dcopy")]
    pub fn dcopy(x: &mut types::VectorF64, y: &mut types::VectorF64) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_dcopy(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function copy the elements of the vector x into the vector y.
    #[doc(alias = "gsl_blas_ccopy")]
    pub fn ccopy(
        x: &mut types::VectorComplexF32,
        y: &mut types::VectorComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_ccopy(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function copy the elements of the vector x into the vector y.
    #[doc(alias = "gsl_blas_zcopy")]
    pub fn zcopy(
        x: &mut types::VectorComplexF64,
        y: &mut types::VectorComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_zcopy(x.unwrap_unique(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    #[doc(alias = "gsl_blas_saxpy")]
    pub fn saxpy(alpha: f32, x: &types::VectorF32, y: &mut types::VectorF32) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_saxpy(alpha, x.unwrap_shared(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    #[doc(alias = "gsl_blas_daxpy")]
    pub fn daxpy(alpha: f64, x: &types::VectorF64, y: &mut types::VectorF64) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_daxpy(alpha, x.unwrap_shared(), y.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    #[doc(alias = "gsl_blas_caxpy")]
    pub fn caxpy(
        alpha: &types::ComplexF32,
        x: &types::VectorComplexF32,
        y: &mut types::VectorComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_caxpy(
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    #[doc(alias = "gsl_blas_zaxpy")]
    pub fn zaxpy(
        alpha: &types::ComplexF64,
        x: &types::VectorComplexF64,
        y: &mut types::VectorComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zaxpy(
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    #[doc(alias = "gsl_blas_sscal")]
    pub fn sscal(alpha: f32, x: &mut types::VectorF32) {
        unsafe { sys::gsl_blas_sscal(alpha, x.unwrap_unique()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    #[doc(alias = "gsl_blas_dscal")]
    pub fn dscal(alpha: f64, x: &mut types::VectorF64) {
        unsafe { sys::gsl_blas_dscal(alpha, x.unwrap_unique()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    #[doc(alias = "gsl_blas_cscal")]
    pub fn cscal(alpha: &types::ComplexF32, x: &mut types::VectorComplexF32) {
        unsafe { sys::gsl_blas_cscal(::std::mem::transmute(*alpha), x.unwrap_unique()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    #[doc(alias = "gsl_blas_zscal")]
    pub fn zscal(alpha: &types::ComplexF64, x: &mut types::VectorComplexF64) {
        unsafe { sys::gsl_blas_zscal(::std::mem::transmute(*alpha), x.unwrap_unique()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    #[doc(alias = "gsl_blas_csscal")]
    pub fn csscal(alpha: f32, x: &mut types::VectorComplexF32) {
        unsafe { sys::gsl_blas_csscal(alpha, x.unwrap_unique()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    #[doc(alias = "gsl_blas_zdscal")]
    pub fn zdscal(alpha: f64, x: &mut types::VectorComplexF64) {
        unsafe { sys::gsl_blas_zdscal(alpha, x.unwrap_unique()) }
    }

    /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
    ///
    /// ```text
    /// [  c  s ] [ a ] = [ r ]
    ///
    /// [ -s  c ] [ b ]   [ 0 ]
    /// ```
    ///
    /// The variables a and b are overwritten by the routine.
    #[doc(alias = "gsl_blas_srotg")]
    pub fn srotg(a: &mut [f32], b: &mut [f32], c: &mut [f32], d: &mut [f32]) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_srotg(
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                c.as_mut_ptr(),
                d.as_mut_ptr(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
    ///
    /// ```text
    /// [  c  s ] [ a ] = [ r ]
    ///
    /// [ -s  c ] [ b ]   [ 0 ]
    /// ```
    ///
    /// The variables a and b are overwritten by the routine.
    #[doc(alias = "gsl_blas_drotg")]
    pub fn drotg(a: &mut [f64], b: &mut [f64], c: &mut [f64], d: &mut [f64]) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_drotg(
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                c.as_mut_ptr(),
                d.as_mut_ptr(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
    #[doc(alias = "gsl_blas_srot")]
    pub fn srot(
        a: &mut types::VectorF32,
        b: &mut types::VectorF32,
        c: f32,
        d: f32,
    ) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_srot(a.unwrap_unique(), b.unwrap_unique(), c, d) };
        result_handler!(ret, ())
    }

    /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
    #[doc(alias = "gsl_blas_drot")]
    pub fn drot(
        a: &mut types::VectorF64,
        b: &mut types::VectorF64,
        c: f64,
        d: f64,
    ) -> Result<(), Value> {
        let ret = unsafe { sys::gsl_blas_drot(a.unwrap_unique(), b.unwrap_unique(), c, d) };
        result_handler!(ret, ())
    }

    /// This function computes a modified Givens transformation.
    /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
    #[doc(alias = "gsl_blas_srotmg")]
    pub fn srotmg(
        d1: &mut [f32],
        d2: &mut [f32],
        b1: &mut [f32],
        b2: f32,
        P: &mut [f32],
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_srotmg(
                d1.as_mut_ptr(),
                d2.as_mut_ptr(),
                b1.as_mut_ptr(),
                b2,
                P.as_mut_ptr(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a modified Givens transformation.
    /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
    #[doc(alias = "gsl_blas_drotmg")]
    pub fn drotmg(
        d1: &mut [f64],
        d2: &mut [f64],
        b1: &mut [f64],
        b2: f64,
        P: &mut [f64],
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_drotmg(
                d1.as_mut_ptr(),
                d2.as_mut_ptr(),
                b1.as_mut_ptr(),
                b2,
                P.as_mut_ptr(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function applies a modified Givens transformation.
    #[doc(alias = "gsl_blas_srotm")]
    pub fn srotm(
        x: &mut types::VectorF32,
        y: &mut types::VectorF32,
        P: &mut [f32],
    ) -> Result<(), Value> {
        let ret =
            unsafe { sys::gsl_blas_srotm(x.unwrap_unique(), y.unwrap_unique(), P.as_mut_ptr()) };
        result_handler!(ret, ())
    }

    /// This function applies a modified Givens transformation.
    #[doc(alias = "gsl_blas_drotm")]
    pub fn drotm(
        x: &mut types::VectorF64,
        y: &mut types::VectorF64,
        P: &mut [f64],
    ) -> Result<(), Value> {
        let ret =
            unsafe { sys::gsl_blas_drotm(x.unwrap_unique(), y.unwrap_unique(), P.as_mut_ptr()) };
        result_handler!(ret, ())
    }
}

pub mod level2 {
    use crate::{enums, types, Value};
    use ffi::FFI;

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_sgemv")]
    pub fn sgemv(
        transA: enums::CblasTranspose,
        alpha: f32,
        A: &types::MatrixF32,
        x: &types::VectorF32,
        beta: f32,
        y: &mut types::VectorF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_sgemv(
                transA.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_dgemv")]
    pub fn dgemv(
        transA: enums::CblasTranspose,
        alpha: f64,
        A: &types::MatrixF64,
        x: &types::VectorF64,
        beta: f64,
        y: &mut types::VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dgemv(
                transA.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_cgemv")]
    pub fn cgemv(
        transA: enums::CblasTranspose,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        x: &types::VectorComplexF32,
        beta: &types::ComplexF32,
        y: &mut types::VectorComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_cgemv(
                transA.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                x.unwrap_shared(),
                ::std::mem::transmute(*beta),
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    #[doc(alias = "gsl_blas_zgemv")]
    pub fn zgemv(
        transA: enums::CblasTranspose,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        x: &types::VectorComplexF64,
        beta: &types::ComplexF64,
        y: &mut types::VectorComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zgemv(
                transA.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                x.unwrap_shared(),
                ::std::mem::transmute(*beta),
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_strmv")]
    pub fn strmv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixF32,
        x: &mut types::VectorF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_strmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_dtrmv")]
    pub fn dtrmv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixF64,
        x: &mut types::VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dtrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ctrmv")]
    pub fn ctrmv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixComplexF32,
        x: &mut types::VectorComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ctrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ztrmv")]
    pub fn ztrmv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixComplexF64,
        x: &mut types::VectorComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ztrmv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_strsv")]
    pub fn strsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixF32,
        x: &mut types::VectorF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_strsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_dtrsv")]
    pub fn dtrsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixF64,
        x: &mut types::VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dtrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ctrsv")]
    pub fn ctrsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixComplexF32,
        x: &mut types::VectorComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ctrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ztrsv")]
    pub fn ztrsv(
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        A: &types::MatrixComplexF64,
        x: &mut types::VectorComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ztrsv(
                uplo.into(),
                transA.into(),
                diag.into(),
                A.unwrap_shared(),
                x.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssymv")]
    pub fn ssymv(
        uplo: enums::CblasUplo,
        alpha: f32,
        A: &types::MatrixF32,
        x: &types::VectorF32,
        beta: f32,
        y: &mut types::VectorF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ssymv(
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsymv")]
    pub fn dsymv(
        uplo: enums::CblasUplo,
        alpha: f64,
        A: &types::MatrixF64,
        x: &types::VectorF64,
        beta: f64,
        y: &mut types::VectorF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dsymv(
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                x.unwrap_shared(),
                beta,
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored. When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically assumed to be zero and are not referenced.
    #[doc(alias = "gsl_blas_chemv")]
    pub fn chemv(
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        x: &types::VectorComplexF32,
        beta: &types::ComplexF32,
        y: &mut types::VectorComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_chemv(
                uplo.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                x.unwrap_shared(),
                ::std::mem::transmute(*beta),
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored. When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically assumed to be zero and are not referenced.
    #[doc(alias = "gsl_blas_zhemv")]
    pub fn zhemv(
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        x: &types::VectorComplexF64,
        beta: &types::ComplexF64,
        y: &mut types::VectorComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zhemv(
                uplo.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                x.unwrap_shared(),
                ::std::mem::transmute(*beta),
                y.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    #[doc(alias = "gsl_blas_sger")]
    pub fn sger(
        alpha: f32,
        x: &types::VectorF32,
        y: &types::VectorF32,
        A: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_sger(
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    #[doc(alias = "gsl_blas_dger")]
    pub fn dger(
        alpha: f64,
        x: &types::VectorF64,
        y: &types::VectorF64,
        A: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dger(
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    #[doc(alias = "gsl_blas_cgeru")]
    pub fn cgeru(
        alpha: &types::ComplexF32,
        x: &types::VectorComplexF32,
        y: &types::VectorComplexF32,
        A: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_cgeru(
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    #[doc(alias = "gsl_blas_zgeru")]
    pub fn zgeru(
        alpha: &types::ComplexF64,
        x: &types::VectorComplexF64,
        y: &types::VectorComplexF64,
        A: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zgeru(
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the conjugate rank-1 update A = \alpha x y^H + A of the matrix A.
    #[doc(alias = "gsl_blas_cgerc")]
    pub fn cgerc(
        alpha: &types::ComplexF32,
        x: &types::VectorComplexF32,
        y: &types::VectorComplexF32,
        A: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_cgerc(
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the conjugate rank-1 update A = \alpha x y^H + A of the matrix A.
    #[doc(alias = "gsl_blas_zgerc")]
    pub fn zgerc(
        alpha: &types::ComplexF64,
        x: &types::VectorComplexF64,
        y: &types::VectorComplexF64,
        A: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zgerc(
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the symmetric rank-1 update A = \alpha x x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssyr")]
    pub fn ssyr(
        uplo: enums::CblasUplo,
        alpha: f32,
        x: &types::VectorF32,
        A: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret =
            unsafe { sys::gsl_blas_ssyr(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function computes the symmetric rank-1 update A = \alpha x x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsyr")]
    pub fn dsyr(
        uplo: enums::CblasUplo,
        alpha: f64,
        x: &types::VectorF64,
        A: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret =
            unsafe { sys::gsl_blas_dsyr(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// These functions compute the hermitian rank-1 update A = \alpha x x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_cher")]
    pub fn cher(
        uplo: enums::CblasUplo,
        alpha: f32,
        x: &types::VectorComplexF32,
        A: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret =
            unsafe { sys::gsl_blas_cher(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// These functions compute the hermitian rank-1 update A = \alpha x x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zher")]
    pub fn zher(
        uplo: enums::CblasUplo,
        alpha: f64,
        x: &types::VectorComplexF64,
        A: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret =
            unsafe { sys::gsl_blas_zher(uplo.into(), alpha, x.unwrap_shared(), A.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// These functions compute the symmetric rank-2 update A = \alpha x y^T + \alpha y x^T + A of the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssyr2")]
    pub fn ssyr2(
        uplo: enums::CblasUplo,
        alpha: f32,
        x: &types::VectorF32,
        y: &types::VectorF32,
        A: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ssyr2(
                uplo.into(),
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute the symmetric rank-2 update A = \alpha x y^T + \alpha y x^T + A of the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsyr2")]
    pub fn dsyr2(
        uplo: enums::CblasUplo,
        alpha: f64,
        x: &types::VectorF64,
        y: &types::VectorF64,
        A: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dsyr2(
                uplo.into(),
                alpha,
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute the hermitian rank-2 update A = \alpha x y^H + \alpha^* y x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_cher2")]
    pub fn cher2(
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF32,
        x: &types::VectorComplexF32,
        y: &types::VectorComplexF32,
        A: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_cher2(
                uplo.into(),
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute the hermitian rank-2 update A = \alpha x y^H + \alpha^* y x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zher2")]
    pub fn zher2(
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF64,
        x: &types::VectorComplexF64,
        y: &types::VectorComplexF64,
        A: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zher2(
                uplo.into(),
                ::std::mem::transmute(*alpha),
                x.unwrap_shared(),
                y.unwrap_shared(),
                A.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }
}

pub mod level3 {
    use crate::{enums, types, Value};
    use ffi::FFI;

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_sgemm")]
    pub fn sgemm(
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: f32,
        A: &types::MatrixF32,
        B: &types::MatrixF32,
        beta: f32,
        C: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_sgemm(
                transA.into(),
                transB.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_dgemm")]
    pub fn dgemm(
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: f64,
        A: &types::MatrixF64,
        B: &types::MatrixF64,
        beta: f64,
        C: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dgemm(
                transA.into(),
                transB.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_cgemm")]
    pub fn cgemm(
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        B: &types::MatrixComplexF32,
        beta: &types::ComplexF32,
        C: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_cgemm(
                transA.into(),
                transB.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    #[doc(alias = "gsl_blas_zgemm")]
    pub fn zgemm(
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        B: &types::MatrixComplexF64,
        beta: &types::ComplexF64,
        C: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zgemm(
                transA.into(),
                transB.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_ssymm")]
    pub fn ssymm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: f32,
        A: &types::MatrixF32,
        B: &types::MatrixF32,
        beta: f32,
        C: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ssymm(
                side.into(),
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_dsymm")]
    pub fn dsymm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: f64,
        A: &types::MatrixF64,
        B: &types::MatrixF64,
        beta: f64,
        C: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dsymm(
                side.into(),
                uplo.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_csymm")]
    pub fn csymm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        B: &types::MatrixComplexF32,
        beta: &types::ComplexF32,
        C: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_csymm(
                side.into(),
                uplo.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    #[doc(alias = "gsl_blas_zsymm")]
    pub fn zsymm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        B: &types::MatrixComplexF64,
        beta: &types::ComplexF64,
        C: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zsymm(
                side.into(),
                uplo.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is Left and C = \alpha B A + \beta C for Side is Right, where the matrix A is hermitian.
    /// When Uplo is Upper then the upper triangle and diagonal of A are used, and when Uplo is Lower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_chemm")]
    pub fn chemm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        B: &types::MatrixComplexF32,
        beta: &types::ComplexF32,
        C: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_chemm(
                side.into(),
                uplo.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is hermitian.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zhemm")]
    pub fn zhemm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        B: &types::MatrixComplexF64,
        beta: &types::ComplexF64,
        C: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zhemm(
                side.into(),
                uplo.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_strmm")]
    pub fn strmm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f32,
        A: &types::MatrixF32,
        B: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_strmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_dtrmm")]
    pub fn dtrmm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f64,
        A: &types::MatrixF64,
        B: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dtrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ctrmm")]
    pub fn ctrmm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        B: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ctrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ztrmm")]
    pub fn ztrmm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        B: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ztrmm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_strsm")]
    pub fn strsm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f32,
        A: &types::MatrixF32,
        B: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_strsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_dtrsm")]
    pub fn dtrsm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: f64,
        A: &types::MatrixF64,
        B: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dtrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ctrsm")]
    pub fn ctrsm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        B: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ctrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    #[doc(alias = "gsl_blas_ztrsm")]
    pub fn ztrsm(
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        B: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ztrsm(
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_ssyrk")]
    pub fn ssyrk(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f32,
        A: &types::MatrixF32,
        beta: f32,
        C: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ssyrk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_dsyrk")]
    pub fn dsyrk(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f64,
        A: &types::MatrixF64,
        beta: f64,
        C: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dsyrk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_csyrk")]
    pub fn csyrk(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        beta: &types::ComplexF32,
        C: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_csyrk(
                uplo.into(),
                trans.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_zsyrk")]
    pub fn zsyrk(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        beta: &types::ComplexF64,
        C: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zsyrk(
                uplo.into(),
                trans.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute a rank-k update of the hermitian matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and C = \alpha A^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_cherk")]
    pub fn cherk(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f32,
        A: &types::MatrixComplexF32,
        beta: f32,
        C: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_cherk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// These functions compute a rank-k update of the hermitian matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and C = \alpha A^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zherk")]
    pub fn zherk(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f64,
        A: &types::MatrixComplexF64,
        beta: f64,
        C: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zherk(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_ssyr2k")]
    pub fn ssyr2k(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f32,
        A: &types::MatrixF32,
        B: &types::MatrixF32,
        beta: f32,
        C: &mut types::MatrixF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_ssyr2k(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_dsyr2k")]
    pub fn dsyr2k(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: f64,
        A: &types::MatrixF64,
        B: &types::MatrixF64,
        beta: f64,
        C: &mut types::MatrixF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_dsyr2k(
                uplo.into(),
                trans.into(),
                alpha,
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_csyr2k")]
    pub fn csyr2k(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        B: &types::MatrixComplexF32,
        beta: &types::ComplexF32,
        C: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_csyr2k(
                uplo.into(),
                trans.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    #[doc(alias = "gsl_blas_zsyr2k")]
    pub fn zsyr2k(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        B: &types::MatrixComplexF64,
        beta: &types::ComplexF64,
        C: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zsyr2k(
                uplo.into(),
                trans.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                ::std::mem::transmute(*beta),
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-2k update of the hermitian matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_cher2k")]
    pub fn cher2k(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &types::ComplexF32,
        A: &types::MatrixComplexF32,
        B: &types::MatrixComplexF32,
        beta: f32,
        C: &mut types::MatrixComplexF32,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_cher2k(
                uplo.into(),
                trans.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    /// This function computes a rank-2k update of the hermitian matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    #[doc(alias = "gsl_blas_zher2k")]
    pub fn zher2k(
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        alpha: &types::ComplexF64,
        A: &types::MatrixComplexF64,
        B: &types::MatrixComplexF64,
        beta: f64,
        C: &mut types::MatrixComplexF64,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::gsl_blas_zher2k(
                uplo.into(),
                trans.into(),
                ::std::mem::transmute(*alpha),
                A.unwrap_shared(),
                B.unwrap_shared(),
                beta,
                C.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }
}
