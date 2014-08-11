//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod level1 {
    use gsl;

    /// This function computes the sum \alpha + x^T y for the vectors x and y, returning the result in result.
    pub fn sdsdot(alpha: f32, x: &gsl::VectorFloat, y: &gsl::VectorFloat, result: &mut f32) -> i32 {
        unsafe { ::ffi::gsl_blas_sdsdot(alpha, x.get_ffi() as *const ::ffi::gsl_vector_float,
            y.get_ffi() as *const ::ffi::gsl_vector_float, result) }
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
    pub fn sdot(x: &gsl::VectorFloat, y: &gsl::VectorFloat, result: &mut f32) -> i32 {
        unsafe { ::ffi::gsl_blas_sdot(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi() as *const ::ffi::gsl_vector_float,
            result) }
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
    pub fn dsdot(x: &gsl::VectorFloat, y: &gsl::VectorFloat, result: &mut f64) -> i32 {
        unsafe { ::ffi::gsl_blas_dsdot(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi() as *const ::ffi::gsl_vector_float,
            result) }
    }

    /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
    pub fn ddot(x: &gsl::Vector, y: &gsl::Vector, result: &mut f64) -> i32 {
        unsafe { ::ffi::gsl_blas_ddot(x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi() as *const ::ffi::gsl_vector, result) }
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning the result in dotu.
    pub fn cdotu(x: &gsl::VectorComplexFloat, y: &gsl::VectorComplexFloat, dotu: &mut gsl::ComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cdotu(x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
            y.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(dotu)) }
    }

    /// This function computes the complex scalar product x^T y for the vectors x and y, returning the result in dotu.
    pub fn zdotu(x: &gsl::VectorComplex, y: &gsl::VectorComplex, dotu: &mut gsl::Complex) -> i32 {
        unsafe { ::ffi::gsl_blas_zdotu(x.get_ffi() as *const ::ffi::gsl_vector_complex,
            y.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(dotu)) }
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y, returning the result in dotc.
    pub fn cdotc(x: &gsl::VectorComplexFloat, y: &gsl::VectorComplexFloat, dotc: &mut gsl::ComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cdotc(x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
            y.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(dotc)) }
    }

    /// This function computes the complex conjugate scalar product x^H y for the vectors x and y, returning the result in dotc.
    pub fn zdotc(x: &gsl::VectorComplex, y: &gsl::VectorComplex, dotc: &mut gsl::Complex) -> i32 {
        unsafe { ::ffi::gsl_blas_zdotc(x.get_ffi() as *const ::ffi::gsl_vector_complex,
            y.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(dotc)) }
    }

    /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
    pub fn snrm2(x: &gsl::VectorFloat) -> f32 {
        unsafe { ::ffi::gsl_blas_snrm2(x.get_ffi() as *const ::ffi::gsl_vector_float) }
    }

    /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
    pub fn dnrm2(x: &gsl::Vector) -> f64 {
        unsafe { ::ffi::gsl_blas_dnrm2(x.get_ffi() as *const ::ffi::gsl_vector) }
    }

    /// This function computes the Euclidean norm of the complex vector x,
    /// 
    /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
    pub fn scnrm2(x: &gsl::VectorComplexFloat) -> f32 {
        unsafe { ::ffi::gsl_blas_scnrm2(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
    }

    /// This function computes the Euclidean norm of the complex vector x,
    /// 
    /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
    pub fn dznrm2(x: &gsl::VectorComplex) -> f64 {
        unsafe { ::ffi::gsl_blas_dznrm2(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
    }

    /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
    pub fn sasum(x: &gsl::VectorFloat) -> f32 {
        unsafe { ::ffi::gsl_blas_sasum(x.get_ffi() as *const ::ffi::gsl_vector_float) }
    }

    /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
    pub fn dasum(x: &gsl::Vector) -> f64 {
        unsafe { ::ffi::gsl_blas_dasum(x.get_ffi() as *const ::ffi::gsl_vector) }
    }

    /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
    pub fn scasum(x: &gsl::VectorComplexFloat) -> f32 {
        unsafe { ::ffi::gsl_blas_scasum(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
    }

    /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
    pub fn dzasum(x: &gsl::VectorComplex) -> f64 {
        unsafe { ::ffi::gsl_blas_dzasum(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn isamax(x: &gsl::VectorFloat) -> u32 {
        unsafe { ::ffi::gsl_blas_isamax(x.get_ffi() as *const ::ffi::gsl_vector_float) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn idamax(x: &gsl::Vector) -> u32 {
        unsafe { ::ffi::gsl_blas_idamax(x.get_ffi() as *const ::ffi::gsl_vector) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn icamax(x: &gsl::VectorComplexFloat) -> u32 {
        unsafe { ::ffi::gsl_blas_icamax(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
    }

    /// This function returns the index of the largest element of the vector x.
    /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
    /// If the largest value occurs several times then the index of the first occurrence is returned.
    pub fn izamax(x: &gsl::VectorComplex) -> u32 {
        unsafe { ::ffi::gsl_blas_izamax(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn sswap(x: &mut gsl::VectorFloat, y: &mut gsl::VectorFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_sswap(x.get_ffi(), y.get_ffi()) }
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn dswap(x: &mut gsl::Vector, y: &mut gsl::Vector) -> i32 {
        unsafe { ::ffi::gsl_blas_dswap(x.get_ffi(), y.get_ffi()) }
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn cswap(x: &mut gsl::VectorComplexFloat, y: &mut gsl::VectorComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cswap(x.get_ffi(), y.get_ffi()) }
    }

    /// This function exchanges the elements of the vectors x and y.
    pub fn zswap(x: &mut gsl::VectorComplex, y: &mut gsl::VectorComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zswap(x.get_ffi(), y.get_ffi()) }
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn scopy(x: &mut gsl::VectorFloat, y: &mut gsl::VectorFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_scopy(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi()) }
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn dcopy(x: &mut gsl::Vector, y: &mut gsl::Vector) -> i32 {
        unsafe { ::ffi::gsl_blas_dcopy(x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi()) }
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn ccopy(x: &mut gsl::VectorComplexFloat, y: &mut gsl::VectorComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ccopy(x.get_ffi() as *const ::ffi::gsl_vector_complex_float, y.get_ffi()) }
    }

    /// This function copy the elements of the vector x into the vector y.
    pub fn zcopy(x: &mut gsl::VectorComplex, y: &mut gsl::VectorComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zcopy(x.get_ffi() as *const ::ffi::gsl_vector_complex, y.get_ffi()) }
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn saxpy(alpha: f32, x: &gsl::VectorFloat, y: &mut gsl::VectorFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_saxpy(alpha, x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi()) }
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn daxpy(alpha: f64, x: &gsl::Vector, y: &mut gsl::Vector) -> i32 {
        unsafe { ::ffi::gsl_blas_daxpy(alpha, x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi()) }
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn caxpy(alpha: &gsl::ComplexFloat, x: &gsl::VectorComplexFloat, y: &mut gsl::VectorComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_caxpy(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex_float, y.get_ffi()) }
    }

    /// This function computes the sum y = \alpha x + y for the vectors x and y.
    pub fn zaxpy(alpha: &gsl::Complex, x: &gsl::VectorComplex, y: &mut gsl::VectorComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zaxpy(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex, y.get_ffi()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn sscal(alpha: f32, x: &mut gsl::VectorFloat) {
        unsafe { ::ffi::gsl_blas_sscal(alpha, x.get_ffi()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn dscal(alpha: f64, x: &mut gsl::Vector) {
        unsafe { ::ffi::gsl_blas_dscal(alpha, x.get_ffi()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn cscal(alpha: &gsl::ComplexFloat, x: &mut gsl::VectorComplexFloat) {
        unsafe { ::ffi::gsl_blas_cscal(::std::mem::transmute(*alpha), x.get_ffi()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn zscal(alpha: &gsl::Complex, x: &mut gsl::VectorComplex) {
        unsafe { ::ffi::gsl_blas_zscal(::std::mem::transmute(*alpha), x.get_ffi()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn csscal(alpha: f32, x: &mut gsl::VectorComplexFloat) {
        unsafe { ::ffi::gsl_blas_csscal(alpha, x.get_ffi()) }
    }

    /// This function rescales the vector x by the multiplicative factor alpha.
    pub fn zdscal(alpha: f64, x: &mut gsl::VectorComplex) {
        unsafe { ::ffi::gsl_blas_zdscal(alpha, x.get_ffi()) }
    }

    /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
    /// 
    /// [  c  s ] [ a ] = [ r ]
    /// 
    /// [ -s  c ] [ b ]   [ 0 ]
    /// 
    /// The variables a and b are overwritten by the routine.
    pub fn srotg(a: &mut [f32], b: &mut [f32], c: &mut [f32], d: &mut [f32]) -> i32 {
        unsafe { ::ffi::gsl_blas_srotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), d.as_mut_ptr()) }
    }

    /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
    /// 
    /// [  c  s ] [ a ] = [ r ]
    /// 
    /// [ -s  c ] [ b ]   [ 0 ]
    /// 
    /// The variables a and b are overwritten by the routine.
    pub fn drotg(a: &mut [f64], b: &mut [f64], c: &mut [f64], d: &mut [f64]) -> i32 {
        unsafe { ::ffi::gsl_blas_drotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), d.as_mut_ptr()) }
    }

    /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
    pub fn srot(a: &mut gsl::VectorFloat, b: &mut gsl::VectorFloat, c: f32, d: f32) -> i32 {
        unsafe { ::ffi::gsl_blas_srot(a.get_ffi(), b.get_ffi(), c, d) }
    }

    /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
    pub fn drot(a: &mut gsl::Vector, b: &mut gsl::Vector, c: f64, d: f64) -> i32 {
        unsafe { ::ffi::gsl_blas_drot(a.get_ffi(), b.get_ffi(), c, d) }
    }

    /// This function computes a modified Givens transformation.
    /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
    pub fn srotmg(d1: &mut [f32], d2: &mut [f32], b1: &mut [f32], b2: f32, P: &mut [f32]) -> i32 {
        unsafe { ::ffi::gsl_blas_srotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2, P.as_mut_ptr()) }
    }

    /// This function computes a modified Givens transformation.
    /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
    pub fn drotmg(d1: &mut [f64], d2: &mut [f64], b1: &mut [f64], b2: f64, P: &mut [f64]) -> i32 {
        unsafe { ::ffi::gsl_blas_drotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2, P.as_mut_ptr()) }
    }

    /// This function applies a modified Givens transformation.
    pub fn srotm(x: &mut gsl::VectorFloat, y: &mut gsl::VectorFloat, P: &mut [f32]) -> i32 {
        unsafe { ::ffi::gsl_blas_srotm(x. get_ffi(), y.get_ffi(), P.as_mut_ptr()) }
    }

    /// This function applies a modified Givens transformation.
    pub fn drotm(x: &mut gsl::Vector, y: &mut gsl::Vector, P: &mut [f64]) -> i32 {
        unsafe { ::ffi::gsl_blas_drotm(x. get_ffi(), y.get_ffi(), P.as_mut_ptr()) }
    }
}

pub mod level2 {
    use gsl;

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn sgemv(transA: ::enums::CblasTranspose, alpha: f32, A: &gsl::MatrixFloat, x: &gsl::VectorFloat, beta: f32, y: &mut gsl::VectorFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_sgemv(transA, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float, x.get_ffi() as *const ::ffi::gsl_vector_float,
            beta, y.get_ffi()) }
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn dgemv(transA: ::enums::CblasTranspose, alpha: f64, A: &gsl::Matrix, x: &gsl::Vector, beta: f64, y: &mut gsl::Vector) -> i32 {
        unsafe { ::ffi::gsl_blas_dgemv(transA, alpha, A.get_ffi() as *const ::ffi::gsl_matrix, x.get_ffi() as *const ::ffi::gsl_vector,
            beta, y.get_ffi()) }
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn cgemv(transA: ::enums::CblasTranspose, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat, x: &gsl::VectorComplexFloat,
        beta: &gsl::ComplexFloat, y: &mut gsl::VectorComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cgemv(transA, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            x.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(*beta), y.get_ffi()) }
    }

    /// This function computes the matrix-vector product and sum y = \alpha op(A) x + \beta y, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    pub fn zgemv(transA: ::enums::CblasTranspose, alpha: &gsl::Complex, A: &gsl::MatrixComplex, x: &gsl::VectorComplex, beta: &gsl::Complex,
        y: &mut gsl::VectorComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zgemv(transA, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            x.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(*beta), y.get_ffi()) }
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strmv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::MatrixFloat,
        x: &mut gsl::VectorFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_strmv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix_float, x.get_ffi()) }
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrmv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::Matrix,
        x: &mut gsl::Vector) -> i32 {
        unsafe { ::ffi::gsl_blas_dtrmv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix, x.get_ffi()) }
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrmv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::MatrixComplexFloat,
        x: &mut gsl::VectorComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ctrmv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix_complex_float, x.get_ffi()) }
    }

    /// This function computes the matrix-vector product x = op(A) x for the triangular matrix A, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrmv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::MatrixComplex,
        x: &mut gsl::VectorComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_ztrmv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix_complex, x.get_ffi()) }
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strsv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::MatrixFloat,
        x: &mut gsl::VectorFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_strsv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix_float, x.get_ffi()) }
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrsv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::Matrix,
        x: &mut gsl::Vector) -> i32 {
        unsafe { ::ffi::gsl_blas_dtrsv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix, x.get_ffi()) }
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrsv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::MatrixComplexFloat,
        x: &mut gsl::VectorComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ctrsv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix_complex_float, x.get_ffi()) }
    }

    /// This function computes inv(op(A)) x for x, where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans.
    /// When Uplo is CblasUpper then the upper triangle of A is used, and when Uplo is CblasLower then the lower triangle of A is used.
    /// If Diag is CblasNonUnit then the diagonal of the matrix is used, but if Diag is CblasUnit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrsv(uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag, A: &gsl::MatrixComplex,
        x: &mut gsl::VectorComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_ztrsv(uplo, transA, diag, A.get_ffi() as *const ::ffi::gsl_matrix_complex, x.get_ffi()) }
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssymv(uplo: ::enums::CblasUplo, alpha: f32, A: &gsl::MatrixFloat, x: &gsl::VectorFloat, beta: f32, y: &mut gsl::VectorFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ssymv(uplo, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float,
            x.get_ffi() as *const ::ffi::gsl_vector_float, beta, y.get_ffi()) }
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsymv(uplo: ::enums::CblasUplo, alpha: f64, A: &gsl::Matrix, x: &gsl::Vector, beta: f64, y: &mut gsl::Vector) -> i32 {
        unsafe { ::ffi::gsl_blas_dsymv(uplo, alpha, A.get_ffi() as *const ::ffi::gsl_matrix, x.get_ffi() as *const ::ffi::gsl_vector,
            beta, y.get_ffi()) }
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored. When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically assumed to be zero and are not referenced.
    pub fn chemv(uplo: ::enums::CblasUplo, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat, x: &gsl::VectorComplexFloat,
        beta: &gsl::ComplexFloat, y: &mut gsl::VectorComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_chemv(uplo, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            x.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(*beta), y.get_ffi()) }
    }

    /// These functions compute the matrix-vector product and sum y = \alpha A x + \beta y for the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored. When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically assumed to be zero and are not referenced.
    pub fn zhemv(uplo: ::enums::CblasUplo, alpha: &gsl::Complex, A: &gsl::MatrixComplex, x: &gsl::VectorComplex, beta: &gsl::Complex,
        y: &mut gsl::VectorComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zhemv(uplo, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            x.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(*beta), y.get_ffi()) }
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn sger(alpha: f32, x: &gsl::VectorFloat, y: &gsl::VectorFloat, A: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_sger(alpha, x.get_ffi() as *const ::ffi::gsl_vector_float,
            y.get_ffi() as *const ::ffi::gsl_vector_float, A.get_ffi()) }
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn dger(alpha: f64, x: &gsl::Vector, y: &gsl::Vector, A: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dger(alpha, x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi() as *const ::ffi::gsl_vector,
            A.get_ffi()) }
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn cgeru(alpha: &gsl::ComplexFloat, x: &gsl::VectorComplexFloat, y: &gsl::VectorComplexFloat, A: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cgeru(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
            y.get_ffi() as *const ::ffi::gsl_vector_complex_float, A.get_ffi()) }
    }

    /// This function computes the rank-1 update A = \alpha x y^T + A of the matrix A.
    pub fn zgeru(alpha: &gsl::Complex, x: &gsl::VectorComplex, y: &gsl::VectorComplex, A: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zgeru(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex,
            y.get_ffi() as *const ::ffi::gsl_vector_complex, A.get_ffi()) }
    }

    /// This function computes the conjugate rank-1 update A = \alpha x y^H + A of the matrix A.
    pub fn cgerc(alpha: &gsl::ComplexFloat, x: &gsl::VectorComplexFloat, y: &gsl::VectorComplexFloat, A: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cgerc(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
            y.get_ffi() as *const ::ffi::gsl_vector_complex_float, A.get_ffi()) }
    }

    /// This function computes the conjugate rank-1 update A = \alpha x y^H + A of the matrix A.
    pub fn zgerc(alpha: &gsl::Complex, x: &gsl::VectorComplex, y: &gsl::VectorComplex, A: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zgerc(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex,
            y.get_ffi() as *const ::ffi::gsl_vector_complex, A.get_ffi()) }
    }

    /// This function computes the symmetric rank-1 update A = \alpha x x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssyr(uplo: ::enums::CblasUplo, alpha: f32, x: &gsl::VectorFloat, A: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ssyr(uplo, alpha, x.get_ffi() as *const ::ffi::gsl_vector_float, A.get_ffi()) }
    }

    /// This function computes the symmetric rank-1 update A = \alpha x x^T + A of the symmetric matrix A. Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsyr(uplo: ::enums::CblasUplo, alpha: f64, x: &gsl::Vector, A: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dsyr(uplo, alpha, x.get_ffi() as *const ::ffi::gsl_vector, A.get_ffi()) }
    }

    /// These functions compute the hermitian rank-1 update A = \alpha x x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cher(uplo: ::enums::CblasUplo, alpha: f32, x: &gsl::VectorComplexFloat, A: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cher(uplo, alpha, x.get_ffi() as *const ::ffi::gsl_vector_complex_float, A.get_ffi()) }
    }

    /// These functions compute the hermitian rank-1 update A = \alpha x x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zher(uplo: ::enums::CblasUplo, alpha: f64, x: &gsl::VectorComplex, A: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zher(uplo, alpha, x.get_ffi() as *const ::ffi::gsl_vector_complex, A.get_ffi()) }
    }

    /// These functions compute the symmetric rank-2 update A = \alpha x y^T + \alpha y x^T + A of the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssyr2(uplo: ::enums::CblasUplo, alpha: f32, x: &gsl::VectorFloat, y: &gsl::VectorFloat, A: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ssyr2(uplo, alpha, x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi() as *const ::ffi::gsl_vector_float,
            A.get_ffi()) }
    }

    /// These functions compute the symmetric rank-2 update A = \alpha x y^T + \alpha y x^T + A of the symmetric matrix A.
    /// Since the matrix A is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsyr2(uplo: ::enums::CblasUplo, alpha: f64, x: &gsl::Vector, y: &gsl::Vector, A: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dsyr2(uplo, alpha, x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi() as *const ::ffi::gsl_vector,
            A.get_ffi()) }
    }

    /// These functions compute the hermitian rank-2 update A = \alpha x y^H + \alpha^* y x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cher2(uplo: ::enums::CblasUplo, alpha: &gsl::ComplexFloat, x: &gsl::VectorComplexFloat, y: &gsl::VectorComplexFloat,
        A: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cher2(uplo, ::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
            y.get_ffi() as *const ::ffi::gsl_vector_complex_float, A.get_ffi()) }
    }

    /// These functions compute the hermitian rank-2 update A = \alpha x y^H + \alpha^* y x^H + A of the hermitian matrix A.
    /// Since the matrix A is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zher2(uplo: ::enums::CblasUplo, alpha: &gsl::Complex, x: &gsl::VectorComplex, y: &gsl::VectorComplex,
        A: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zher2(uplo, ::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex,
            y.get_ffi() as *const ::ffi::gsl_vector_complex, A.get_ffi()) }
    }
}

pub mod level3 {
    use gsl;

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn sgemm(transA: ::enums::CblasTranspose, transB: ::enums::CblasTranspose, alpha: f32, A: &gsl::MatrixFloat,
        B: &gsl::MatrixFloat, beta: f32, C: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_sgemm(transA, transB, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float,
            B.get_ffi() as *const ::ffi::gsl_matrix_float, beta, C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn dgemm(transA: ::enums::CblasTranspose, transB: ::enums::CblasTranspose, alpha: f64, A: &gsl::Matrix, B: &gsl::Matrix,
        beta: f64, C: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dgemm(transA, transB, alpha, A.get_ffi() as *const ::ffi::gsl_matrix,
            B.get_ffi() as *const ::ffi::gsl_matrix, beta, C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn cgemm(transA: ::enums::CblasTranspose, transB: ::enums::CblasTranspose, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat,
        B: &gsl::MatrixComplexFloat, beta: &gsl::ComplexFloat, C: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cgemm(transA, transB, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex_float, ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
    pub fn zgemm(transA: ::enums::CblasTranspose, transB: ::enums::CblasTranspose, alpha: &gsl::Complex, A: &gsl::MatrixComplex,
        B: &gsl::MatrixComplex, beta: &gsl::Complex, C: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zgemm(transA, transB, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex, ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn ssymm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, alpha: f32, A: &gsl::MatrixFloat, B: &gsl::MatrixFloat, beta: f32,
        C: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ssymm(side, uplo, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float, B.get_ffi() as *const ::ffi::gsl_matrix_float,
            beta, C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn dsymm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, alpha: f64, A: &gsl::Matrix, B: &gsl::Matrix, beta: f64,
        C: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dsymm(side, uplo, alpha, A.get_ffi() as *const ::ffi::gsl_matrix, B.get_ffi() as *const ::ffi::gsl_matrix,
            beta, C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn csymm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat, B: &gsl::MatrixComplexFloat,
        beta: &gsl::ComplexFloat, C: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_csymm(side, uplo, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex_float, ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is symmetric.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    pub fn zsymm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, alpha: &gsl::Complex, A: &gsl::MatrixComplex, B: &gsl::MatrixComplex,
        beta: &gsl::Complex, C: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zsymm(side, uplo, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex, ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is Left and C = \alpha B A + \beta C for Side is Right, where the matrix A is hermitian.
    /// When Uplo is Upper then the upper triangle and diagonal of A are used, and when Uplo is Lower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn chemm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat, B: &gsl::MatrixComplexFloat,
        beta: &gsl::ComplexFloat, C: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_chemm(side, uplo, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex_float, ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product and sum C = \alpha A B + \beta C for Side is CblasLeft and C = \alpha B A + \beta C for Side is CblasRight, where the matrix A is hermitian.
    /// When Uplo is CblasUpper then the upper triangle and diagonal of A are used, and when Uplo is CblasLower then the lower triangle and diagonal of A are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zhemm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, alpha: &gsl::Complex, A: &gsl::MatrixComplex, B: &gsl::MatrixComplex,
        beta: &gsl::Complex, C: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zhemm(side, uplo, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex, ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strmm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: f32, A: &gsl::MatrixFloat, B: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_strmm(side, uplo, transA, diag, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float, B.get_ffi()) }
    }

    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrmm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: f64, A: &gsl::Matrix, B: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dtrmm(side, uplo, transA, diag, alpha, A.get_ffi() as *const ::ffi::gsl_matrix, B.get_ffi()) }
    }
    
    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrmm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat, B: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ctrmm(side, uplo, transA, diag, ::std::mem::transmute(*alpha),
            A.get_ffi() as *const ::ffi::gsl_matrix_complex_float, B.get_ffi()) }
    }
    
    /// This function computes the matrix-matrix product B = \alpha op(A) B for Side is Left and B = \alpha B op(A) for Side is CblasRight.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrmm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: &gsl::Complex, A: &gsl::MatrixComplex, B: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_ztrmm(side, uplo, transA, diag, ::std::mem::transmute(*alpha),
            A.get_ffi() as *const ::ffi::gsl_matrix_complex, B.get_ffi()) }
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn strsm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: f32, A: &gsl::MatrixFloat, B: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_strsm(side, uplo, transA, diag, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float, B.get_ffi()) }
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn dtrsm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: f64, A: &gsl::Matrix, B: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dtrsm(side, uplo, transA, diag, alpha, A.get_ffi() as *const ::ffi::gsl_matrix, B.get_ffi()) }
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ctrsm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat, B: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ctrsm(side, uplo, transA, diag, ::std::mem::transmute(*alpha),
            A.get_ffi() as *const ::ffi::gsl_matrix_complex_float, B.get_ffi()) }
    }

    /// This function computes the inverse-matrix matrix product B = \alpha op(inv(A))B for Side is Left and B = \alpha B op(inv(A)) for Side is Right.
    /// The matrix A is triangular and op(A) = A, A^T, A^H for TransA = NoTrans, Trans, ConjTrans.
    /// When Uplo is Upper then the upper triangle of A is used, and when Uplo is Lower then the lower triangle of A is used.
    /// If Diag is NonUnit then the diagonal of A is used, but if Diag is Unit then the diagonal elements of the matrix A are taken as unity and are not referenced.
    pub fn ztrsm(side: ::enums::CblasSide, uplo: ::enums::CblasUplo, transA: ::enums::CblasTranspose, diag: ::enums::CblasDiag,
        alpha: &gsl::Complex, A: &gsl::MatrixComplex, B: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_ztrsm(side, uplo, transA, diag, ::std::mem::transmute(*alpha),
            A.get_ffi() as *const ::ffi::gsl_matrix_complex, B.get_ffi()) }
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn ssyrk(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: f32, A: &gsl::MatrixFloat, beta: f32,
        C: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ssyrk(uplo, trans, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float, beta, C.get_ffi()) }
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn dsyrk(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: f64, A: &gsl::Matrix, beta: f64,
        C: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dsyrk(uplo, trans, alpha, A.get_ffi() as *const ::ffi::gsl_matrix, beta, C.get_ffi()) }
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn csyrk(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat,
        beta: &gsl::ComplexFloat, C: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_csyrk(uplo, trans, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes a rank-k update of the symmetric matrix C, C = \alpha A A^T + \beta C when Trans is NoTrans and C = \alpha A^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn zsyrk(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: &gsl::Complex, A: &gsl::MatrixComplex,
        beta: &gsl::Complex, C: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zsyrk(uplo, trans, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// These functions compute a rank-k update of the hermitian matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and C = \alpha A^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cherk(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: f32, A: &gsl::MatrixComplexFloat,
        beta: f32, C: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cherk(uplo, trans, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_complex_float, beta, C.get_ffi()) }
    }

    /// These functions compute a rank-k update of the hermitian matrix C, C = \alpha A A^H + \beta C when Trans is NoTrans and C = \alpha A^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zherk(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: f64, A: &gsl::MatrixComplex,
        beta: f64, C: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zherk(uplo, trans, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_complex, beta, C.get_ffi()) }
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn ssyr2k(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: f32, A: &gsl::MatrixFloat, B: &gsl::MatrixFloat,
        beta: f32,  C: &mut gsl::MatrixFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_ssyr2k(uplo, trans, alpha, A.get_ffi() as *const ::ffi::gsl_matrix_float,
            B.get_ffi() as *const ::ffi::gsl_matrix_float, beta, C.get_ffi()) }
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn dsyr2k(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: f64, A: &gsl::Matrix, B: &gsl::Matrix, beta: f64,
        C: &mut gsl::Matrix) -> i32 {
        unsafe { ::ffi::gsl_blas_dsyr2k(uplo, trans, alpha, A.get_ffi() as *const ::ffi::gsl_matrix,
            B.get_ffi() as *const ::ffi::gsl_matrix, beta, C.get_ffi()) }
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn csyr2k(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat,
        B: &gsl::MatrixComplexFloat, beta: &gsl::ComplexFloat, C: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_csyr2k(uplo, trans, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex_float, ::std::mem::transmute(*beta), C.get_ffi()) }
    }

    /// This function computes a rank-2k update of the symmetric matrix C, C = \alpha A B^T + \alpha B A^T + \beta C when Trans is NoTrans and C = \alpha A^T B + \alpha B^T A + \beta C when Trans is Trans.
    /// Since the matrix C is symmetric only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    pub fn zsyr2k(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: &gsl::Complex, A: &gsl::MatrixComplex, B: &gsl::MatrixComplex,
        beta: &gsl::Complex, C: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zsyr2k(uplo, trans, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex, ::std::mem::transmute(*beta), C.get_ffi()) }
    }
    
    /// This function computes a rank-2k update of the hermitian matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn cher2k(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: &gsl::ComplexFloat, A: &gsl::MatrixComplexFloat,
        B: &gsl::MatrixComplexFloat, beta: f32, C: &mut gsl::MatrixComplexFloat) -> i32 {
        unsafe { ::ffi::gsl_blas_cher2k(uplo, trans, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex_float,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex_float, beta, C.get_ffi()) }
    }

    /// This function computes a rank-2k update of the hermitian matrix C, C = \alpha A B^H + \alpha^* B A^H + \beta C when Trans is NoTrans and C = \alpha A^H B + \alpha^* B^H A + \beta C when Trans is ConjTrans.
    /// Since the matrix C is hermitian only its upper half or lower half need to be stored.
    /// When Uplo is Upper then the upper triangle and diagonal of C are used, and when Uplo is Lower then the lower triangle and diagonal of C are used.
    /// The imaginary elements of the diagonal are automatically set to zero.
    pub fn zher2k(uplo: ::enums::CblasUplo, trans: ::enums::CblasTranspose, alpha: &gsl::Complex, A: &gsl::MatrixComplex, B: &gsl::MatrixComplex,
        beta: f64, C: &mut gsl::MatrixComplex) -> i32 {
        unsafe { ::ffi::gsl_blas_zher2k(uplo, trans, ::std::mem::transmute(*alpha), A.get_ffi() as *const ::ffi::gsl_matrix_complex,
            B.get_ffi() as *const ::ffi::gsl_matrix_complex, beta, C.get_ffi()) }
    }
}