//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Linear Algebra

This chapter describes functions for solving linear systems. The library provides linear algebra operations which operate directly on the 
gsl_vector and gsl_matrix objects. These routines use the standard algorithms from Golub & Van Loan’s Matrix Computations with Level-1 and 
Level-2 BLAS calls for efficiency.

##LU Decomposition

A general square matrix A has an LU decomposition into upper and lower triangular matrices,

P A = L U
where P is a permutation matrix, L is unit lower triangular matrix and U is upper triangular matrix. For square matrices this decomposition 
can be used to convert the linear system A x = b into a pair of triangular systems (L y = P b, U x = y), which can be solved by forward and 
back-substitution. Note that the LU decomposition is valid for singular matrices.
!*/

use ffi;
use enums;

/// Factorise a general N x N matrix A into,
///
///  P A = L U
///
/// where P is a permutation matrix, L is unit lower triangular and U is upper triangular.
///
/// L is stored in the strict lower triangular part of the input matrix. The diagonal elements of L are unity and are not stored.
///
/// U is stored in the diagonal and upper triangular part of the input matrix.  
/// 
/// P is stored in the permutation p. Column j of P is column k of the identity matrix, where k = permutation->data[j]
///
/// signum gives the sign of the permutation, (-1)^n, where n is the  number of interchanges in the permutation. 
///
/// See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss Elimination with Partial Pivoting).
pub fn LU_decomp(a: &::MatrixF64, p: &::Permutation, signum: &mut i32) -> enums::Value {
    unsafe { ffi::gsl_linalg_LU_decomp(ffi::FFI::unwrap(a), ffi::FFI::unwrap(p), signum) }
}

/// Factorise a general N x N complex matrix A into,
///
///   P A = L U
///
/// where P is a permutation matrix, L is unit lower triangular and U is upper triangular.
///
/// L is stored in the strict lower triangular part of the input matrix. The diagonal elements of L are unity and are not stored.
///
/// U is stored in the diagonal and upper triangular part of the input matrix.  
/// 
/// P is stored in the permutation p. Column j of P is column k of the identity matrix, where k = permutation->data[j]
///
/// signum gives the sign of the permutation, (-1)^n, where n is the number of interchanges in the permutation. 
///
/// See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss Elimination with Partial Pivoting).
pub fn complex_LU_decomp(a: &::MatrixComplexF64, p: &::Permutation, signum: &mut i32) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_LU_decomp(ffi::FFI::unwrap(a), ffi::FFI::unwrap(p), signum) }
}

/// This function solves the square system A x = b using the LU decomposition of A into (LU, p) given by LU_decomp or LU_decomp as input.
pub fn LU_solve(lu: &::MatrixF64, p: &::Permutation, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_LU_solve(ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x)) }
}

/// This function solves the square system A x = b using the LU decomposition of A into (LU, p) given by LU_decomp or LU_decomp as input.
pub fn complex_LU_solve(lu: &::MatrixComplexF64, p: &::Permutation, b: &::VectorComplexF64, x: &::VectorComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_LU_solve(ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix_complex, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(b) as *const ffi::gsl_vector_complex, ffi::FFI::unwrap(x)) }
}

/// This function solves the square system A x = b in-place using the precomputed LU decomposition of A into (LU,p). On input x should contain
/// the right-hand side b, which is replaced by the solution on output.
pub fn LU_svx(lu: &::MatrixF64, p: &::Permutation, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_LU_svx(ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(x)) }
}

/// This function solves the square system A x = b in-place using the precomputed LU decomposition of A into (LU,p). On input x should contain
/// the right-hand side b, which is replaced by the solution on output.
pub fn complex_LU_svx(lu: &::MatrixComplexF64, p: &::Permutation, x: &::VectorComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_LU_svx(ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix_complex, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(x)) }
}

/// This function applies an iterative improvement to x, the solution of A x = b, from the precomputed LU decomposition of A into (LU,p). The
/// initial residual r = A x - b is also computed and stored in residual.
pub fn LU_refine(a: &::MatrixF64, lu: &::MatrixF64, p: &::Permutation, b: &::VectorF64, x: &::VectorF64, residual: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_LU_refine(ffi::FFI::unwrap(a) as *const ffi::gsl_matrix, ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix,
        ffi::FFI::unwrap(p) as *const ffi::gsl_permutation, ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x),
        ffi::FFI::unwrap(residual)) }
}

/// This function applies an iterative improvement to x, the solution of A x = b, from the precomputed LU decomposition of A into (LU,p). The
/// initial residual r = A x - b is also computed and stored in residual.
pub fn complex_LU_refine(a: &::MatrixComplexF64, lu: &::MatrixComplexF64, p: &::Permutation, b: &::VectorComplexF64, x: &::VectorComplexF64,
    residual: &::VectorComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_LU_refine(ffi::FFI::unwrap(a) as *const ffi::gsl_matrix_complex, ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix_complex,
        ffi::FFI::unwrap(p) as *const ffi::gsl_permutation, ffi::FFI::unwrap(b) as *const ffi::gsl_vector_complex, ffi::FFI::unwrap(x),
        ffi::FFI::unwrap(residual)) }
}

/// This function computes the inverse of a matrix A from its LU decomposition (LU,p), storing the result in the matrix inverse. The inverse
/// is computed by solving the system A x = b for each column of the identity matrix. It is preferable to avoid direct use of the inverse
/// whenever possible, as the linear solver functions can obtain the same result more efficiently and reliably (consult any introductory
/// textbook on numerical linear algebra for details).
pub fn LU_invert(lu: &::MatrixF64, p: &::Permutation, inverse: &::MatrixF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_LU_invert(ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(inverse)) }
}

/// This function computes the inverse of a matrix A from its LU decomposition (LU,p), storing the result in the matrix inverse. The inverse
/// is computed by solving the system A x = b for each column of the identity matrix. It is preferable to avoid direct use of the inverse
/// whenever possible, as the linear solver functions can obtain the same result more efficiently and reliably (consult any introductory
/// textbook on numerical linear algebra for details).
pub fn complex_LU_invert(lu: &::MatrixComplexF64, p: &::Permutation, inverse: &::MatrixComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_LU_invert(ffi::FFI::unwrap(lu) as *const ffi::gsl_matrix_complex, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(inverse)) }
}

/// This function computes the determinant of a matrix A from its LU decomposition, LU. The determinant is computed as the product of the
/// diagonal elements of U and the sign of the row permutation signum.
pub fn LU_det(lu: &::MatrixF64, signum: i32) -> f64 {
    unsafe { ffi::gsl_linalg_LU_det(ffi::FFI::unwrap(lu), signum) }
}

/// This function computes the determinant of a matrix A from its LU decomposition, LU. The determinant is computed as the product of the
/// diagonal elements of U and the sign of the row permutation signum.
pub fn complex_LU_det(lu: &::MatrixComplexF64, signum: i32) -> ::ComplexF64 {
    ::ComplexF64 {
        data: unsafe { ffi::gsl_linalg_complex_LU_det(ffi::FFI::unwrap(lu), signum).data }
    }
}

/// These functions compute the logarithm of the absolute value of the determinant of a matrix A, \ln|\det(A)|, from its LU decomposition,
/// LU. This function may be useful if the direct computation of the determinant would overflow or underflow.
pub fn LU_lndet(lu: &::MatrixF64) -> f64 {
    unsafe { ffi::gsl_linalg_LU_lndet(ffi::FFI::unwrap(lu)) }
}

/// These functions compute the logarithm of the absolute value of the determinant of a matrix A, \ln|\det(A)|, from its LU decomposition,
/// LU. This function may be useful if the direct computation of the determinant would overflow or underflow.
pub fn complex_LU_lndet(lu: &::MatrixComplexF64) -> f64 {
    unsafe { ffi::gsl_linalg_complex_LU_lndet(ffi::FFI::unwrap(lu)) }
}

/// This function computes the sign or phase factor of the determinant of a matrix A, \det(A)/|\det(A)|, from its LU decomposition, LU.
pub fn LU_sgndet(lu: &::MatrixF64, signum: i32) -> f64 {
    unsafe { ffi::gsl_linalg_LU_sgndet(ffi::FFI::unwrap(lu), signum) }
}

/// This function computes the sign or phase factor of the determinant of a matrix A, \det(A)/|\det(A)|, from its LU decomposition, LU.
pub fn complex_LU_sgndet(lu: &::MatrixComplexF64, signum: i32) -> ::ComplexF64 {
    ::ComplexF64 {
        data: unsafe { ffi::gsl_linalg_complex_LU_sgndet(ffi::FFI::unwrap(lu), signum).data }
    }
}