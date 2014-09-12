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

##QR Decomposition

A general rectangular M-by-N matrix A has a QR decomposition into the product of an orthogonal M-by-M square matrix Q (where Q^T Q = I) and 
an M-by-N right-triangular matrix R,

A = Q R
This decomposition can be used to convert the linear system A x = b into the triangular system R x = Q^T b, which can be solved by back-substitution. 
Another use of the QR decomposition is to compute an orthonormal basis for a set of vectors. The first N columns of Q form an orthonormal 
basis for the range of A, ran(A), when A has full column rank.

##QR Decomposition with Column Pivoting

The QR decomposition can be extended to the rank deficient case by introducing a column permutation P,

A P = Q R
The first r columns of Q form an orthonormal basis for the range of A for a matrix with column rank r. This decomposition can also be used 
to convert the linear system A x = b into the triangular system R y = Q^T b, x = P y, which can be solved by back-substitution and permutation. 
We denote the QR decomposition with column pivoting by QRP^T since A = Q R P^T.

##Singular Value Decomposition

A general rectangular M-by-N matrix A has a singular value decomposition (SVD) into the product of an M-by-N orthogonal matrix U, an N-by-N 
diagonal matrix of singular values S and the transpose of an N-by-N orthogonal square matrix V,

A = U S V^T

The singular values \sigma_i = S_{ii} are all non-negative and are generally chosen to form a non-increasing sequence \sigma_1 >= \sigma_2 >= 
... >= \sigma_N >= 0.

The singular value decomposition of a matrix has many practical uses. The condition number of the matrix is given by the ratio of the largest 
singular value to the smallest singular value. The presence of a zero singular value indicates that the matrix is singular. The number of 
non-zero singular values indicates the rank of the matrix. In practice singular value decomposition of a rank-deficient matrix will not produce 
exact zeroes for singular values, due to finite numerical precision. Small singular values should be edited by choosing a suitable tolerance.

For a rank-deficient matrix, the null space of A is given by the columns of V corresponding to the zero singular values. Similarly, the range 
of A is given by columns of U corresponding to the non-zero singular values.

Note that the routines here compute the “thin” version of the SVD with U as M-by-N orthogonal matrix. This allows in-place computation and is 
the most commonly-used form in practice. Mathematically, the “full” SVD is defined with U as an M-by-M orthogonal matrix and S as an M-by-N 
diagonal matrix (with additional rows of zeros).

##Cholesky Decomposition

A symmetric, positive definite square matrix A has a Cholesky decomposition into a product of a lower triangular matrix L and its transpose L^T,

A = L L^T

This is sometimes referred to as taking the square-root of a matrix. The Cholesky decomposition can only be carried out when all the eigenvalues 
of the matrix are positive. This decomposition can be used to convert the linear system A x = b into a pair of triangular systems (L y = b, 
L^T x = y), which can be solved by forward and back-substitution.

##Tridiagonal Decomposition of Real Symmetric Matrices

A symmetric matrix A can be factorized by similarity transformations into the form,

A = Q T Q^T

where Q is an orthogonal matrix and T is a symmetric tridiagonal matrix.
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

/// This function factorizes the M-by-N matrix A into the QR decomposition A = Q R. On output the diagonal and upper triangular part of the
/// input matrix contain the matrix R. The vector tau and the columns of the lower triangular part of the matrix A contain the Householder
/// coefficients and Householder vectors which encode the orthogonal matrix Q. The vector tau must be of length k=\min(M,N). The matrix Q
/// is related to these components by, Q = Q_k ... Q_2 Q_1 where Q_i = I - \tau_i v_i v_i^T and v_i is the Householder vector v_i =
/// (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i)). This is the same storage scheme as used by LAPACK.
/// 
/// The algorithm used to perform the decomposition is Householder QR (Golub & Van Loan, Matrix Computations, Algorithm 5.2.1).
pub fn QR_decomp(a: &::MatrixF64, tau: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_decomp(ffi::FFI::unwrap(a), ffi::FFI::unwrap(tau)) }
}

/// This function solves the square system A x = b using the QR decomposition of A held in (QR, tau) which must have been computed previously
/// with gsl_linalg_QR_decomp. The least-squares solution for rectangular systems can be found using QR_lssolve.
pub fn QR_solve(qr: &::MatrixF64, tau: &::VectorF64, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_solve(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x)) }
}

/// This function solves the square system A x = b in-place using the QR decomposition of A held in (QR,tau) which must have been computed
/// previously by gsl_linalg_QR_decomp. On input x should contain the right-hand side b, which is replaced by the solution on output.
pub fn QR_svx(qr: &::MatrixF64, tau: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_svx(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(x)) }
}

/// This function finds the least squares solution to the overdetermined system A x = b where the matrix A has more rows than columns. The
/// least squares solution minimizes the Euclidean norm of the residual, ||Ax - b||.The routine requires as input the QR decomposition of
/// A into (QR, tau) given by gsl_linalg_QR_decomp. The solution is returned in x. The residual is computed as a by-product and stored in
/// residual.
pub fn QR_lssolve(qr: &::MatrixF64, tau: &::VectorF64, b: &::VectorF64, x: &::VectorF64, residual: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_lssolve(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x), ffi::FFI::unwrap(residual)) }
}

/// This function applies the matrix Q^T encoded in the decomposition (QR,tau) to the vector v, storing the result Q^T v in v. The matrix
/// multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q^T.
pub fn QR_QTvec(qr: &::MatrixF64, tau: &::VectorF64, v: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_QTvec(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(v)) }
}

/// This function applies the matrix Q encoded in the decomposition (QR,tau) to the vector v, storing the result Q v in v. The matrix
/// multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q.
pub fn QR_Qvec(qr: &::MatrixF64, tau: &::VectorF64, v: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_Qvec(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(v)) }
}

/// This function applies the matrix Q^T encoded in the decomposition (QR,tau) to the matrix A, storing the result Q^T A in A. The matrix
/// multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q^T.
pub fn QR_QTmat(qr: &::MatrixF64, tau: &::VectorF64, v: &::MatrixF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_QTmat(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(v)) }
}

/// This function solves the triangular system R x = b for x. It may be useful if the product b' = Q^T b has already been computed using
/// gsl_linalg_QR_QTvec.
pub fn QR_Rsolve(qr: &::MatrixF64, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_Rsolve(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(b) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(x)) }
}

/// This function solves the triangular system R x = b for x in-place. On input x should contain the right-hand side b and is replaced by
/// the solution on output. This function may be useful if the product b' = Q^T b has already been computed using gsl_linalg_QR_QTvec.
pub fn QR_Rsvx(qr: &::MatrixF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_Rsvx(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(x)) }
}

/// This function unpacks the encoded QR decomposition (QR,tau) into the matrices Q and R, where Q is M-by-M and R is M-by-N.
pub fn QR_unpack(qr: &::MatrixF64, tau: &::VectorF64, q: &::MatrixF64, r: &::MatrixF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_unpack(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(q), ffi::FFI::unwrap(r)) }
}

/// This function solves the system R x = Q^T b for x. It can be used when the QR decomposition of a matrix is available in unpacked
/// form as (Q, R).
pub fn QR_QRsolve(q: &::MatrixF64, r: &::MatrixF64, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_QRsolve(ffi::FFI::unwrap(q), ffi::FFI::unwrap(r), ffi::FFI::unwrap(b) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(x)) }
}

/// This function performs a rank-1 update w v^T of the QR decomposition (Q, R). The update is given by Q'R' = Q (R + w v^T) where the
/// output matrices Q' and R' are also orthogonal and right triangular. Note that w is destroyed by the update.
pub fn QR_update(q: &::MatrixF64, r: &::MatrixF64, w: &::VectorF64, v: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QR_update(ffi::FFI::unwrap(q), ffi::FFI::unwrap(r), ffi::FFI::unwrap(w),
        ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
}

/// This function solves the triangular system R x = b for the N-by-N matrix R.
pub fn R_solve(r: &::MatrixF64, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_R_solve(ffi::FFI::unwrap(r) as *const ffi::gsl_matrix, ffi::FFI::unwrap(b) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(x)) }
}

/// This function solves the triangular system R x = b in-place. On input x should contain the right-hand side b, which is replaced by
/// the solution on output.
pub fn R_svx(r: &::MatrixF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_R_svx(ffi::FFI::unwrap(r) as *const ffi::gsl_matrix, ffi::FFI::unwrap(x)) }
}

/// This function factorizes the M-by-N matrix A into the QRP^T decomposition A = Q R P^T. On output the diagonal and upper triangular part
/// of the input matrix contain the matrix R. The permutation matrix P is stored in the permutation p. The sign of the permutation is given
/// by signum. It has the value (-1)^n, where n is the number of interchanges in the permutation. The vector tau and the columns of the lower
/// triangular part of the matrix A contain the Householder coefficients and vectors which encode the orthogonal matrix Q. The vector tau must
/// be of length k=\min(M,N). The matrix Q is related to these components by, Q = Q_k ... Q_2 Q_1 where Q_i = I - \tau_i v_i v_i^T and v_i is
/// the Householder vector v_i = (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i)). This is the same storage scheme as used by LAPACK. The vector norm is
/// a workspace of length N used for column pivoting.
/// 
/// The algorithm used to perform the decomposition is Householder QR with column pivoting (Golub & Van Loan, Matrix Computations, Algorithm 5.4.1).
pub fn QRPT_decomp(a: &::MatrixF64, tau: &::VectorF64, p: &::Permutation, signum: &mut i32, norm: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_decomp(ffi::FFI::unwrap(a), ffi::FFI::unwrap(tau), ffi::FFI::unwrap(p), signum,
        ffi::FFI::unwrap(norm)) }
}

/// This function factorizes the matrix A into the decomposition A = Q R P^T without modifying A itself and storing the output in the separate
/// matrices q and r.
pub fn QRPT_decomp2(a: &::MatrixF64, q: &::MatrixF64, r: &::MatrixF64, tau: &::VectorF64, p: &::Permutation, signum: &mut i32, norm: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_decomp2(ffi::FFI::unwrap(a) as *const ffi::gsl_matrix, ffi::FFI::unwrap(q), ffi::FFI::unwrap(r),
        ffi::FFI::unwrap(tau), ffi::FFI::unwrap(p), signum, ffi::FFI::unwrap(norm)) }
}

/// This function solves the square system A x = b using the QRP^T decomposition of A held in (QR, tau, p) which must have been computed previously
/// by QRPT_decomp.
pub fn QRPT_solve(qr: &::MatrixF64, tau: &::VectorF64, p: &::Permutation, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_solve(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(p) as *const ffi::gsl_permutation, ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x)) }
}

/// This function solves the square system A x = b in-place using the QRP^T decomposition of A held in (QR,tau,p). On input x should contain the
/// right-hand side b, which is replaced by the solution on output.
pub fn QRPT_svx(qr: &::MatrixF64, tau: &::VectorF64, p: &::Permutation, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_svx(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(p) as *const ffi::gsl_permutation, ffi::FFI::unwrap(x)) }
}

/// This function solves the square system R P^T x = Q^T b for x. It can be used when the QR decomposition of a matrix is available in unpacked
/// form as (Q, R).
pub fn QRPT_QRsolve(q: &::MatrixF64, r: &::MatrixF64, p: &::Permutation, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_QRsolve(ffi::FFI::unwrap(q) as *const ffi::gsl_matrix, ffi::FFI::unwrap(r) as *const ffi::gsl_matrix,
        ffi::FFI::unwrap(p) as *const ffi::gsl_permutation, ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x)) }
}

/// This function performs a rank-1 update w v^T of the QRP^T decomposition (Q, R, p). The update is given by Q'R' = Q (R + w v^T P) where the
/// output matrices Q' and R' are also orthogonal and right triangular. Note that w is destroyed by the update. The permutation p is not changed.
pub fn QRPT_update(q: &::MatrixF64, r: &::MatrixF64, p: &::Permutation, w: &::VectorF64, v: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_update(ffi::FFI::unwrap(q) as *const ffi::gsl_matrix, ffi::FFI::unwrap(r) as *const ffi::gsl_matrix,
        ffi::FFI::unwrap(p) as *const ffi::gsl_permutation, ffi::FFI::unwrap(w), ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
}

/// This function solves the triangular system R P^T x = b for the N-by-N matrix R contained in QR.
pub fn QRPT_Rsolve(qr: &::MatrixF64, p: &::Permutation, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_Rsolve(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x)) }
}

/// This function solves the triangular system R P^T x = b in-place for the N-by-N matrix R contained in QR. On input x should contain the
/// right-hand side b, which is replaced by the solution on output.
pub fn QRPT_Rsvx(qr: &::MatrixF64, p: &::Permutation, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_QRPT_Rsvx(ffi::FFI::unwrap(qr) as *const ffi::gsl_matrix, ffi::FFI::unwrap(p) as *const ffi::gsl_permutation,
        ffi::FFI::unwrap(x)) }
}

/// This function factorizes the M-by-N matrix A into the singular value decomposition A = U S V^T for M >= N. On output the matrix A is replaced
/// by U. The diagonal elements of the singular value matrix S are stored in the vector S. The singular values are non-negative and form a
/// non-increasing sequence from S_1 to S_N. The matrix V contains the elements of V in untransposed form. To form the product U S V^T it is
/// necessary to take the transpose of V. A workspace of length N is required in work.
/// 
/// This routine uses the Golub-Reinsch SVD algorithm.
pub fn SV_decomp(a: &::MatrixF64, v: &::MatrixF64, s: &::VectorF64, work: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_SV_decomp(ffi::FFI::unwrap(a), ffi::FFI::unwrap(v), ffi::FFI::unwrap(s), ffi::FFI::unwrap(work)) }
}

/// This function computes the SVD using the modified Golub-Reinsch algorithm, which is faster for M>>N. It requires the vector work of length
/// N and the N-by-N matrix X as additional working space.
pub fn SV_decomp_mod(a: &::MatrixF64, x: &::MatrixF64, v: &::MatrixF64, s: &::VectorF64, work: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_SV_decomp_mod(ffi::FFI::unwrap(a), ffi::FFI::unwrap(x), ffi::FFI::unwrap(v), ffi::FFI::unwrap(s),
        ffi::FFI::unwrap(work)) }
}

/// This function computes the SVD of the M-by-N matrix A using one-sided Jacobi orthogonalization for M >= N. The Jacobi method can compute
/// singular values to higher relative accuracy than Golub-Reinsch algorithms (see references for details).
pub fn SV_decomp_jacobi(a: &::MatrixF64, v: &::MatrixF64, s: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_SV_decomp_jacobi(ffi::FFI::unwrap(a), ffi::FFI::unwrap(v), ffi::FFI::unwrap(s)) }
}

/// This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously
/// with gsl_linalg_SV_decomp.
/// 
/// Only non-zero singular values are used in computing the solution. The parts of the solution corresponding to singular values of zero are
/// ignored. Other singular values can be edited out by setting them to zero before calling this function.
/// 
/// In the over-determined case where A has more rows than columns the system is solved in the least squares sense, returning the solution
/// x which minimizes ||A x - b||_2.
pub fn SV_solve(u: &::MatrixF64, v: &::MatrixF64, s: &::VectorF64, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_SV_solve(ffi::FFI::unwrap(u) as *const ffi::gsl_matrix, ffi::FFI::unwrap(v) as *const ffi::gsl_matrix,
        ffi::FFI::unwrap(s) as *const ffi::gsl_vector, ffi::FFI::unwrap(b) as *const ffi::gsl_vector, ffi::FFI::unwrap(x)) }
}

/// This function computes the statistical leverage values h_i of a matrix A using its singular value decomposition (U, S, V) previously computed
/// with gsl_linalg_SV_decomp. h_i are the diagonal values of the matrix A (A^T A)^{-1} A^T and depend only on the matrix U which is the input to
/// this function.
pub fn SV_leverage(u: &::MatrixF64, h: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_SV_leverage(ffi::FFI::unwrap(u) as *const ffi::gsl_matrix, ffi::FFI::unwrap(h)) }
}

/// This function factorizes the symmetric, positive-definite square matrix A into the Cholesky decomposition A = L L^T (or A = L L^H for
/// the complex case). On input, the values from the diagonal and lower-triangular part of the matrix A are used (the upper triangular part
/// is ignored). On output the diagonal and lower triangular part of the input matrix A contain the matrix L, while the upper triangular part
/// of the input matrix is overwritten with L^T (the diagonal terms being identical for both L and L^T). If the matrix is not positive-definite
/// then the decomposition will fail, returning the error code enums::Dom.
/// 
/// When testing whether a matrix is positive-definite, disable the error handler first to avoid triggering an error.
pub fn cholesky_decomp(a: &::MatrixF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_cholesky_decomp(ffi::FFI::unwrap(a)) }
}

/// This function factorizes the symmetric, positive-definite square matrix A into the Cholesky decomposition A = L L^T (or A = L L^H for
/// the complex case). On input, the values from the diagonal and lower-triangular part of the matrix A are used (the upper triangular part
/// is ignored). On output the diagonal and lower triangular part of the input matrix A contain the matrix L, while the upper triangular part
/// of the input matrix is overwritten with L^T (the diagonal terms being identical for both L and L^T). If the matrix is not positive-definite
/// then the decomposition will fail, returning the error code enums::Dom.
/// 
/// When testing whether a matrix is positive-definite, disable the error handler first to avoid triggering an error.
pub fn complex_cholesky_decomp(a: &::MatrixComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_cholesky_decomp(ffi::FFI::unwrap(a)) }
}

/// This function solves the system A x = b using the Cholesky decomposition of A held in the matrix cholesky which must have been previously
/// computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp.
pub fn cholesky_solve(cholesky: &::MatrixF64, b: &::VectorF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_cholesky_solve(ffi::FFI::unwrap(cholesky) as *const ffi::gsl_matrix, ffi::FFI::unwrap(b) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(x)) }
}

/// This function solves the system A x = b using the Cholesky decomposition of A held in the matrix cholesky which must have been previously
/// computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp.
pub fn complex_cholesky_solve(cholesky: &::MatrixComplexF64, b: &::VectorComplexF64, x: &::VectorComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_cholesky_solve(ffi::FFI::unwrap(cholesky) as *const ffi::gsl_matrix_complex,
        ffi::FFI::unwrap(b) as *const ffi::gsl_vector_complex, ffi::FFI::unwrap(x)) }
}

/// This function solves the system A x = b in-place using the Cholesky decomposition of A held in the matrix cholesky which must have been
/// previously computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On input x should contain the right-hand side
/// b, which is replaced by the solution on output.
pub fn cholesky_svx(cholesky: &::MatrixF64, x: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_cholesky_svx(ffi::FFI::unwrap(cholesky) as *const ffi::gsl_matrix, ffi::FFI::unwrap(x)) }
}

/// This function solves the system A x = b in-place using the Cholesky decomposition of A held in the matrix cholesky which must have been
/// previously computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On input x should contain the right-hand side
/// b, which is replaced by the solution on output.
pub fn complex_cholesky_svx(cholesky: &::MatrixComplexF64, x: &::VectorComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_cholesky_svx(ffi::FFI::unwrap(cholesky) as *const ffi::gsl_matrix_complex, ffi::FFI::unwrap(x)) }
}

/// This function computes the inverse of a matrix from its Cholesky decomposition cholesky, which must have been previously computed by
/// gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On output, the inverse is stored in-place in cholesky.
pub fn cholesky_invert(cholesky: &::MatrixF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_cholesky_invert(ffi::FFI::unwrap(cholesky)) }
}

/// This function computes the inverse of a matrix from its Cholesky decomposition cholesky, which must have been previously computed by
/// gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On output, the inverse is stored in-place in cholesky.
pub fn complex_cholesky_invert(cholesky: &::MatrixComplexF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_complex_cholesky_invert(ffi::FFI::unwrap(cholesky)) }
}

/// This function factorizes the symmetric square matrix A into the symmetric tridiagonal decomposition Q T Q^T. On output the diagonal and
/// subdiagonal part of the input matrix A contain the tridiagonal matrix T. The remaining lower triangular part of the input matrix contains
/// the Householder vectors which, together with the Householder coefficients tau, encode the orthogonal matrix Q. This storage scheme is
/// the same as used by LAPACK. The upper triangular part of A is not referenced.
pub fn symmtd_decomp(a: &::MatrixF64, tau: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_symmtd_decomp(ffi::FFI::unwrap(a), ffi::FFI::unwrap(tau)) }
}

/// This function unpacks the encoded symmetric tridiagonal decomposition (A, tau) obtained from gsl_linalg_symmtd_decomp into the orthogonal
/// matrix Q, the vector of diagonal elements diag and the vector of subdiagonal elements subdiag.
pub fn symmtd_unpack(a: &::MatrixF64, tau: &::VectorF64, q: &::MatrixF64, diag: &::VectorF64, subdiag: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_symmtd_unpack(ffi::FFI::unwrap(a) as *const ffi::gsl_matrix, ffi::FFI::unwrap(tau) as *const ffi::gsl_vector,
        ffi::FFI::unwrap(q), ffi::FFI::unwrap(diag), ffi::FFI::unwrap(subdiag)) }
}

/// This function unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal decomposition (A, tau) obtained from
/// gsl_linalg_symmtd_decomp into the vectors diag and subdiag.
pub fn symmtd_unpack_T(a: &::MatrixF64, diag: &::VectorF64, subdiag: &::VectorF64) -> enums::Value {
    unsafe { ffi::gsl_linalg_symmtd_unpack_T(ffi::FFI::unwrap(a) as *const ffi::gsl_matrix, ffi::FFI::unwrap(diag), ffi::FFI::unwrap(subdiag)) }
}