//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Linear Algebra

This chapter describes functions for solving linear systems. The library provides linear algebra operations which operate directly on the
gsl_vector and gsl_matrix objects. These routines use the standard algorithms from Golub & Van Loan’s Matrix Computations with Level-1 and
Level-2 BLAS calls for efficiency.

## LU Decomposition

A general square matrix A has an LU decomposition into upper and lower triangular matrices,

P A = L U
where P is a permutation matrix, L is unit lower triangular matrix and U is upper triangular matrix. For square matrices this decomposition
can be used to convert the linear system A x = b into a pair of triangular systems (L y = P b, U x = y), which can be solved by forward and
back-substitution. Note that the LU decomposition is valid for singular matrices.

## QR Decomposition

A general rectangular M-by-N matrix A has a QR decomposition into the product of an orthogonal M-by-M square matrix Q (where Q^T Q = I) and
an M-by-N right-triangular matrix R,

A = Q R
This decomposition can be used to convert the linear system A x = b into the triangular system R x = Q^T b, which can be solved by back-substitution.
Another use of the QR decomposition is to compute an orthonormal basis for a set of vectors. The first N columns of Q form an orthonormal
basis for the range of A, ran(A), when A has full column rank.

## QR Decomposition with Column Pivoting

The QR decomposition can be extended to the rank deficient case by introducing a column permutation P,

A P = Q R
The first r columns of Q form an orthonormal basis for the range of A for a matrix with column rank r. This decomposition can also be used
to convert the linear system A x = b into the triangular system R y = Q^T b, x = P y, which can be solved by back-substitution and permutation.
We denote the QR decomposition with column pivoting by QRP^T since A = Q R P^T.

## Singular Value Decomposition

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

## Cholesky Decomposition

A symmetric, positive definite square matrix A has a Cholesky decomposition into a product of a lower triangular matrix L and its transpose L^T,

A = L L^T

This is sometimes referred to as taking the square-root of a matrix. The Cholesky decomposition can only be carried out when all the eigenvalues
of the matrix are positive. This decomposition can be used to convert the linear system A x = b into a pair of triangular systems (L y = b,
L^T x = y), which can be solved by forward and back-substitution.

## Tridiagonal Decomposition of Real Symmetric Matrices

A symmetric matrix A can be factorized by similarity transformations into the form,

A = Q T Q^T

where Q is an orthogonal matrix and T is a symmetric tridiagonal matrix.

## Tridiagonal Decomposition of Hermitian Matrices

A hermitian matrix A can be factorized by similarity transformations into the form,

A = U T U^T

where U is a unitary matrix and T is a real symmetric tridiagonal matrix.

## Hessenberg Decomposition of Real Matrices

A general real matrix A can be decomposed by orthogonal similarity transformations into the form

A = U H U^T

where U is orthogonal and H is an upper Hessenberg matrix, meaning that it has zeros below the first subdiagonal. The Hessenberg reduction
is the first step in the Schur decomposition for the nonsymmetric eigenvalue problem, but has applications in other areas as well.

## Hessenberg-Triangular Decomposition of Real Matrices

A general real matrix pair (A, B) can be decomposed by orthogonal similarity transformations into the form

A = U H V^T
B = U R V^T

where U and V are orthogonal, H is an upper Hessenberg matrix, and R is upper triangular. The Hessenberg-Triangular reduction is the first
step in the generalized Schur decomposition for the generalized eigenvalue problem.

## Bidiagonalization

A general matrix A can be factorized by similarity transformations into the form,

A = U B V^T
where U and V are orthogonal matrices and B is a N-by-N bidiagonal matrix with non-zero entries only on the diagonal and superdiagonal. The
size of U is M-by-N and the size of V is N-by-N.

## Householder Transformations

A Householder transformation is a rank-1 modification of the identity matrix which can be used to zero out selected elements of a vector.
A Householder matrix P takes the form,

P = I - \tau v v^T

where v is a vector (called the Householder vector) and \tau = 2/(v^T v). The functions described in this section use the rank-1 structure
of the Householder matrix to create and apply Householder transformations efficiently.

## Tridiagonal Systems

The functions described in this section efficiently solve symmetric, non-symmetric and cyclic tridiagonal systems with minimal storage. Note
that the current implementations of these functions use a variant of Cholesky decomposition, so the tridiagonal matrix must be positive definite.
For non-positive definite matrices, the functions return the error code ::Sing.

## Balancing

The process of balancing a matrix applies similarity transformations to make the rows and columns have comparable norms. This is useful, for
example, to reduce roundoff errors in the solution of eigenvalue problems. Balancing a matrix A consists of replacing A with a similar matrix

A' = D^(-1) A D

where D is a diagonal matrix whose entries are powers of the floating point radix.

##14.16 References and Further Reading

Further information on the algorithms described in this section can be found in the following book,

G. H. Golub, C. F. Van Loan, Matrix Computations (3rd Ed, 1996), Johns Hopkins University Press, ISBN 0-8018-5414-8.
The LAPACK library is described in the following manual,

LAPACK Users’ Guide (Third Edition, 1999), Published by SIAM, ISBN 0-89871-447-8.
http://www.netlib.org/lapack

The LAPACK source code can be found at the website above, along with an online copy of the users guide.

The Modified Golub-Reinsch algorithm is described in the following paper,

T.F. Chan, “An Improved Algorithm for Computing the Singular Value Decomposition”, ACM Transactions on Mathematical Software, 8 (1982), pp 72–83.
The Jacobi algorithm for singular value decomposition is described in the following papers,

J.C. Nash, “A one-sided transformation method for the singular value decomposition and algebraic eigenproblem”, Computer Journal, Volume 18, Number
1 (1975), p 74–76
J.C. Nash and S. Shlien “Simple algorithms for the partial singular value decomposition”, Computer Journal, Volume 30 (1987), p 268–275.
James Demmel, Krešimir Veselić, “Jacobi’s Method is more accurate than QR”, Lapack Working Note 15 (LAWN-15), October 1989. Available from netlib,
http://www.netlib.org/lapack/ in the lawns or lawnspdf directories.
!*/

use crate::enums;
use crate::Value;
use ffi::FFI;

use types::complex::FFFI;

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
/// P is stored in the permutation p. Column j of P is column k of the identity matrix, where `k = permutation->data[j]`
///
/// signum gives the sign of the permutation, (-1)^n, where n is the  number of interchanges in the permutation.
///
/// See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss Elimination with Partial Pivoting).
#[doc(alias = "gsl_linalg_LU_decomp")]
pub fn LU_decomp(a: &mut ::MatrixF64, p: &mut ::Permutation, signum: &mut i32) -> Value {
    Value::from(unsafe { sys::gsl_linalg_LU_decomp(a.unwrap_unique(), p.unwrap_unique(), signum) })
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
/// P is stored in the permutation p. Column j of P is column k of the identity matrix, where `k = permutation->data[j]`
///
/// signum gives the sign of the permutation, (-1)^n, where n is the number of interchanges in the permutation.
///
/// See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss Elimination with Partial Pivoting).
#[doc(alias = "gsl_linalg_complex_LU_decomp")]
pub fn complex_LU_decomp(
    a: &mut ::MatrixComplexF64,
    p: &mut ::Permutation,
    signum: &mut i32,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_LU_decomp(a.unwrap_unique(), p.unwrap_unique(), signum)
    })
}

/// This function solves the square system A x = b using the LU decomposition of A into (LU, p) given by LU_decomp or LU_decomp as input.
#[doc(alias = "gsl_linalg_LU_solve")]
pub fn LU_solve(
    lu: &::MatrixF64,
    p: &::Permutation,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_LU_solve(
            lu.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the square system A x = b using the LU decomposition of A into (LU, p) given by LU_decomp or LU_decomp as input.
#[doc(alias = "gsl_linalg_complex_LU_solve")]
pub fn complex_LU_solve(
    lu: &::MatrixComplexF64,
    p: &::Permutation,
    b: &::VectorComplexF64,
    x: &mut ::VectorComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_LU_solve(
            lu.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the square system A x = b in-place using the precomputed LU decomposition of A into (LU,p). On input x should contain
/// the right-hand side b, which is replaced by the solution on output.
#[doc(alias = "gsl_linalg_LU_svx")]
pub fn LU_svx(lu: &::MatrixF64, p: &::Permutation, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_LU_svx(lu.unwrap_shared(), p.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function solves the square system A x = b in-place using the precomputed LU decomposition of A into (LU,p). On input x should contain
/// the right-hand side b, which is replaced by the solution on output.
#[doc(alias = "gsl_linalg_complex_LU_svx")]
pub fn complex_LU_svx(
    lu: &::MatrixComplexF64,
    p: &::Permutation,
    x: &mut ::VectorComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_LU_svx(lu.unwrap_shared(), p.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function applies an iterative improvement to x, the solution of A x = b, from the precomputed LU decomposition of A into (LU,p). The
/// initial residual r = A x - b is also computed and stored in residual.
#[doc(alias = "gsl_linalg_LU_refine")]
pub fn LU_refine(
    a: &::MatrixF64,
    lu: &::MatrixF64,
    p: &::Permutation,
    b: &::VectorF64,
    x: &mut ::VectorF64,
    residual: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_LU_refine(
            a.unwrap_shared(),
            lu.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
            residual.unwrap_unique(),
        )
    })
}

/// This function applies an iterative improvement to x, the solution of A x = b, from the precomputed LU decomposition of A into (LU,p). The
/// initial residual r = A x - b is also computed and stored in residual.
#[doc(alias = "gsl_linalg_complex_LU_refine")]
pub fn complex_LU_refine(
    a: &mut ::MatrixComplexF64,
    lu: &::MatrixComplexF64,
    p: &::Permutation,
    b: &::VectorComplexF64,
    x: &mut ::VectorComplexF64,
    residual: &mut ::VectorComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_LU_refine(
            a.unwrap_unique(),
            lu.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
            residual.unwrap_unique(),
        )
    })
}

/// This function computes the inverse of a matrix A from its LU decomposition (LU,p), storing the result in the matrix inverse. The inverse
/// is computed by solving the system A x = b for each column of the identity matrix. It is preferable to avoid direct use of the inverse
/// whenever possible, as the linear solver functions can obtain the same result more efficiently and reliably (consult any introductory
/// textbook on numerical linear algebra for details).
#[doc(alias = "gsl_linalg_LU_invert")]
pub fn LU_invert(lu: &::MatrixF64, p: &::Permutation, inverse: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_LU_invert(
            lu.unwrap_shared(),
            p.unwrap_shared(),
            inverse.unwrap_unique(),
        )
    })
}

/// This function computes the inverse of a matrix A from its LU decomposition (LU,p), storing the result in the matrix inverse. The inverse
/// is computed by solving the system A x = b for each column of the identity matrix. It is preferable to avoid direct use of the inverse
/// whenever possible, as the linear solver functions can obtain the same result more efficiently and reliably (consult any introductory
/// textbook on numerical linear algebra for details).
#[doc(alias = "gsl_linalg_complex_LU_invert")]
pub fn complex_LU_invert(
    lu: &::MatrixComplexF64,
    p: &::Permutation,
    inverse: &mut ::MatrixComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_LU_invert(
            lu.unwrap_shared(),
            p.unwrap_shared(),
            inverse.unwrap_unique(),
        )
    })
}

/// This function computes the determinant of a matrix A from its LU decomposition, LU. The determinant is computed as the product of the
/// diagonal elements of U and the sign of the row permutation signum.
#[doc(alias = "gsl_linalg_LU_det")]
pub fn LU_det(lu: &mut ::MatrixF64, signum: i32) -> f64 {
    unsafe { sys::gsl_linalg_LU_det(lu.unwrap_unique(), signum) }
}

/// This function computes the determinant of a matrix A from its LU decomposition, LU. The determinant is computed as the product of the
/// diagonal elements of U and the sign of the row permutation signum.
#[doc(alias = "gsl_linalg_complex_LU_det")]
pub fn complex_LU_det(lu: &mut ::MatrixComplexF64, signum: i32) -> ::ComplexF64 {
    unsafe { sys::gsl_linalg_complex_LU_det(lu.unwrap_unique(), signum).wrap() }
}

/// These functions compute the logarithm of the absolute value of the determinant of a matrix A, \ln|\det(A)|, from its LU decomposition,
/// LU. This function may be useful if the direct computation of the determinant would overflow or underflow.
#[doc(alias = "gsl_linalg_LU_lndet")]
pub fn LU_lndet(lu: &mut ::MatrixF64) -> f64 {
    unsafe { sys::gsl_linalg_LU_lndet(lu.unwrap_unique()) }
}

/// This function computes the sign or phase factor of the determinant of a matrix A, \det(A)/|\det(A)|, from its LU decomposition, LU.
#[doc(alias = "gsl_linalg_complex_LU_lndet")]
pub fn complex_LU_lndet(lu: &mut ::MatrixComplexF64) -> f64 {
    unsafe { sys::gsl_linalg_complex_LU_lndet(lu.unwrap_unique()) }
}

/// This function computes the sign or phase factor of the determinant of a matrix A, \det(A)/|\det(A)|, from its LU decomposition, LU.
#[doc(alias = "gsl_linalg_LU_sgndet")]
pub fn LU_sgndet(lu: &mut ::MatrixF64, signum: i32) -> i32 {
    unsafe { sys::gsl_linalg_LU_sgndet(lu.unwrap_unique(), signum) }
}

/// This function computes the sign or phase factor of the determinant of a matrix A, \det(A)/|\det(A)|, from its LU decomposition, LU.
#[doc(alias = "gsl_linalg_complex_LU_sgndet")]
pub fn complex_LU_sgndet(lu: &mut ::MatrixComplexF64, signum: i32) -> ::ComplexF64 {
    unsafe { sys::gsl_linalg_complex_LU_sgndet(lu.unwrap_unique(), signum).wrap() }
}

/// This function factorizes the M-by-N matrix A into the QR decomposition A = Q R. On output the diagonal and upper triangular part of the
/// input matrix contain the matrix R. The vector tau and the columns of the lower triangular part of the matrix A contain the Householder
/// coefficients and Householder vectors which encode the orthogonal matrix Q. The vector tau must be of length k=\min(M,N). The matrix Q
/// is related to these components by, Q = Q_k ... Q_2 Q_1 where Q_i = I - \tau_i v_i v_i^T and v_i is the Householder vector v_i =
/// (0,...,1,A(i+1,i),A(i+2,i),...,A(m,i)). This is the same storage scheme as used by LAPACK.
///
/// The algorithm used to perform the decomposition is Householder QR (Golub & Van Loan, Matrix Computations, Algorithm 5.2.1).
#[doc(alias = "gsl_linalg_QR_decomp")]
pub fn QR_decomp(a: &mut ::MatrixF64, tau: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_QR_decomp(a.unwrap_unique(), tau.unwrap_unique()) })
}

/// This function solves the square system A x = b using the QR decomposition of A held in (QR, tau) which must have been computed previously
/// with gsl_linalg_QR_decomp. The least-squares solution for rectangular systems can be found using QR_lssolve.
#[doc(alias = "gsl_linalg_QR_solve")]
pub fn QR_solve(
    qr: &::MatrixF64,
    tau: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_solve(
            qr.unwrap_shared(),
            tau.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the square system A x = b in-place using the QR decomposition of A held in (QR,tau) which must have been computed
/// previously by gsl_linalg_QR_decomp. On input x should contain the right-hand side b, which is replaced by the solution on output.
#[doc(alias = "gsl_linalg_QR_svx")]
pub fn QR_svx(qr: &::MatrixF64, tau: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_svx(qr.unwrap_shared(), tau.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function finds the least squares solution to the overdetermined system A x = b where the matrix A has more rows than columns. The
/// least squares solution minimizes the Euclidean norm of the residual, ||Ax - b||.The routine requires as input the QR decomposition of
/// A into (QR, tau) given by gsl_linalg_QR_decomp. The solution is returned in x. The residual is computed as a by-product and stored in
/// residual.
#[doc(alias = "gsl_linalg_QR_lssolve")]
pub fn QR_lssolve(
    qr: &::MatrixF64,
    tau: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
    residual: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_lssolve(
            qr.unwrap_shared(),
            tau.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
            residual.unwrap_unique(),
        )
    })
}

/// This function applies the matrix Q^T encoded in the decomposition (QR,tau) to the vector v, storing the result Q^T v in v. The matrix
/// multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q^T.
#[doc(alias = "gsl_linalg_QR_QTvec")]
pub fn QR_QTvec(qr: &::MatrixF64, tau: &::VectorF64, v: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_QTvec(qr.unwrap_shared(), tau.unwrap_shared(), v.unwrap_unique())
    })
}

/// This function applies the matrix Q encoded in the decomposition (QR,tau) to the vector v, storing the result Q v in v. The matrix
/// multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q.
#[doc(alias = "gsl_linalg_QR_Qvec")]
pub fn QR_Qvec(qr: &::MatrixF64, tau: &::VectorF64, v: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_Qvec(qr.unwrap_shared(), tau.unwrap_shared(), v.unwrap_unique())
    })
}

/// This function applies the matrix Q^T encoded in the decomposition (QR,tau) to the matrix A, storing the result Q^T A in A. The matrix
/// multiplication is carried out directly using the encoding of the Householder vectors without needing to form the full matrix Q^T.
#[doc(alias = "gsl_linalg_QR_QTmat")]
pub fn QR_QTmat(qr: &::MatrixF64, tau: &::VectorF64, v: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_QTmat(qr.unwrap_shared(), tau.unwrap_shared(), v.unwrap_unique())
    })
}

/// This function solves the triangular system R x = b for x. It may be useful if the product b' = Q^T b has already been computed using
/// gsl_linalg_QR_QTvec.
#[doc(alias = "gsl_linalg_QR_Rsolve")]
pub fn QR_Rsolve(qr: &::MatrixF64, b: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_Rsolve(qr.unwrap_shared(), b.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function solves the triangular system R x = b for x in-place. On input x should contain the right-hand side b and is replaced by
/// the solution on output. This function may be useful if the product b' = Q^T b has already been computed using gsl_linalg_QR_QTvec.
#[doc(alias = "gsl_linalg_QR_Rsvx")]
pub fn QR_Rsvx(qr: &::MatrixF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_QR_Rsvx(qr.unwrap_shared(), x.unwrap_unique()) })
}

/// This function unpacks the encoded QR decomposition (QR,tau) into the matrices Q and R, where Q is M-by-M and R is M-by-N.
#[doc(alias = "gsl_linalg_QR_unpack")]
pub fn QR_unpack(
    qr: &::MatrixF64,
    tau: &::VectorF64,
    q: &mut ::MatrixF64,
    r: &mut ::MatrixF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_unpack(
            qr.unwrap_shared(),
            tau.unwrap_shared(),
            q.unwrap_unique(),
            r.unwrap_unique(),
        )
    })
}

/// This function solves the system R x = Q^T b for x. It can be used when the QR decomposition of a matrix is available in unpacked
/// form as (Q, R).
#[doc(alias = "gsl_linalg_QR_QRsolve")]
pub fn QR_QRsolve(
    q: &mut ::MatrixF64,
    r: &mut ::MatrixF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_QRsolve(
            q.unwrap_unique(),
            r.unwrap_unique(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function performs a rank-1 update w v^T of the QR decomposition (Q, R). The update is given by Q'R' = Q (R + w v^T) where the
/// output matrices Q' and R' are also orthogonal and right triangular. Note that w is destroyed by the update.
#[doc(alias = "gsl_linalg_QR_update")]
pub fn QR_update(
    q: &mut ::MatrixF64,
    r: &mut ::MatrixF64,
    mut w: ::VectorF64,
    v: &::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QR_update(
            q.unwrap_unique(),
            r.unwrap_unique(),
            w.unwrap_unique(),
            v.unwrap_shared(),
        )
    })
}

/// This function solves the triangular system R x = b for the N-by-N matrix R.
#[doc(alias = "gsl_linalg_R_solve")]
pub fn R_solve(r: &::MatrixF64, b: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_R_solve(r.unwrap_shared(), b.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function solves the triangular system R x = b in-place. On input x should contain the right-hand side b, which is replaced by
/// the solution on output.
#[doc(alias = "gsl_linalg_R_svx")]
pub fn R_svx(r: &::MatrixF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_R_svx(r.unwrap_shared(), x.unwrap_unique()) })
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
#[doc(alias = "gsl_linalg_QRPT_decomp")]
pub fn QRPT_decomp(
    a: &mut ::MatrixF64,
    tau: &mut ::VectorF64,
    p: &mut ::Permutation,
    signum: &mut i32,
    norm: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_decomp(
            a.unwrap_unique(),
            tau.unwrap_unique(),
            p.unwrap_unique(),
            signum,
            norm.unwrap_unique(),
        )
    })
}

/// This function factorizes the matrix A into the decomposition A = Q R P^T without modifying A itself and storing the output in the separate
/// matrices q and r.
#[doc(alias = "gsl_linalg_QRPT_decomp2")]
pub fn QRPT_decomp2(
    a: &::MatrixF64,
    q: &mut ::MatrixF64,
    r: &mut ::MatrixF64,
    tau: &mut ::VectorF64,
    p: &mut ::Permutation,
    signum: &mut i32,
    norm: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_decomp2(
            a.unwrap_shared(),
            q.unwrap_unique(),
            r.unwrap_unique(),
            tau.unwrap_unique(),
            p.unwrap_unique(),
            signum,
            norm.unwrap_unique(),
        )
    })
}

/// This function solves the square system A x = b using the QRP^T decomposition of A held in (QR, tau, p) which must have been computed previously
/// by QRPT_decomp.
#[doc(alias = "gsl_linalg_QRPT_solve")]
pub fn QRPT_solve(
    qr: &::MatrixF64,
    tau: &::VectorF64,
    p: &::Permutation,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_solve(
            qr.unwrap_shared(),
            tau.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the square system A x = b in-place using the QRP^T decomposition of A held in (QR,tau,p). On input x should contain the
/// right-hand side b, which is replaced by the solution on output.
#[doc(alias = "gsl_linalg_QRPT_svx")]
pub fn QRPT_svx(
    qr: &::MatrixF64,
    tau: &::VectorF64,
    p: &::Permutation,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_svx(
            qr.unwrap_shared(),
            tau.unwrap_shared(),
            p.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the square system R P^T x = Q^T b for x. It can be used when the QR decomposition of a matrix is available in unpacked
/// form as (Q, R).
#[doc(alias = "gsl_linalg_QRPT_QRsolve")]
pub fn QRPT_QRsolve(
    q: &::MatrixF64,
    r: &::MatrixF64,
    p: &::Permutation,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_QRsolve(
            q.unwrap_shared(),
            r.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function performs a rank-1 update w v^T of the QRP^T decomposition (Q, R, p). The update is given by Q'R' = Q (R + w v^T P) where the
/// output matrices Q' and R' are also orthogonal and right triangular. Note that w is destroyed by the update. The permutation p is not changed.
#[doc(alias = "gsl_linalg_QRPT_update")]
pub fn QRPT_update(
    q: &mut ::MatrixF64,
    r: &mut ::MatrixF64,
    p: &::Permutation,
    w: &mut ::VectorF64,
    v: &::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_update(
            q.unwrap_unique(),
            r.unwrap_unique(),
            p.unwrap_shared(),
            w.unwrap_unique(),
            v.unwrap_shared(),
        )
    })
}

/// This function solves the triangular system R P^T x = b for the N-by-N matrix R contained in QR.
#[doc(alias = "gsl_linalg_QRPT_Rsolve")]
pub fn QRPT_Rsolve(
    qr: &::MatrixF64,
    p: &::Permutation,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_Rsolve(
            qr.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the triangular system R P^T x = b in-place for the N-by-N matrix R contained in QR. On input x should contain the
/// right-hand side b, which is replaced by the solution on output.
#[doc(alias = "gsl_linalg_QRPT_Rsvx")]
pub fn QRPT_Rsvx(qr: &::MatrixF64, p: &::Permutation, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_QRPT_Rsvx(qr.unwrap_shared(), p.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function factorizes the M-by-N matrix A into the singular value decomposition A = U S V^T for M >= N. On output the matrix A is replaced
/// by U. The diagonal elements of the singular value matrix S are stored in the vector S. The singular values are non-negative and form a
/// non-increasing sequence from S_1 to S_N. The matrix V contains the elements of V in untransposed form. To form the product U S V^T it is
/// necessary to take the transpose of V. A workspace of length N is required in work.
///
/// This routine uses the Golub-Reinsch SVD algorithm.
#[doc(alias = "gsl_linalg_SV_decomp")]
pub fn SV_decomp(
    a: &mut ::MatrixF64,
    v: &mut ::MatrixF64,
    s: &mut ::VectorF64,
    work: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_SV_decomp(
            a.unwrap_unique(),
            v.unwrap_unique(),
            s.unwrap_unique(),
            work.unwrap_unique(),
        )
    })
}

/// This function computes the SVD using the modified Golub-Reinsch algorithm, which is faster for M>>N. It requires the vector work of length
/// N and the N-by-N matrix X as additional working space.
#[doc(alias = "gsl_linalg_SV_decomp_mod")]
pub fn SV_decomp_mod(
    a: &mut ::MatrixF64,
    x: &mut ::MatrixF64,
    v: &mut ::MatrixF64,
    s: &mut ::VectorF64,
    work: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_SV_decomp_mod(
            a.unwrap_unique(),
            x.unwrap_unique(),
            v.unwrap_unique(),
            s.unwrap_unique(),
            work.unwrap_unique(),
        )
    })
}

/// This function computes the SVD of the M-by-N matrix A using one-sided Jacobi orthogonalization for M >= N. The Jacobi method can compute
/// singular values to higher relative accuracy than Golub-Reinsch algorithms (see references for details).
#[doc(alias = "gsl_linalg_SV_decomp_jacobi")]
pub fn SV_decomp_jacobi(a: &mut ::MatrixF64, v: &mut ::MatrixF64, s: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_SV_decomp_jacobi(a.unwrap_unique(), v.unwrap_unique(), s.unwrap_unique())
    })
}

/// This function solves the system A x = b using the singular value decomposition (U, S, V) of A which must have been computed previously
/// with gsl_linalg_SV_decomp.
///
/// Only non-zero singular values are used in computing the solution. The parts of the solution corresponding to singular values of zero are
/// ignored. Other singular values can be edited out by setting them to zero before calling this function.
///
/// In the over-determined case where A has more rows than columns the system is solved in the least squares sense, returning the solution
/// x which minimizes ||A x - b||_2.
#[doc(alias = "gsl_linalg_SV_solve")]
pub fn SV_solve(
    u: &::MatrixF64,
    v: &::MatrixF64,
    s: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_SV_solve(
            u.unwrap_shared(),
            v.unwrap_shared(),
            s.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function computes the statistical leverage values h_i of a matrix A using its singular value decomposition (U, S, V) previously computed
/// with gsl_linalg_SV_decomp. h_i are the diagonal values of the matrix A (A^T A)^{-1} A^T and depend only on the matrix U which is the input to
/// this function.
#[doc(alias = "gsl_linalg_SV_leverage")]
pub fn SV_leverage(u: &::MatrixF64, h: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_SV_leverage(u.unwrap_shared(), h.unwrap_unique()) })
}

/// This function factorizes the symmetric, positive-definite square matrix A into the Cholesky decomposition A = L L^T (or A = L L^H for
/// the complex case). On input, the values from the diagonal and lower-triangular part of the matrix A are used (the upper triangular part
/// is ignored). On output the diagonal and lower triangular part of the input matrix A contain the matrix L, while the upper triangular part
/// of the input matrix is overwritten with L^T (the diagonal terms being identical for both L and L^T). If the matrix is not positive-definite
/// then the decomposition will fail, returning the error code ::Dom.
///
/// When testing whether a matrix is positive-definite, disable the error handler first to avoid triggering an error.
#[doc(alias = "gsl_linalg_cholesky_decomp")]
pub fn cholesky_decomp(a: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_cholesky_decomp(a.unwrap_unique()) })
}

/// This function factorizes the symmetric, positive-definite square matrix A into the Cholesky decomposition A = L L^T (or A = L L^H for
/// the complex case). On input, the values from the diagonal and lower-triangular part of the matrix A are used (the upper triangular part
/// is ignored). On output the diagonal and lower triangular part of the input matrix A contain the matrix L, while the upper triangular part
/// of the input matrix is overwritten with L^T (the diagonal terms being identical for both L and L^T). If the matrix is not positive-definite
/// then the decomposition will fail, returning the error code ::Dom.
///
/// When testing whether a matrix is positive-definite, disable the error handler first to avoid triggering an error.
#[doc(alias = "gsl_linalg_complex_cholesky_decomp")]
pub fn complex_cholesky_decomp(a: &mut ::MatrixComplexF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_complex_cholesky_decomp(a.unwrap_unique()) })
}

/// This function solves the system A x = b using the Cholesky decomposition of A held in the matrix cholesky which must have been previously
/// computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp.
#[doc(alias = "gsl_linalg_cholesky_solve")]
pub fn cholesky_solve(cholesky: &::MatrixF64, b: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_solve(
            cholesky.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the system A x = b using the Cholesky decomposition of A held in the matrix cholesky which must have been previously
/// computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp.
#[doc(alias = "gsl_linalg_complex_cholesky_solve")]
pub fn complex_cholesky_solve(
    cholesky: &::MatrixComplexF64,
    b: &::VectorComplexF64,
    x: &mut ::VectorComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_cholesky_solve(
            cholesky.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the system A x = b in-place using the Cholesky decomposition of A held in the matrix cholesky which must have been
/// previously computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On input x should contain the right-hand side
/// b, which is replaced by the solution on output.
#[doc(alias = "gsl_linalg_cholesky_svx")]
pub fn cholesky_svx(cholesky: &::MatrixF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_svx(cholesky.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function solves the system A x = b in-place using the Cholesky decomposition of A held in the matrix cholesky which must have been
/// previously computed by gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On input x should contain the right-hand side
/// b, which is replaced by the solution on output.
#[doc(alias = "gsl_linalg_complex_cholesky_svx")]
pub fn complex_cholesky_svx(cholesky: &::MatrixComplexF64, x: &mut ::VectorComplexF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_cholesky_svx(cholesky.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function computes the inverse of a matrix from its Cholesky decomposition cholesky, which must have been previously computed by
/// gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On output, the inverse is stored in-place in cholesky.
#[doc(alias = "gsl_linalg_cholesky_invert")]
pub fn cholesky_invert(cholesky: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_cholesky_invert(cholesky.unwrap_unique()) })
}

/// This function computes the inverse of a matrix from its Cholesky decomposition cholesky, which must have been previously computed by
/// gsl_linalg_cholesky_decomp or gsl_linalg_complex_cholesky_decomp. On output, the inverse is stored in-place in cholesky.
#[doc(alias = "gsl_linalg_complex_cholesky_invert")]
pub fn complex_cholesky_invert(cholesky: &mut ::MatrixComplexF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_complex_cholesky_invert(cholesky.unwrap_unique()) })
}

/// This function factorizes the symmetric square matrix A into the symmetric tridiagonal decomposition Q T Q^T. On output the diagonal and
/// subdiagonal part of the input matrix A contain the tridiagonal matrix T. The remaining lower triangular part of the input matrix contains
/// the Householder vectors which, together with the Householder coefficients tau, encode the orthogonal matrix Q. This storage scheme is
/// the same as used by LAPACK. The upper triangular part of A is not referenced.
#[doc(alias = "gsl_linalg_symmtd_decomp")]
pub fn symmtd_decomp(a: &mut ::MatrixF64, tau: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_symmtd_decomp(a.unwrap_unique(), tau.unwrap_unique()) })
}

/// This function unpacks the encoded symmetric tridiagonal decomposition (A, tau) obtained from gsl_linalg_symmtd_decomp into the orthogonal
/// matrix Q, the vector of diagonal elements diag and the vector of subdiagonal elements subdiag.
#[doc(alias = "gsl_linalg_symmtd_unpack")]
pub fn symmtd_unpack(
    a: &::MatrixF64,
    tau: &::VectorF64,
    q: &mut ::MatrixF64,
    diag: &mut ::VectorF64,
    subdiag: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_symmtd_unpack(
            a.unwrap_shared(),
            tau.unwrap_shared(),
            q.unwrap_unique(),
            diag.unwrap_unique(),
            subdiag.unwrap_unique(),
        )
    })
}

/// This function unpacks the diagonal and subdiagonal of the encoded symmetric tridiagonal decomposition (A, tau) obtained from
/// gsl_linalg_symmtd_decomp into the vectors diag and subdiag.
#[doc(alias = "gsl_linalg_symmtd_unpack_T")]
pub fn symmtd_unpack_T(
    a: &::MatrixF64,
    diag: &mut ::VectorF64,
    subdiag: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_symmtd_unpack_T(
            a.unwrap_shared(),
            diag.unwrap_unique(),
            subdiag.unwrap_unique(),
        )
    })
}

/// This function factorizes the hermitian matrix A into the symmetric tridiagonal decomposition U T U^T. On output the real parts of the
/// diagonal and subdiagonal part of the input matrix A contain the tridiagonal matrix T. The remaining lower triangular part of the input
/// matrix contains the Householder vectors which, together with the Householder coefficients tau, encode the unitary matrix U. This storage
/// scheme is the same as used by LAPACK. The upper triangular part of A and imaginary parts of the diagonal are not referenced.
#[doc(alias = "gsl_linalg_hermtd_decomp")]
pub fn hermtd_decomp(a: &mut ::MatrixComplexF64, tau: &mut ::VectorComplexF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_hermtd_decomp(a.unwrap_unique(), tau.unwrap_unique()) })
}

/// This function unpacks the encoded tridiagonal decomposition (A, tau) obtained from gsl_linalg_hermtd_decomp into the unitary matrix U,
/// the real vector of diagonal elements diag and the real vector of subdiagonal elements subdiag.
#[doc(alias = "gsl_linalg_hermtd_unpack")]
pub fn hermtd_unpack(
    a: &::MatrixComplexF64,
    tau: &::VectorComplexF64,
    u: &mut ::MatrixComplexF64,
    diag: &mut ::VectorF64,
    subdiag: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_hermtd_unpack(
            a.unwrap_shared(),
            tau.unwrap_shared(),
            u.unwrap_unique(),
            diag.unwrap_unique(),
            subdiag.unwrap_unique(),
        )
    })
}

/// This function unpacks the diagonal and subdiagonal of the encoded tridiagonal decomposition (A, tau) obtained from the
/// gsl_linalg_hermtd_decomp into the real vectors diag and subdiag.
#[doc(alias = "gsl_linalg_hermtd_unpack_T")]
pub fn hermtd_unpack_T(
    a: &::MatrixComplexF64,
    diag: &mut ::VectorF64,
    subdiag: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_hermtd_unpack_T(
            a.unwrap_shared(),
            diag.unwrap_unique(),
            subdiag.unwrap_unique(),
        )
    })
}

/// This function computes the Hessenberg decomposition of the matrix A by applying the similarity transformation H = U^T A U. On output, H
/// is stored in the upper portion of A. The information required to construct the matrix U is stored in the lower triangular portion of A.
/// U is a product of N - 2 Householder matrices. The Householder vectors are stored in the lower portion of A (below the subdiagonal) and
/// the Householder coefficients are stored in the vector tau. tau must be of length N.
#[doc(alias = "gsl_linalg_hessenberg_decomp")]
pub fn hessenberg_decomp(a: &mut ::MatrixF64, tau: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_hessenberg_decomp(a.unwrap_unique(), tau.unwrap_unique())
    })
}

/// This function constructs the orthogonal matrix U from the information stored in the Hessenberg matrix H along with the vector tau. H and
/// tau are outputs from gsl_linalg_hessenberg_decomp.
#[doc(alias = "gsl_linalg_hessenberg_unpack")]
pub fn hessenberg_unpack(h: &mut ::MatrixF64, tau: &mut ::VectorF64, u: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_hessenberg_unpack(h.unwrap_unique(), tau.unwrap_unique(), u.unwrap_unique())
    })
}

/// This function is similar to gsl_linalg_hessenberg_unpack, except it accumulates the matrix U into V, so that V' = VU. The matrix V must
/// be initialized prior to calling this function. Setting V to the identity matrix provides the same result as gsl_linalg_hessenberg_unpack.
/// If H is order N, then V must have N columns but may have any number of rows.
#[doc(alias = "gsl_linalg_hessenberg_unpack_accum")]
pub fn hessenberg_unpack_accum(
    h: &mut ::MatrixF64,
    tau: &mut ::VectorF64,
    v: &mut ::MatrixF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_hessenberg_unpack_accum(
            h.unwrap_unique(),
            tau.unwrap_unique(),
            v.unwrap_unique(),
        )
    })
}

/// This function sets the lower triangular portion of H, below the subdiagonal, to zero. It is useful for clearing out the Householder
/// vectors after calling gsl_linalg_hessenberg_decomp.
#[doc(alias = "gsl_linalg_hessenberg_set_zero")]
pub fn hessenberg_set_zero(h: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_hessenberg_set_zero(h.unwrap_unique()) })
}

/// This function computes the Hessenberg-Triangular decomposition of the matrix pair (A, B). On output, H is stored in A, and R is stored
/// in B. If U and V are provided (they may be null), the similarity transformations are stored in them. Additional workspace of length N
/// is needed in work.
#[doc(alias = "gsl_linalg_hesstri_decomp")]
pub fn hesstri_decomp(
    a: &mut ::MatrixF64,
    b: &mut ::MatrixF64,
    u: &mut ::MatrixF64,
    v: &mut ::MatrixF64,
    work: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_hesstri_decomp(
            a.unwrap_unique(),
            b.unwrap_unique(),
            u.unwrap_unique(),
            v.unwrap_unique(),
            work.unwrap_unique(),
        )
    })
}

/// This function factorizes the M-by-N matrix A into bidiagonal form U B V^T. The diagonal and superdiagonal of the matrix B are stored in
/// the diagonal and superdiagonal of A. The orthogonal matrices U and V are stored as compressed Householder vectors in the remaining elements
/// of A. The Householder coefficients are stored in the vectors tau_U and tau_V. The length of tau_U must equal the number of elements in
/// the diagonal of A and the length of tau_V should be one element shorter.
#[doc(alias = "gsl_linalg_bidiag_decomp")]
pub fn bidiag_decomp(
    a: &mut ::MatrixF64,
    tau_u: &mut ::VectorF64,
    tau_v: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_bidiag_decomp(
            a.unwrap_unique(),
            tau_u.unwrap_unique(),
            tau_v.unwrap_unique(),
        )
    })
}

/// This function unpacks the bidiagonal decomposition of A produced by gsl_linalg_bidiag_decomp, (A, tau_U, tau_V) into the separate orthogonal
/// matrices U, V and the diagonal vector diag and superdiagonal superdiag. Note that U is stored as a compact M-by-N orthogonal matrix satisfying
/// U^T U = I for efficiency.
#[doc(alias = "gsl_linalg_bidiag_unpack")]
pub fn bidiag_unpack(
    a: &mut ::MatrixF64,
    tau_u: &::VectorF64,
    u: &mut ::MatrixF64,
    tau_v: &::VectorF64,
    v: &mut ::MatrixF64,
    diag: &mut ::VectorF64,
    superdiag: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_bidiag_unpack(
            a.unwrap_unique(),
            tau_u.unwrap_shared(),
            u.unwrap_unique(),
            tau_v.unwrap_shared(),
            v.unwrap_unique(),
            diag.unwrap_unique(),
            superdiag.unwrap_unique(),
        )
    })
}

/// This function unpacks the bidiagonal decomposition of A produced by gsl_linalg_bidiag_decomp, (A, tau_U, tau_V) into the separate orthogonal
/// matrices U, V and the diagonal vector diag and superdiagonal superdiag. The matrix U is stored in-place in A.
#[doc(alias = "gsl_linalg_bidiag_unpack2")]
pub fn bidiag_unpack2(
    a: &mut ::MatrixF64,
    tau_u: &mut ::VectorF64,
    tau_v: &mut ::VectorF64,
    v: &mut ::MatrixF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_bidiag_unpack2(
            a.unwrap_unique(),
            tau_u.unwrap_unique(),
            tau_v.unwrap_unique(),
            v.unwrap_unique(),
        )
    })
}

/// This function unpacks the diagonal and superdiagonal of the bidiagonal decomposition of A from gsl_linalg_bidiag_decomp, into the diagonal
/// vector diag and superdiagonal vector superdiag.
#[doc(alias = "gsl_linalg_bidiag_unpack_B")]
pub fn bidiag_unpack_B(
    a: &::MatrixF64,
    diag: &mut ::VectorF64,
    superdiag: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_bidiag_unpack_B(
            a.unwrap_shared(),
            diag.unwrap_unique(),
            superdiag.unwrap_unique(),
        )
    })
}

/// This function prepares a Householder transformation P = I - \tau v v^T which can be used to zero all the elements of the input vector except
/// the first. On output the transformation is stored in the vector v and the scalar \tau is returned.
#[doc(alias = "gsl_linalg_householder_transform")]
pub fn householder_transform(v: &mut ::VectorF64) -> f64 {
    unsafe { sys::gsl_linalg_householder_transform(v.unwrap_unique()) }
}

/// This function prepares a Householder transformation P = I - \tau v v^T which can be used to zero all the elements of the input vector except
/// the first. On output the transformation is stored in the vector v and the scalar \tau is returned.
#[doc(alias = "gsl_linalg_complex_householder_transform")]
pub fn complex_householder_transform(v: &mut ::VectorComplexF64) -> ::ComplexF64 {
    unsafe {
        ::std::mem::transmute(sys::gsl_linalg_complex_householder_transform(
            v.unwrap_unique(),
        ))
    }
}

/// This function applies the Householder matrix P defined by the scalar tau and the vector v to the left-hand side of the matrix A. On output
/// the result P A is stored in A.
#[doc(alias = "gsl_linalg_householder_hm")]
pub fn householder_hm(tau: f64, v: &::VectorF64, a: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_householder_hm(tau, v.unwrap_shared(), a.unwrap_unique())
    })
}

/// This function applies the Householder matrix P defined by the scalar tau and the vector v to the left-hand side of the matrix A. On output
/// the result P A is stored in A.
#[doc(alias = "gsl_linalg_complex_householder_hm")]
pub fn complex_householder_hm(
    tau: &::ComplexF64,
    v: &::VectorComplexF64,
    a: &mut ::MatrixComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_householder_hm(
            ::std::mem::transmute(*tau),
            v.unwrap_shared(),
            a.unwrap_unique(),
        )
    })
}

/// This function applies the Householder matrix P defined by the scalar tau and the vector v to the right-hand side of the matrix A. On output
/// the result A P is stored in A.
#[doc(alias = "gsl_linalg_householder_mh")]
pub fn householder_mh(tau: f64, v: &::VectorF64, a: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_householder_mh(tau, v.unwrap_shared(), a.unwrap_unique())
    })
}

/// This function applies the Householder matrix P defined by the scalar tau and the vector v to the right-hand side of the matrix A. On output
/// the result A P is stored in A.
#[doc(alias = "gsl_linalg_complex_householder_mh")]
pub fn complex_householder_mh(
    tau: &::ComplexF64,
    v: &::VectorComplexF64,
    a: &mut ::MatrixComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_householder_mh(
            ::std::mem::transmute(*tau),
            v.unwrap_shared(),
            a.unwrap_unique(),
        )
    })
}

/// This function applies the Householder transformation P defined by the scalar tau and the vector v to the vector w. On output the result P
/// w is stored in w.
#[doc(alias = "gsl_linalg_householder_hv")]
pub fn householder_hv(tau: f64, v: &::VectorF64, w: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_householder_hv(tau, v.unwrap_shared(), w.unwrap_unique())
    })
}

/// This function applies the Householder transformation P defined by the scalar tau and the vector v to the vector w. On output the result P
/// w is stored in w.
#[doc(alias = "gsl_linalg_complex_householder_hv")]
pub fn complex_householder_hv(
    tau: &::ComplexF64,
    v: &::VectorComplexF64,
    w: &mut ::VectorComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_householder_hv(
            ::std::mem::transmute(*tau),
            v.unwrap_shared(),
            w.unwrap_unique(),
        )
    })
}

/// This function solves the system A x = b directly using Householder transformations. On output the solution is stored in x and b is not
/// modified. The matrix A is destroyed by the Householder transformations.
#[doc(alias = "gsl_linalg_HH_solve")]
pub fn HH_solve(mut a: ::MatrixF64, b: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_HH_solve(a.unwrap_unique(), b.unwrap_shared(), x.unwrap_unique())
    })
}

/// This function solves the system A x = b in-place using Householder transformations. On input x should contain the right-hand side b,
/// which is replaced by the solution on output. The matrix A is destroyed by the Householder transformations.
#[doc(alias = "gsl_linalg_HH_svx")]
pub fn HH_svx(mut a: ::MatrixF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_HH_svx(a.unwrap_unique(), x.unwrap_unique()) })
}

/// This function solves the general N-by-N system A x = b where A is tridiagonal (N >= 2). The super-diagonal and sub-diagonal vectors
/// e and f must be one element shorter than the diagonal vector diag. The form of A for the 4-by-4 case is shown below,
///
/// ```text
/// A = ( d_0 e_0  0   0  )
///     ( f_0 d_1 e_1  0  )
///     (  0  f_1 d_2 e_2 )
///     (  0   0  f_2 d_3 )
/// ```
#[doc(alias = "gsl_linalg_solve_tridiag")]
pub fn solve_tridiag(
    diag: &::VectorF64,
    e: &::VectorF64,
    f: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_solve_tridiag(
            diag.unwrap_shared(),
            e.unwrap_shared(),
            f.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the general N-by-N system A x = b where A is symmetric tridiagonal (N >= 2). The off-diagonal vector e must be one
/// element shorter than the diagonal vector diag. The form of A for the 4-by-4 case is shown below,
///
/// ```text
/// A = ( d_0 e_0  0   0  )
///     ( e_0 d_1 e_1  0  )
///     (  0  e_1 d_2 e_2 )
///     (  0   0  e_2 d_3 )
/// ```
#[doc(alias = "gsl_linalg_solve_symm_tridiag")]
pub fn solve_symm_tridiag(
    diag: &::VectorF64,
    e: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_solve_symm_tridiag(
            diag.unwrap_shared(),
            e.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the general N-by-N system A x = b where A is cyclic tridiagonal (N >= 3). The cyclic super-diagonal and sub-diagonal
/// vectors e and f must have the same number of elements as the diagonal vector diag. The form of A for the 4-by-4 case is shown below,
///
/// ```text
/// A = ( d_0 e_0  0  f_3 )
///     ( f_0 d_1 e_1  0  )
///     (  0  f_1 d_2 e_2 )
///     ( e_3  0  f_2 d_3 )
/// ```
#[doc(alias = "gsl_linalg_solve_cyc_tridiag")]
pub fn solve_cyc_tridiag(
    diag: &::VectorF64,
    e: &::VectorF64,
    f: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_solve_cyc_tridiag(
            diag.unwrap_shared(),
            e.unwrap_shared(),
            f.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function solves the general N-by-N system A x = b where A is symmetric cyclic tridiagonal (N >= 3). The cyclic off-diagonal vector
/// e must have the same number of elements as the diagonal vector diag. The form of A for the 4-by-4 case is shown below,
///
/// ```text
/// A = ( d_0 e_0  0  e_3 )
///     ( e_0 d_1 e_1  0  )
///     (  0  e_1 d_2 e_2 )
///     ( e_3  0  e_2 d_3 )
/// ```
#[doc(alias = "gsl_linalg_solve_symm_cyc_tridiag")]
pub fn solve_symm_cyc_tridiag(
    diag: &::VectorF64,
    e: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_solve_symm_cyc_tridiag(
            diag.unwrap_shared(),
            e.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

/// This function replaces the matrix A with its balanced counterpart and stores the diagonal elements of the similarity transformation into
/// the vector D.
#[doc(alias = "gsl_linalg_balance_matrix")]
pub fn balance_matrix(a: &mut ::MatrixF64, d: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_balance_matrix(a.unwrap_unique(), d.unwrap_unique()) })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_decomp")]
pub fn pcholesky_decomp(a: &mut ::MatrixF64, p: &mut ::Permutation) -> Value {
    Value::from(unsafe { sys::gsl_linalg_pcholesky_decomp(a.unwrap_unique(), p.unwrap_unique()) })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_solve")]
pub fn pcholesky_solve(
    LDLT: &::MatrixF64,
    p: &::Permutation,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_pcholesky_solve(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_svx")]
pub fn pcholesky_svx(LDLT: &::MatrixF64, p: &::Permutation, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_pcholesky_svx(LDLT.unwrap_shared(), p.unwrap_shared(), x.unwrap_unique())
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_decomp2")]
pub fn pcholesky_decomp2(A: &mut ::MatrixF64, p: &mut ::Permutation, S: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_pcholesky_decomp2(A.unwrap_unique(), p.unwrap_unique(), S.unwrap_unique())
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_solve2")]
pub fn pcholesky_solve2(
    LDLT: &::MatrixF64,
    p: &::Permutation,
    S: &::VectorF64,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_pcholesky_solve2(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            S.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_svx2")]
pub fn pcholesky_svx2(
    LDLT: &::MatrixF64,
    p: &::Permutation,
    S: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_pcholesky_svx2(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            S.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_invert")]
pub fn pcholesky_invert(LDLT: &::MatrixF64, p: &::Permutation, Ainv: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_pcholesky_invert(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            Ainv.unwrap_unique(),
        )
    })
}

/// Returns `(Value, rcond)`.
#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_pcholesky_rcond")]
pub fn pcholesky_rcond(
    LDLT: &::MatrixF64,
    p: &::Permutation,
    work: &mut ::VectorF64,
) -> (Value, f64) {
    let mut rcond = 0.;
    let ret = unsafe {
        sys::gsl_linalg_pcholesky_rcond(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            &mut rcond,
            work.unwrap_unique(),
        )
    };
    (Value::from(ret), rcond)
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_mcholesky_decomp")]
pub fn mcholesky_decomp(A: &mut ::MatrixF64, p: &mut ::Permutation, E: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_mcholesky_decomp(A.unwrap_unique(), p.unwrap_unique(), E.unwrap_unique())
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_mcholesky_solve")]
pub fn mcholesky_solve(
    LDLT: &::MatrixF64,
    p: &::Permutation,
    b: &::VectorF64,
    x: &mut ::VectorF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_mcholesky_solve(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_mcholesky_svx")]
pub fn mcholesky_svx(LDLT: &::MatrixF64, p: &::Permutation, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_mcholesky_svx(LDLT.unwrap_shared(), p.unwrap_shared(), x.unwrap_unique())
    })
}

/// Returns `(Value, rcond)`.
#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_mcholesky_rcond")]
pub fn mcholesky_rcond(
    LDLT: &::MatrixF64,
    p: &::Permutation,
    work: &mut ::VectorF64,
) -> (Value, f64) {
    let mut rcond = 0.;
    let ret = unsafe {
        sys::gsl_linalg_mcholesky_rcond(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            &mut rcond,
            work.unwrap_unique(),
        )
    };
    (Value::from(ret), rcond)
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_mcholesky_invert")]
pub fn mcholesky_invert(LDLT: &::MatrixF64, p: &::Permutation, Ainv: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_mcholesky_invert(
            LDLT.unwrap_shared(),
            p.unwrap_shared(),
            Ainv.unwrap_unique(),
        )
    })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_cholesky_band_decomp")]
pub fn cholesky_band_decomp(A: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_cholesky_band_decomp(A.unwrap_unique()) })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_cholesky_band_solve")]
pub fn cholesky_band_solve(LLT: &::MatrixF64, b: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_band_solve(
            LLT.unwrap_shared(),
            b.unwrap_shared(),
            x.unwrap_unique(),
        )
    })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_cholesky_band_svx")]
pub fn cholesky_band_svx(LLT: &::MatrixF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_band_svx(LLT.unwrap_shared(), x.unwrap_unique())
    })
}

#[cfg(feature = "v2_7")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_7")))]
#[doc(alias = "gsl_linalg_cholesky_band_solvem")]
pub fn cholesky_band_solvem(LLT: &::MatrixF64, B: &::MatrixF64, X: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_band_solvem(
            LLT.unwrap_shared(),
            B.unwrap_shared(),
            X.unwrap_unique(),
        )
    })
}

#[cfg(feature = "v2_7")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_7")))]
#[doc(alias = "gsl_linalg_cholesky_band_svxm")]
pub fn cholesky_band_svxm(LLT: &::MatrixF64, X: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_band_svxm(LLT.unwrap_shared(), X.unwrap_unique())
    })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_cholesky_band_invert")]
pub fn cholesky_band_invert(LLT: &::MatrixF64, Ainv: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_band_invert(LLT.unwrap_shared(), Ainv.unwrap_unique())
    })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_cholesky_band_unpack")]
pub fn cholesky_band_unpack(LLT: &::MatrixF64, L: &mut ::MatrixF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_cholesky_band_unpack(LLT.unwrap_shared(), L.unwrap_unique())
    })
}

/// Returns `(Value, rcond)`.
#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_cholesky_band_rcond")]
pub fn cholesky_band_rcond(LLT: &::MatrixF64, work: &mut ::VectorF64) -> (Value, f64) {
    let mut rcond = 0.;
    let ret = unsafe {
        sys::gsl_linalg_cholesky_band_rcond(LLT.unwrap_shared(), &mut rcond, work.unwrap_unique())
    };
    (Value::from(ret), rcond)
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_decomp")]
pub fn ldlt_decomp(A: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_ldlt_decomp(A.unwrap_unique()) })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_solve")]
pub fn ldlt_solve(LDLT: &::MatrixF64, b: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_ldlt_solve(LDLT.unwrap_shared(), b.unwrap_shared(), x.unwrap_unique())
    })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_svx")]
pub fn ldlt_svx(LDLT: &::MatrixF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_ldlt_svx(LDLT.unwrap_shared(), x.unwrap_unique()) })
}

/// Returns `(Value, rcond)`.
#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_rcond")]
pub fn ldlt_rcond(LDLT: &::MatrixF64, work: &mut ::VectorF64) -> (Value, f64) {
    let mut rcond = 0.;
    let ret = unsafe {
        sys::gsl_linalg_ldlt_rcond(LDLT.unwrap_shared(), &mut rcond, work.unwrap_unique())
    };
    (Value::from(ret), rcond)
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_band_decomp")]
pub fn ldlt_band_decomp(A: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_ldlt_band_decomp(A.unwrap_unique()) })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_band_solve")]
pub fn ldlt_band_solve(LDLT: &::MatrixF64, b: &::VectorF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_ldlt_band_solve(LDLT.unwrap_shared(), b.unwrap_shared(), x.unwrap_unique())
    })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_band_svx")]
pub fn ldlt_band_svx(LDLT: &::MatrixF64, x: &mut ::VectorF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_ldlt_band_svx(LDLT.unwrap_shared(), x.unwrap_unique()) })
}

#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_band_unpack")]
pub fn ldlt_band_unpack(LDLT: &::MatrixF64, L: &mut ::MatrixF64, D: &mut ::VectorF64) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_ldlt_band_unpack(LDLT.unwrap_shared(), L.unwrap_unique(), D.unwrap_unique())
    })
}

/// Returns `(Value, rcond)`.
#[cfg(feature = "v2_6")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_6")))]
#[doc(alias = "gsl_linalg_ldlt_band_rcond")]
pub fn ldlt_band_rcond(LDLT: &::MatrixF64, work: &mut ::VectorF64) -> (Value, f64) {
    let mut rcond = 0.;
    let ret = unsafe {
        sys::gsl_linalg_ldlt_band_rcond(LDLT.unwrap_shared(), &mut rcond, work.unwrap_unique())
    };
    (Value::from(ret), rcond)
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_tri_upper_invert")]
pub fn tri_upper_invert(T: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_tri_upper_invert(T.unwrap_unique()) })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_tri_lower_invert")]
pub fn tri_lower_invert(T: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_tri_lower_invert(T.unwrap_unique()) })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_tri_upper_unit_invert")]
pub fn tri_upper_unit_invert(T: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_tri_upper_unit_invert(T.unwrap_unique()) })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_tri_lower_unit_invert")]
pub fn tri_lower_unit_invert(T: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_tri_lower_unit_invert(T.unwrap_unique()) })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_complex_tri_invert")]
pub fn tri_invert(
    Uplo: enums::CblasUplo,
    Diag: enums::CblasDiag,
    T: &mut ::MatrixComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_tri_invert(Uplo.into(), Diag.into(), T.unwrap_unique())
    })
}

#[doc(alias = "gsl_linalg_complex_tri_invert")]
pub fn complex_tri_invert(
    Uplo: enums::CblasUplo,
    Diag: enums::CblasDiag,
    T: &mut ::MatrixComplexF64,
) -> Value {
    Value::from(unsafe {
        sys::gsl_linalg_complex_tri_invert(Uplo.into(), Diag.into(), T.unwrap_unique())
    })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_tri_LTL")]
pub fn tri_LTL(L: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_tri_LTL(L.unwrap_unique()) })
}

#[cfg(feature = "v2_2")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
#[doc(alias = "gsl_linalg_tri_UL")]
pub fn tri_UL(LU: &mut ::MatrixF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_tri_UL(LU.unwrap_unique()) })
}

#[doc(alias = "gsl_linalg_complex_tri_LHL")]
pub fn complex_tri_LHL(L: &mut ::MatrixComplexF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_complex_tri_LHL(L.unwrap_unique()) })
}

#[doc(alias = "gsl_linalg_complex_tri_UL")]
pub fn complex_tri_UL(LU: &mut ::MatrixComplexF64) -> Value {
    Value::from(unsafe { sys::gsl_linalg_complex_tri_UL(LU.unwrap_unique()) })
}

/// Returns `(c, s)`.
#[doc(alias = "gsl_linalg_givens")]
pub fn givens(a: f64, b: f64) -> (f64, f64) {
    let mut c = 0.;
    let mut s = 0.;
    unsafe { sys::gsl_linalg_givens(a, b, &mut c, &mut s) };
    (c, s)
}

#[doc(alias = "gsl_linalg_givens_gv")]
pub fn givens_gv(v: &mut ::VectorF64, i: usize, j: usize, c: f64, s: f64) {
    unsafe { sys::gsl_linalg_givens_gv(v.unwrap_unique(), i, j, c, s) }
}
