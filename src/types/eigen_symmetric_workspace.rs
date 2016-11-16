//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Real Symmetric Matrices

For real symmetric matrices, the library uses the symmetric bidiagonalization and QR reduction
method. This is described in Golub & van Loan, section 8.3. The computed eigenvalues are accurate to
an absolute accuracy of \epsilon ||A||_2, where \epsilon is the machine precision.

#Complex Hermitian Matrices

For hermitian matrices, the library uses the complex form of the symmetric bidiagonalization and QR
reduction method.

#Real Nonsymmetric Matrices

The solution of the real nonsymmetric eigensystem problem for a matrix A involves computing the
Schur decomposition

A = Z T Z^T

where Z is an orthogonal matrix of Schur vectors and T, the Schur form, is quasi upper triangular
with diagonal 1-by-1 blocks which are real eigenvalues of A, and diagonal 2-by-2 blocks whose
eigenvalues are complex conjugate eigenvalues of A. The algorithm used is the double-shift Francis
method.

#Real Generalized Symmetric-Definite Eigensystems

The real generalized symmetric-definite eigenvalue problem is to find eigenvalues \lambda and
eigenvectors x such that

A x = lambda B x

where A and B are symmetric matrices, and B is positive-definite. This problem reduces to the
standard symmetric eigenvalue problem by applying the Cholesky decomposition to B:

```latex
                      A x = lambda B x
                      A x = lambda L L^t x
( L^{-1} A L^{-t} ) L^t x = lambda L^t x
```

Therefore, the problem becomes C y = lambda y where C = L^{-1} A L^{-t} is symmetric, and y = L^t x.
The standard symmetric eigensolver can be applied to the matrix C. The resulting eigenvectors are
backtransformed to find the vectors of the original problem. The eigenvalues and eigenvectors of the
generalized symmetric-definite eigenproblem are always real.

#Complex Generalized Hermitian-Definite Eigensystems

The complex generalized hermitian-definite eigenvalue problem is to find eigenvalues \lambda and
eigenvectors x such that

A x = \lambda B x

where A and B are hermitian matrices, and B is positive-definite. Similarly to the real case, this
can be reduced to C y = \lambda y where C = L^{-1} A L^{-H} is hermitian, and y = L^H x. The
standard hermitian eigensolver can be applied to the matrix C. The resulting eigenvectors are
backtransformed to find the vectors of the original problem. The eigenvalues of the generalized
hermitian-definite eigenproblem are always real.

#Real Generalized Nonsymmetric Eigensystems

Given two square matrices (A, B), the generalized nonsymmetric eigenvalue problem is to find
eigenvalues \lambda and eigenvectors x such that

A x = \lambda B x

We may also define the problem as finding eigenvalues \mu and eigenvectors y such that

\mu A y = B y
Note that these two problems are equivalent (with \lambda = 1/\mu) if neither \lambda nor \mu is
zero. If say, \lambda is zero, then it is still a well defined eigenproblem, but its alternate
problem involving \mu is not. Therefore, to allow for zero (and infinite) eigenvalues, the problem
which is actually solved is

\beta A x = \alpha B x
The eigensolver routines below will return two values \alpha and \beta and leave it to the user to
perform the divisions \lambda = \alpha / \beta and \mu = \beta / \alpha.

If the determinant of the matrix pencil A - \lambda B is zero for all \lambda, the problem is said
to be singular; otherwise it is called regular. Singularity normally leads to some
\alpha = \beta = 0 which means the eigenproblem is ill-conditioned and generally does not have well
defined eigenvalue solutions. The routines below are intended for regular matrix pencils and could
yield unpredictable results when applied to singular pencils.

The solution of the real generalized nonsymmetric eigensystem problem for a matrix pair (A, B)
involves computing the generalized Schur decomposition

A = Q S Z^T
B = Q T Z^T
where Q and Z are orthogonal matrices of left and right Schur vectors respectively, and (S, T) is
the generalized Schur form whose diagonal elements give the \alpha and \beta values. The algorithm
used is the QZ method due to Moler and Stewart (see references).
!*/

use ffi;
use enums;
use types::{MatrixF64, MatrixComplexF64, VectorF64, VectorComplexF64};

pub struct EigenSymmetricWorkspace {
    w: *mut ffi::gsl_eigen_symm_workspace
}

impl EigenSymmetricWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n real symmetric
    /// matrices. The size of the workspace is O(2n).
    pub fn new(n: usize) -> Option<EigenSymmetricWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_symm_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenSymmetricWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues of the real symmetric matrix `A`. The diagonal and
    /// lower triangular part of `A` are destroyed during the computation, but the strict upper
    /// triangular part is not referenced. The eigenvalues are stored in the vector `eval` and are
    /// unordered.
    pub fn symm(&self, A: &MatrixF64, eval: &mut VectorF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_symm(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval), self.w) }
    }
}

impl Drop for EigenSymmetricWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_symm_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_symm_workspace> for EigenSymmetricWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_symm_workspace) -> EigenSymmetricWorkspace {
        EigenSymmetricWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenSymmetricWorkspace) -> *mut ffi::gsl_eigen_symm_workspace {
        t.w
    }
}

pub struct EigenSymmetricVWorkspace {
    w: *mut ffi::gsl_eigen_symmv_workspace
}

impl EigenSymmetricVWorkspace {
    /// This function allocates a workspace for computing eigenvalues and eigenvectors of n-by-n
    /// real symmetric matrices. The size of the workspace is O(4n).
    pub fn new(n: usize) -> Option<EigenSymmetricVWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_symmv_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenSymmetricVWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues and eigenvectors of the real symmetric matrix `A`.
    /// The diagonal and lower triangular part of `A` are destroyed during the computation, but the
    /// strict upper triangular part is not referenced. The eigenvalues are stored in the vector
    /// `eval` and are unordered. The corresponding eigenvectors are stored in the columns of the
    /// matrix `evec`. For example, the eigenvector in the first column corresponds to the first
    /// eigenvalue. The eigenvectors are guaranteed to be mutually orthogonal and normalised to unit
    /// magnitude.
    pub fn symmv(&self, A: &MatrixF64, eval: &mut VectorF64, evec: &mut MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_symmv(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval),
                                      ffi::FFI::unwrap(evec), self.w) }
    }
}

impl Drop for EigenSymmetricVWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_symmv_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_symmv_workspace> for EigenSymmetricVWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_symmv_workspace) -> EigenSymmetricVWorkspace {
        EigenSymmetricVWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenSymmetricVWorkspace) -> *mut ffi::gsl_eigen_symmv_workspace {
        t.w
    }
}

pub struct EigenHermitianWorkspace {
    w: *mut ffi::gsl_eigen_herm_workspace,
}

impl EigenHermitianWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n complex hermitian
    /// matrices. The size of the workspace is O(3n).
    pub fn new(n: usize) -> Option<EigenHermitianWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_herm_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenHermitianWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues of the complex hermitian matrix `A`. Additional
    /// workspace of the appropriate size must be provided in `self`. The diagonal and lower
    /// triangular part of `A` are destroyed during the computation, but the strict upper triangular
    /// part is not referenced. The imaginary parts of the diagonal are assumed to be zero and are
    /// not referenced. The eigenvalues are stored in the vector `eval` and are unordered.
    pub fn herm(&self, A: &MatrixComplexF64, eval: &mut VectorF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_herm(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval), self.w) }
    }
}

impl Drop for EigenHermitianWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_herm_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_herm_workspace> for EigenHermitianWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_herm_workspace) -> EigenHermitianWorkspace {
        EigenHermitianWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenHermitianWorkspace) -> *mut ffi::gsl_eigen_herm_workspace {
        t.w
    }
}

pub struct EigenHermitianVWorkspace {
    w: *mut ffi::gsl_eigen_hermv_workspace
}

impl EigenHermitianVWorkspace {
    /// This function allocates a workspace for computing eigenvalues and eigenvectors of n-by-n
    /// complex hermitian matrices. The size of the workspace is O(5n).
    pub fn new(n: usize) -> Option<EigenHermitianVWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_hermv_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenHermitianVWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues and eigenvectors of the complex hermitian matrix `A`.
    /// Additional workspace of the appropriate size must be provided in `self`. The diagonal and
    /// lower triangular part of `A` are destroyed during the computation, but the strict upper
    /// triangular part is not referenced. The imaginary parts of the diagonal are assumed to be
    /// zero and are not referenced. The eigenvalues are stored in the vector `eval` and are
    /// unordered. The corresponding complex eigenvectors are stored in the columns of the matrix
    /// `evec`. For example, the eigenvector in the first column corresponds to the first
    /// eigenvalue. The eigenvectors are guaranteed to be mutually orthogonal and normalised to unit
    /// magnitude.
    pub fn hermv(&self, A: &MatrixComplexF64, eval: &mut VectorF64,
                 evec: &mut MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_hermv(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval),
                                      ffi::FFI::unwrap(evec), self.w) }
    }
}

impl Drop for EigenHermitianVWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_hermv_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_hermv_workspace> for EigenHermitianVWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_hermv_workspace) -> EigenHermitianVWorkspace {
        EigenHermitianVWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenHermitianVWorkspace) -> *mut ffi::gsl_eigen_hermv_workspace {
        t.w
    }
}

pub struct EigenNonSymmWorkspace {
    w: *mut ffi::gsl_eigen_nonsymm_workspace
}

impl EigenNonSymmWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n complex hermitian
    /// matrices. The size of the workspace is O(3n).
    pub fn new(n: usize) -> Option<EigenNonSymmWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_nonsymm_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenNonSymmWorkspace {
                w: tmp,
            })
        }
    }

    /// This function sets some parameters which determine how the eigenvalue problem is solved in
    /// subsequent calls to gsl_eigen_nonsymm.
    ///
    /// If compute_t is set to 1, the full Schur form T will be computed by gsl_eigen_nonsymm. If it
    /// is set to 0, T will not be computed (this is the default setting). Computing the full Schur
    /// form T requires approximately 1.5–2 times the number of flops.
    ///
    /// If balance is set to 1, a balancing transformation is applied to the matrix prior to
    /// computing eigenvalues. This transformation is designed to make the rows and columns of the
    /// matrix have comparable norms, and can result in more accurate eigenvalues for matrices whose
    /// entries vary widely in magnitude. See
    /// [`Balancing`](http://www.gnu.org/software/gsl/manual/html_node/Balancing.html#Balancing) for
    /// more information. Note that the balancing transformation does not preserve the orthogonality
    /// of the Schur vectors, so if you wish to compute the Schur vectors with gsl_eigen_nonsymm_Z
    /// you will obtain the Schur vectors of the balanced matrix instead of the original matrix.
    /// The relationship will be
    ///
    /// T = Q^t D^(-1) A D Q
    ///
    /// where Q is the matrix of Schur vectors for the balanced matrix, and D is the balancing
    /// transformation. Then gsl_eigen_nonsymm_Z will compute a matrix Z which satisfies
    ///
    /// T = Z^(-1) A Z
    ///
    /// with Z = D Q. Note that Z will not be orthogonal. For this reason, balancing is not
    /// performed by default.
    pub fn params(&self, compute_t: i32, balance: i32) {
        unsafe { ffi::gsl_eigen_nonsymm_params(compute_t, balance, self.w) }
    }

    /// This function computes the eigenvalues of the real nonsymmetric matrix `A` and stores them
    /// in the vector `eval`. If T is desired, it is stored in the upper portion of `A` on output.
    /// Otherwise, on output, the diagonal of `A` will contain the 1-by-1 real eigenvalues and
    /// 2-by-2 complex conjugate eigenvalue systems, and the rest of `A` is destroyed. In rare
    /// cases, this function may fail to find all eigenvalues. If this happens, an error code is
    /// returned and the number of converged eigenvalues is stored in w->n_evals. The converged
    /// eigenvalues are stored in the beginning of `eval`.
    pub fn nonsymm(&self, A: &MatrixF64, eval: &mut VectorComplexF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_nonsymm(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval), self.w) }
    }

    /// This function is identical to gsl_eigen_nonsymm except that it also computes the Schur
    /// vectors and stores them into `Z`.
    pub fn nonsymm_Z(&self, A: &MatrixF64, eval: &mut VectorComplexF64,
                     Z: &mut MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_nonsymm_Z(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval),
                                          ffi::FFI::unwrap(Z), self.w) }
    }

    pub fn n_evals(&self) -> usize {
        unsafe { (*self.w).n_evals }
    }
}

impl Drop for EigenNonSymmWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_nonsymm_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_nonsymm_workspace> for EigenNonSymmWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_nonsymm_workspace) -> EigenNonSymmWorkspace {
        EigenNonSymmWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenNonSymmWorkspace) -> *mut ffi::gsl_eigen_nonsymm_workspace {
        t.w
    }
}

pub struct EigenNonSymmVWorkspace {
    w: *mut ffi::gsl_eigen_nonsymmv_workspace
}

impl EigenNonSymmVWorkspace {
    /// This function allocates a workspace for computing eigenvalues and eigenvectors of n-by-n
    /// real nonsymmetric matrices. The size of the workspace is O(5n).
    pub fn new(n: usize) -> Option<EigenNonSymmVWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_nonsymmv_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenNonSymmVWorkspace {
                w: tmp,
            })
        }
    }

    /// This function sets parameters which determine how the eigenvalue problem is solved in
    /// subsequent calls to gsl_eigen_nonsymmv. If `balance` is set to 1, a balancing transformation
    /// is applied to the matrix. See gsl_eigen_nonsymm_params for more information. Balancing is
    /// turned off by default since it does not preserve the orthogonality of the Schur vectors.
    pub fn params(&self, balance: i32) {
        unsafe { ffi::gsl_eigen_nonsymmv_params(balance, self.w) }
    }

    /// This function computes eigenvalues and right eigenvectors of the n-by-n real nonsymmetric
    /// matrix `A`. It first calls gsl_eigen_nonsymm to compute the eigenvalues, Schur form T, and
    /// Schur vectors. Then it finds eigenvectors of T and backtransforms them using the Schur
    /// vectors. The Schur vectors are destroyed in the process, but can be saved by using
    /// gsl_eigen_nonsymmv_Z. The computed eigenvectors are normalized to have unit magnitude. On
    /// output, the upper portion of `A` contains the Schur form T. If gsl_eigen_nonsymm fails, no
    /// eigenvectors are computed, and an error code is returned.
    pub fn nonsymmv(&self, A: &MatrixF64, eval: &mut VectorComplexF64,
                    evec: &mut MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_nonsymmv(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval),
                                         ffi::FFI::unwrap(evec), self.w) }
    }

    /// This function is identical to gsl_eigen_nonsymmv except that it also saves the Schur vectors
    /// into `Z`.
    pub fn nonsymmv_Z(&self, A: &MatrixF64, eval: &mut VectorComplexF64,
                      evec: &mut MatrixComplexF64, Z: &mut MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_nonsymmv_Z(ffi::FFI::unwrap(A), ffi::FFI::unwrap(eval),
                                           ffi::FFI::unwrap(evec), ffi::FFI::unwrap(Z), self.w) }
    }
}

impl Drop for EigenNonSymmVWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_nonsymmv_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_nonsymmv_workspace> for EigenNonSymmVWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_nonsymmv_workspace) -> EigenNonSymmVWorkspace {
        EigenNonSymmVWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenNonSymmVWorkspace) -> *mut ffi::gsl_eigen_nonsymmv_workspace {
        t.w
    }
}

pub struct EigenGenSymmWorkspace {
    w: *mut ffi::gsl_eigen_gensymm_workspace
}

impl EigenGenSymmWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n real generalized
    /// symmetric-definite eigensystems. The size of the workspace is O(2n).
    pub fn new(n: usize) -> Option<EigenGenSymmWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_gensymm_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenGenSymmWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues of the real generalized symmetric-definite matrix
    /// pair (A, B), and stores them in `eval`, using the method outlined above. On output, `B`
    /// contains its Cholesky decomposition and `A` is destroyed.
    pub fn gensymm(&self, A: MatrixF64, B: &MatrixF64, eval: &mut VectorF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_gensymm(ffi::FFI::unwrap(&A), ffi::FFI::unwrap(B),
                                        ffi::FFI::unwrap(eval), self.w) }
    }
}

impl Drop for EigenGenSymmWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_gensymm_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_gensymm_workspace> for EigenGenSymmWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_gensymm_workspace) -> EigenGenSymmWorkspace {
        EigenGenSymmWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenGenSymmWorkspace) -> *mut ffi::gsl_eigen_gensymm_workspace {
        t.w
    }
}

pub struct EigenGenSymmVWorkspace {
    w: *mut ffi::gsl_eigen_gensymmv_workspace
}

impl EigenGenSymmVWorkspace {
    /// This function allocates a workspace for computing eigenvalues and eigenvectors of n-by-n
    /// real generalized symmetric-definite eigensystems. The size of the workspace is O(4n).
    pub fn new(n: usize) -> Option<EigenGenSymmVWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_gensymmv_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenGenSymmVWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues and eigenvectors of the real generalized
    /// symmetric-definite matrix pair (A, B), and stores them in `eval` and `evec` respectively.
    /// The computed eigenvectors are normalized to have unit magnitude. On output, `B` contains its
    /// Cholesky decomposition and `A` is destroyed.
    pub fn gensymmv(&self, A: MatrixF64, B: &MatrixF64, eval: &mut VectorF64,
                    evec: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_gensymmv(ffi::FFI::unwrap(&A), ffi::FFI::unwrap(B),
                                         ffi::FFI::unwrap(eval), ffi::FFI::unwrap(evec), self.w) }
    }
}

impl Drop for EigenGenSymmVWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_gensymmv_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_gensymmv_workspace> for EigenGenSymmVWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_gensymmv_workspace) -> EigenGenSymmVWorkspace {
        EigenGenSymmVWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenGenSymmVWorkspace) -> *mut ffi::gsl_eigen_gensymmv_workspace {
        t.w
    }
}

pub struct EigenGenHermWorkspace {
    w: *mut ffi::gsl_eigen_genherm_workspace
}

impl EigenGenHermWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n complex generalized
    /// hermitian-definite eigensystems. The size of the workspace is O(3n).
    pub fn new(n: usize) -> Option<EigenGenHermWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_genherm_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenGenHermWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues of the complex generalized hermitian-definite matrix
    /// pair (A, B), and stores them in `eval`, using the method outlined above. On output, `B`
    /// contains its Cholesky decomposition and `A` is destroyed.
    pub fn genherm(&self, A: MatrixComplexF64, B: &MatrixComplexF64,
                   eval: &mut VectorF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_genherm(ffi::FFI::unwrap(&A), ffi::FFI::unwrap(B),
                                        ffi::FFI::unwrap(eval), self.w) }
    }
}

impl Drop for EigenGenHermWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_genherm_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_genherm_workspace> for EigenGenHermWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_genherm_workspace) -> EigenGenHermWorkspace {
        EigenGenHermWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenGenHermWorkspace) -> *mut ffi::gsl_eigen_genherm_workspace {
        t.w
    }
}

pub struct EigenGenHermVWorkspace {
    w: *mut ffi::gsl_eigen_genhermv_workspace
}

impl EigenGenHermVWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n complex generalized
    /// hermitian-definite eigensystems. The size of the workspace is O(3n).
    pub fn new(n: usize) -> Option<EigenGenHermVWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_genhermv_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenGenHermVWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes the eigenvalues of the complex generalized hermitian-definite matrix
    /// pair (A, B), and stores them in `eval`, using the method outlined above. On output, `B`
    /// contains its Cholesky decomposition and `A` is destroyed.
    pub fn genhermv(&self, A: &MatrixComplexF64, B: &MatrixComplexF64, eval: &mut VectorF64,
                    evec: &mut MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_genhermv(ffi::FFI::unwrap(A), ffi::FFI::unwrap(B),
                                         ffi::FFI::unwrap(eval), ffi::FFI::unwrap(evec), self.w) }
    }
}

impl Drop for EigenGenHermVWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_genhermv_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_genhermv_workspace> for EigenGenHermVWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_genhermv_workspace) -> EigenGenHermVWorkspace {
        EigenGenHermVWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenGenHermVWorkspace) -> *mut ffi::gsl_eigen_genhermv_workspace {
        t.w
    }
}

pub struct EigenGenWorkspace {
    w: *mut ffi::gsl_eigen_gen_workspace
}

impl EigenGenWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n real generalized
    /// nonsymmetric eigensystems. The size of the workspace is O(n).
    pub fn new(n: usize) -> Option<EigenGenWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_gen_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenGenWorkspace {
                w: tmp,
            })
        }
    }

    /// This function sets some parameters which determine how the eigenvalue problem is solved in
    /// subsequent calls to gsl_eigen_gen.
    ///
    /// If compute_s is set to 1, the full Schur form S will be computed by gsl_eigen_gen. If it is
    /// set to 0, S will not be computed (this is the default setting). S is a quasi upper
    /// triangular matrix with 1-by-1 and 2-by-2 blocks on its diagonal. 1-by-1 blocks correspond to
    /// real eigenvalues, and 2-by-2 blocks correspond to complex eigenvalues.
    ///
    /// If compute_t is set to 1, the full Schur form T will be computed by gsl_eigen_gen. If it is
    /// set to 0, T will not be computed (this is the default setting). T is an upper triangular
    /// matrix with non-negative elements on its diagonal. Any 2-by-2 blocks in S will correspond to
    /// a 2-by-2 diagonal block in T.
    ///
    /// The balance parameter is currently ignored, since generalized balancing is not yet
    /// implemented.
    pub fn params(&self, compute_s: i32, compute_t: i32, balance: i32) {
        unsafe { ffi::gsl_eigen_gen_params(compute_s, compute_t, balance, self.w) }
    }

    /// This function computes the eigenvalues of the real generalized nonsymmetric matrix pair
    /// (A, B), and stores them as pairs in (alpha, beta), where alpha is complex and beta is real.
    /// If \beta_i is non-zero, then \lambda = \alpha_i / \beta_i is an eigenvalue. Likewise, if
    /// \alpha_i is non-zero, then \mu = \beta_i / \alpha_i is an eigenvalue of the alternate
    /// problem \mu A y = B y. The elements of beta are normalized to be non-negative.
    ///
    /// If S is desired, it is stored in A on output. If T is desired, it is stored in B on output.
    /// The ordering of eigenvalues in (alpha, beta) follows the ordering of the diagonal blocks in
    /// the Schur forms S and T. In rare cases, this function may fail to find all eigenvalues. If
    /// this occurs, an error code is returned.
    pub fn gen(&self, A: &MatrixF64, B: &MatrixF64, alpha: &mut VectorComplexF64,
               beta: &mut VectorF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_gen(ffi::FFI::unwrap(A), ffi::FFI::unwrap(B),
                                    ffi::FFI::unwrap(alpha), ffi::FFI::unwrap(beta), self.w) }
    }

    /// This function is identical to gsl_eigen_gen except that it also computes the left and right
    /// Schur vectors and stores them into `Q` and `Z` respectively.
    pub fn gen_QZ(&self, A: &MatrixF64, B: &MatrixF64, alpha: &mut VectorComplexF64,
                  beta: &mut VectorF64, Q: &mut MatrixF64, Z: &mut MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_gen_QZ(ffi::FFI::unwrap(A), ffi::FFI::unwrap(B),
                                       ffi::FFI::unwrap(alpha), ffi::FFI::unwrap(beta),
                                       ffi::FFI::unwrap(Q), ffi::FFI::unwrap(Z), self.w) }
    }
}

impl Drop for EigenGenWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_gen_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_gen_workspace> for EigenGenWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_gen_workspace) -> EigenGenWorkspace {
        EigenGenWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenGenWorkspace) -> *mut ffi::gsl_eigen_gen_workspace {
        t.w
    }
}

pub struct EigenGenVWorkspace {
    w: *mut ffi::gsl_eigen_genv_workspace
}

impl EigenGenVWorkspace {
    /// This function allocates a workspace for computing eigenvalues of n-by-n real generalized
    /// nonsymmetric eigensystems. The size of the workspace is O(n).
    pub fn new(n: usize) -> Option<EigenGenVWorkspace> {
        let tmp = unsafe { ffi::gsl_eigen_genv_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(EigenGenVWorkspace {
                w: tmp,
            })
        }
    }

    /// This function computes eigenvalues and right eigenvectors of the n-by-n real generalized
    /// nonsymmetric matrix pair (A, B). The eigenvalues are stored in (alpha, beta) and the
    /// eigenvectors are stored in evec. It first calls gsl_eigen_gen to compute the eigenvalues,
    /// Schur forms, and Schur vectors. Then it finds eigenvectors of the Schur forms and
    /// backtransforms them using the Schur vectors. The Schur vectors are destroyed in the process,
    /// but can be saved by using gsl_eigen_genv_QZ. The computed eigenvectors are normalized to
    /// have unit magnitude. On output, (A, B) contains the generalized Schur form (S, T). If
    /// gsl_eigen_gen fails, no eigenvectors are computed, and an error code is returned.
    pub fn genv(&self, A: &MatrixF64, B: &MatrixF64, alpha: &mut VectorComplexF64,
                beta: &mut VectorF64, evec: &mut MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_genv(ffi::FFI::unwrap(A), ffi::FFI::unwrap(B),
                                     ffi::FFI::unwrap(alpha), ffi::FFI::unwrap(beta),
                                     ffi::FFI::unwrap(evec), self.w) }
    }

    /// This function is identical to gsl_eigen_genv except that it also computes the left and right
    /// Schur vectors and stores them into `Q` and `Z` respectively.
    pub fn genv_QZ(&self, A: &MatrixF64, B: &MatrixF64, alpha: &mut VectorComplexF64,
                   beta: &mut VectorF64, evec: &mut MatrixComplexF64, Q: &mut MatrixF64,
                   Z: &mut MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_eigen_genv_QZ(ffi::FFI::unwrap(A), ffi::FFI::unwrap(B),
                                        ffi::FFI::unwrap(alpha), ffi::FFI::unwrap(beta),
                                        ffi::FFI::unwrap(evec), ffi::FFI::unwrap(Q),
                                        ffi::FFI::unwrap(Z), self.w) }
    }
}

impl Drop for EigenGenVWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_eigen_genv_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_eigen_genv_workspace> for EigenGenVWorkspace {
    fn wrap(t: *mut ffi::gsl_eigen_genv_workspace) -> EigenGenVWorkspace {
        EigenGenVWorkspace {
            w: t
        }
    }

    fn unwrap(t: &EigenGenVWorkspace) -> *mut ffi::gsl_eigen_genv_workspace {
        t.w
    }
}

#[test]
fn eigen_symmetric_workspace() {
    use MatrixF64;
    use VectorF64;

    let e = EigenSymmetricWorkspace::new(2).unwrap();
    let m = MatrixF64::new(2, 2).unwrap();

    let data = [5., 5., 1., 6.];
    m.set(0, 0, data[0]);
    m.set(0, 1, data[1]);
    m.set(1, 0, data[2]);
    m.set(1, 1, data[3]);
    let mut v = VectorF64::new(2).unwrap();
    e.symm(&m, &mut v);
    assert_eq!(&format!("{:.4} {:.4}", v.get(0), v.get(1)), "4.3820 6.6180");
}

// C code:
//
// ```ignore
// #include <gsl/gsl_eigen.h>
// #include <gsl/gsl_matrix.h>
//
// int main() {
//     gsl_eigen_symmv_workspace *t = gsl_eigen_symmv_alloc(3);
//     double data[] = {5., 5., 1., 6.};
//     gsl_matrix *m = gsl_matrix_calloc(2, 2);
//     gsl_matrix_set(m, 0, 0, data[0]);
//     gsl_matrix_set(m, 0, 1, data[1]);
//     gsl_matrix_set(m, 1, 0, data[2]);
//     gsl_matrix_set(m, 1, 1, data[3]);
//     gsl_vector *v = gsl_vector_calloc(2);
//     gsl_matrix *m2 = gsl_matrix_calloc(2, 2);
//     gsl_eigen_symmv(m, v, m2, t);
//     printf("%f %f\n", gsl_vector_get(v, 0), gsl_vector_get(v, 1));
//     printf("%f %f\n", gsl_matrix_get(m2, 0, 0), gsl_matrix_get(m2, 0, 1));
//     printf("%f %f\n", gsl_matrix_get(m2, 1, 0), gsl_matrix_get(m2, 1, 1));
//     gsl_eigen_symmv_free(t);
//     gsl_vector_free(v);
//     gsl_matrix_free(m);
//     gsl_matrix_free(m2);
//     return 0;
// }
// ```
#[test]
fn eigen_symmetric_vworkspace() {
    use MatrixF64;
    use VectorF64;

    let e = EigenSymmetricVWorkspace::new(3).unwrap();
    let data = [5., 5., 1., 6.];
    let m = MatrixF64::new(2, 2).unwrap();

    m.set(0, 0, data[0]);
    m.set(0, 1, data[1]);
    m.set(1, 0, data[2]);
    m.set(1, 1, data[3]);
    let mut m2 = MatrixF64::new(2, 2).unwrap();
    let mut v = VectorF64::new(2).unwrap();
    e.symmv(&m, &mut v, &mut m2);
    assert_eq!(&format!("{:.4} {:.4}", v.get(0), v.get(1)), "4.3820 6.6180");
    assert_eq!(&format!("{:.4} {:.4}", m2.get(0, 0), m2.get(0, 1)), "0.8507 0.5257");
    assert_eq!(&format!("{:.4} {:.4}", m2.get(1, 0), m2.get(1, 1)), "-0.5257 0.8507");
}

// C code:
//
// ```ignore
// #include <gsl/gsl_eigen.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_complex_math.h>
//
// int main() {
//     gsl_eigen_herm_workspace *t = gsl_eigen_herm_alloc(3);
//     gsl_matrix_complex *m = gsl_matrix_complex_calloc(2, 2);
//     gsl_complex c1 = gsl_complex_rect(5., 5.);
//     gsl_complex c2 = gsl_complex_rect(1., 4.);
//     gsl_complex c3 = gsl_complex_rect(2., 3.);
//     gsl_complex c4 = gsl_complex_rect(5., 7.);
//     gsl_matrix_complex_set(m, 0, 0, c1);
//     gsl_matrix_complex_set(m, 0, 1, c2);
//     gsl_matrix_complex_set(m, 1, 0, c3);
//     gsl_matrix_complex_set(m, 1, 1, c4);
//     gsl_vector *v = gsl_vector_calloc(2);
//     gsl_eigen_herm(m, v, t);
//     printf("%f %f\n", gsl_vector_get(v, 0), gsl_vector_get(v, 1));
//     gsl_eigen_herm_free(t);
//     gsl_vector_free(v);
//     gsl_matrix_complex_free(m);
//     return 0;
// }
// ```
#[test]
fn eigen_hermitian_workspace() {
    use MatrixComplexF64;
    use VectorF64;
    use ComplexF64;

    let e = EigenHermitianWorkspace::new(3).unwrap();
    let m = MatrixComplexF64::new(2, 2).unwrap();

    m.set(0, 0, &ComplexF64::rect(5., 5.));
    m.set(0, 1, &ComplexF64::rect(1., 4.));
    m.set(1, 0, &ComplexF64::rect(2., 3.));
    m.set(1, 1, &ComplexF64::rect(5., 7.));

    let mut v = VectorF64::new(2).unwrap();
    e.herm(&m, &mut v);
    assert_eq!(&format!("{:.4} {:.4}", v.get(0), v.get(1)), "8.6056 1.3944");
}

// C code:
//
// ```ignore
// #include <gsl/gsl_eigen.h>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_complex_math.h>
//
// int main() {
//     gsl_eigen_hermv_workspace *t = gsl_eigen_hermv_alloc(3);
//     gsl_matrix_complex *m = gsl_matrix_complex_calloc(2, 2);
//     gsl_complex c1 = gsl_complex_rect(5., 5.);
//     gsl_complex c2 = gsl_complex_rect(1., 4.);
//     gsl_complex c3 = gsl_complex_rect(2., 3.);
//     gsl_complex c4 = gsl_complex_rect(5., 7.);
//     gsl_matrix_complex_set(m, 0, 0, c1);
//     gsl_matrix_complex_set(m, 0, 1, c2);
//     gsl_matrix_complex_set(m, 1, 0, c3);
//     gsl_matrix_complex_set(m, 1, 1, c4);
//     gsl_vector *v = gsl_vector_calloc(2);
//     gsl_matrix_complex *m2 = gsl_matrix_complex_calloc(2, 2);
//     gsl_eigen_hermv(m, v, m2, t);
//     printf("%f %f\n", gsl_vector_get(v, 0), gsl_vector_get(v, 1));
//     printf("(%f, %f) (%f, %f)\n",
//            gsl_matrix_complex_get(m2, 0, 0).dat[0], gsl_matrix_complex_get(m2, 0, 0).dat[1],
//            gsl_matrix_complex_get(m2, 0, 1).dat[0], gsl_matrix_complex_get(m2, 0, 1).dat[1]);
//     printf("(%f, %f) (%f, %f)\n",
//            gsl_matrix_complex_get(m2, 1, 0).dat[0], gsl_matrix_complex_get(m2, 1, 0).dat[1],
//            gsl_matrix_complex_get(m2, 1, 1).dat[0], gsl_matrix_complex_get(m2, 1, 1).dat[1]);
//     gsl_eigen_hermv_free(t);
//     gsl_vector_free(v);
//     gsl_matrix_complex_free(m);
//     gsl_matrix_complex_free(m2);
//     return 0;
// }
// ```
#[test]
fn eigen_hermitian_vworkspace() {
    use ComplexF64;

    let e = EigenHermitianVWorkspace::new(3).unwrap();
    let m = MatrixComplexF64::new(2, 2).unwrap();

    m.set(0, 0, &ComplexF64::rect(5., 5.));
    m.set(0, 1, &ComplexF64::rect(1., 4.));
    m.set(1, 0, &ComplexF64::rect(2., 3.));
    m.set(1, 1, &ComplexF64::rect(5., 7.));

    let mut v = VectorF64::new(2).unwrap();
    let mut m2 = MatrixComplexF64::new(2, 2).unwrap();
    e.hermv(&m, &mut v, &mut m2);
    assert_eq!(&format!("{:.4} {:.4}", v.get(0), v.get(1)), "8.6056 1.3944");
    assert_eq!(&format!("({:.4}, {:.4}) ({:.4}, {:.4})",
                        m2.get(0, 0).dat[0], m2.get(0, 0).dat[1],
                        m2.get(0, 1).dat[0], m2.get(0, 1).dat[1]),
               "(0.7071, 0.0000) (0.7071, 0.0000)");
    assert_eq!(&format!("({:.4}, {:.4}) ({:.4}, {:.4})",
                        m2.get(1, 0).dat[0], m2.get(1, 0).dat[1],
                        m2.get(1, 1).dat[0], m2.get(1, 1).dat[1]),
               "(0.3922, 0.5883) (-0.3922, -0.5883)");
}
