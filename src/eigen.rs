//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
##References and Further Reading

Further information on the algorithms described in this section can be found in the following book,

G. H. Golub, C. F. Van Loan, Matrix Computations (3rd Ed, 1996), Johns Hopkins University Press, ISBN 0-8018-5414-8.
Further information on the generalized eigensystems QZ algorithm can be found in this paper,

C. Moler, G. Stewart, “An Algorithm for Generalized Matrix Eigenvalue Problems”, SIAM J. Numer. Anal., Vol 10, No 2, 1973.
Eigensystem routines for very large matrices can be found in the Fortran library LAPACK. The LAPACK library is described in,

LAPACK Users’ Guide (Third Edition, 1999), Published by SIAM, ISBN 0-89871-447-8.
http://www.netlib.org/lapack

The LAPACK source code can be found at the website above along with an online copy of the users guide.
!*/

use enums;
use ffi;
use types::{VectorF64, MatrixF64, MatrixComplexF64, VectorComplexF64};

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding real eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type
pub fn symmv_sort(eval: &mut VectorF64, evec: &mut MatrixF64, sort_type: ::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_symmv_sort(ffi::FFI::unwrap_unique(eval), ffi::FFI::unwrap_unique(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding complex eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type.
pub fn hermv_sort(eval: &mut VectorF64, evec: &mut MatrixComplexF64, sort_type: ::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_hermv_sort(ffi::FFI::unwrap_unique(eval), ffi::FFI::unwrap_unique(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding complex eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type. Only EigenSort::AbsAsc and
/// EigenSort::AbsDesc are supported due to the eigenvalues being complex.
pub fn nonsymmv_sort(eval: &mut VectorComplexF64, evec: &mut MatrixComplexF64, sort_type: ::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_nonsymmv_sort(ffi::FFI::unwrap_unique(eval), ffi::FFI::unwrap_unique(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding real eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type.
pub fn gensymmv_sort(eval: &mut VectorF64, evec: &mut MatrixF64, sort_type: ::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_gensymmv_sort(ffi::FFI::unwrap_unique(eval), ffi::FFI::unwrap_unique(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding complex eigenvectors stored in the
/// columns of the matrix evec into ascending or descending order according to the value of the parameter sort_type.
pub fn genhermv_sort(eval: &mut VectorF64, evec: &mut MatrixComplexF64, sort_type: ::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_genhermv_sort(ffi::FFI::unwrap_unique(eval), ffi::FFI::unwrap_unique(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vectors (alpha, beta) and the corresponding complex eigenvectors stored
/// in the columns of the matrix evec into ascending or descending order according to the value of the parameter sort_type. Only
/// EigenSort::AbsAsc and EigenSort::AbsDesc are supported due to the eigenvalues being complex.
pub fn genv_sort(alpha: &mut VectorComplexF64, beta: &mut VectorF64, evec: &mut MatrixComplexF64, sort_type: ::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_genv_sort(ffi::FFI::unwrap_unique(alpha), ffi::FFI::unwrap_unique(beta), ffi::FFI::unwrap_unique(evec), sort_type) }
}
