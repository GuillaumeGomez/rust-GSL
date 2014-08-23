//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use enums;
use ffi;
use types::{VectorF64, MatrixF64, MatrixComplexF64, VectorComplexF64};

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding real eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type
pub fn symmv_sort(eval: &VectorF64, evec: &MatrixF64, sort_type: enums::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_symmv_sort(ffi::FFI::unwrap(eval), ffi::FFI::unwrap(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding complex eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type.
pub fn hermv_sort(eval: &VectorF64, evec: &MatrixComplexF64, sort_type: enums::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_hermv_sort(ffi::FFI::unwrap(eval), ffi::FFI::unwrap(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding complex eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type. Only EigenSort::AbsAsc and
/// EigenSort::AbsDesc are supported due to the eigenvalues being complex.
pub fn nonsymmv_sort(eval: &VectorComplexF64, evec: &MatrixComplexF64, sort_type: enums::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_nonsymmv_sort(ffi::FFI::unwrap(eval), ffi::FFI::unwrap(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding real eigenvectors stored in the columns
/// of the matrix evec into ascending or descending order according to the value of the parameter sort_type.
pub fn gensymmv_sort(eval: &VectorF64, evec: &MatrixF64, sort_type: enums::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_gensymmv_sort(ffi::FFI::unwrap(eval), ffi::FFI::unwrap(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vector eval and the corresponding complex eigenvectors stored in the
/// columns of the matrix evec into ascending or descending order according to the value of the parameter sort_type.
pub fn genhermv_sort(eval: &VectorF64, evec: &MatrixComplexF64, sort_type: enums::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_genhermv_sort(ffi::FFI::unwrap(eval), ffi::FFI::unwrap(evec), sort_type) }
}

/// This function simultaneously sorts the eigenvalues stored in the vectors (alpha, beta) and the corresponding complex eigenvectors stored
/// in the columns of the matrix evec into ascending or descending order according to the value of the parameter sort_type. Only
/// EigenSort::AbsAsc and EigenSort::AbsDesc are supported due to the eigenvalues being complex.
pub fn genv_sort(alpha: &VectorComplexF64, beta: &VectorF64, evec: &MatrixComplexF64, sort_type: enums::EigenSort) -> enums::Value {
    unsafe { ffi::gsl_eigen_genv_sort(ffi::FFI::unwrap(alpha), ffi::FFI::unwrap(beta), ffi::FFI::unwrap(evec), sort_type) }
}