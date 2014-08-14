//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::fmt;
use std::fmt::{Formatter,Show};
use types::{ComplexF32, ComplexF64};
use types::{VectorComplexF64, VectorComplexF32};
use ffi;
use enums;

pub struct MatrixComplexF64 {
    mat: *mut ffi::gsl_matrix_complex
}

impl MatrixComplexF64 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_matrix_complex {
        self.mat
    }

    /// Creates a new MatrixF64 with all elements set to zero
    /// 
    /// Example with n1 = 2 and n2 = 3 :
    /// 
    /// XX XX XX
    /// 
    /// XX XX XX
    pub fn new(n1: u64, n2: u64) -> Option<MatrixComplexF64> {
        let tmp = unsafe { ffi::gsl_matrix_complex_calloc(n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(MatrixComplexF64 {
                mat: tmp
            })
        }
    }

    /// This function returns the (i,j)-th element of the matrix.
    /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, y: u64, x: u64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(ffi::gsl_matrix_complex_get(self.mat as *const ffi::gsl_matrix_complex, y, x)) }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
    pub fn set(&self, y: u64, x: u64, value: &ComplexF64) -> &MatrixComplexF64 {
        unsafe { ffi::gsl_matrix_complex_set(self.mat, y, x, ::std::mem::transmute(*value)) };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    pub fn set_all(&self, x: &ComplexF64) -> &MatrixComplexF64 {
        unsafe { ffi::gsl_matrix_complex_set_all(self.mat, ::std::mem::transmute(*x)) };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    pub fn set_zero(&self) -> &MatrixComplexF64 {
        unsafe { ffi::gsl_matrix_complex_set_zero(self.mat) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    pub fn set_identity(&self) -> &MatrixComplexF64 {
        unsafe { ffi::gsl_matrix_complex_set_identity(self.mat) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
    pub fn copy_from(&self, other: &MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_memcpy(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_memcpy(other.mat, self.mat as *const ffi::gsl_matrix_complex) }
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&self, other: &MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_swap(self.mat, other.mat) }
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: u64) -> Option<(VectorComplexF64, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_complex_alloc((*self.mat).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_complex_get_row(tmp, self.mat as *const ffi::gsl_matrix_complex, y) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: u64) -> Option<(VectorComplexF64, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_complex_alloc((*self.mat).size1) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_complex_get_col(tmp, self.mat as *const ffi::gsl_matrix_complex, x) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    pub fn set_row(&self, y: u64, v: &VectorComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_set_row(self.mat, y, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_complex) }
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&self, x: u64, v: &VectorComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_set_col(self.mat, x, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_complex) }
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&self, y1: u64, y2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_swap_rows(self.mat, y1, y2) }
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&self, x1: u64, x2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_swap_columns(self.mat, x1, x2) }
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&self, i: u64, j: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_swap_rowcol(self.mat, i, j) }
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(MatrixComplexF64, enums::Value)> {
        let dest = unsafe { ffi::gsl_matrix_complex_alloc((*self.mat).size1, (*self.mat).size2) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_complex_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix_complex) };

            Some((MatrixComplexF64{mat: dest}, ret))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&self) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_transpose(self.mat) }
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&self, other: &MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_add(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&self, other: &MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_sub(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&self, other: &MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&self, other: &MatrixComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_div_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&self, x: &ComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_scale(self.mat, ::std::mem::transmute(*x)) }
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&self, x: &ComplexF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_add_constant(self.mat, ::std::mem::transmute(*x)) }
    }

    /// This function returns true if all the elements of the self matrix are stricly zero.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_isnull(self.mat as *const ffi::gsl_matrix_complex) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_ispos(self.mat as *const ffi::gsl_matrix_complex) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_isneg(self.mat as *const ffi::gsl_matrix_complex) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_isnonneg(self.mat as *const ffi::gsl_matrix_complex) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all elements of the two matrix are equal.
    pub fn equal(&self, other: &MatrixComplexF64) -> bool {
        match unsafe { ffi::gsl_matrix_complex_equal(self.mat as *const ffi::gsl_matrix_complex, other.mat as *const ffi::gsl_matrix_complex) } {
            1 => true,
            _ => false
        }
    }

    pub fn clone(&self) -> Option<MatrixComplexF64> {
        unsafe {
            if self.mat.is_null() {
                None
            } else {
                match MatrixComplexF64::new((*self.mat).size1, (*self.mat).size2) {
                    Some(m) => {
                        m.copy_from(self);
                        Some(m)
                    }
                    None => None
                }
            }
        }
    }
}

impl Drop for MatrixComplexF64 {
    fn drop(&mut self) {
        unsafe { ffi::gsl_matrix_complex_free(self.mat) };
        self.mat = ::std::ptr::mut_null();
    }
}

impl Show for MatrixComplexF64 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        unsafe {
            for y in range(0u64, (*self.mat).size1) {
                write!(f, "[");
                for x in range(0u64, (*self.mat).size2) {
                    if x < (*self.mat).size2 - 1 {
                        write!(f, "{}, ", self.get(y, x));
                    } else {
                        write!(f, "{}", self.get(y, x));
                    }
                }
                if y < (*self.mat).size1 - 1 {
                    write!(f, "]\n");
                }
            }
        }
        write!(f, "]")
    }
}

pub struct MatrixComplexF32 {
    mat: *mut ffi::gsl_matrix_complex_float
}

impl MatrixComplexF32 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_matrix_complex_float {
        self.mat
    }

    /// Creates a new MatrixF64 with all elements set to zero
    /// 
    /// Example with n1 = 2 and n2 = 3 :
    /// 
    /// XX XX XX
    /// 
    /// XX XX XX
    pub fn new(n1: u64, n2: u64) -> Option<MatrixComplexF32> {
        let tmp = unsafe { ffi::gsl_matrix_complex_float_calloc(n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(MatrixComplexF32 {
                mat: tmp
            })
        }
    }

    /// This function returns the (i,j)-th element of the matrix.
    /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, y: u64, x: u64) -> ComplexF32 {
        unsafe { ::std::mem::transmute(ffi::gsl_matrix_complex_float_get(self.mat as *const ffi::gsl_matrix_complex_float, y, x)) }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
    pub fn set(&self, y: u64, x: u64, value: &ComplexF32) -> &MatrixComplexF32 {
        unsafe { ffi::gsl_matrix_complex_float_set(self.mat, y, x, ::std::mem::transmute(*value)) };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    pub fn set_all(&self, x: &ComplexF32) -> &MatrixComplexF32 {
        unsafe { ffi::gsl_matrix_complex_float_set_all(self.mat, ::std::mem::transmute(*x)) };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    pub fn set_zero(&self) -> &MatrixComplexF32 {
        unsafe { ffi::gsl_matrix_complex_float_set_zero(self.mat) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    pub fn set_identity(&self) -> &MatrixComplexF32 {
        unsafe { ffi::gsl_matrix_complex_float_set_identity(self.mat) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
    pub fn copy_from(&self, other: &MatrixComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_memcpy(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &MatrixComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_memcpy(other.mat, self.mat as *const ffi::gsl_matrix_complex_float) }
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&self, other: &MatrixComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_swap(self.mat, other.mat) }
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: u64) -> Option<(VectorComplexF32, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_complex_float_alloc((*self.mat).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_complex_float_get_row(tmp, self.mat as *const ffi::gsl_matrix_complex_float, y) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: u64) -> Option<(VectorComplexF32, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_complex_float_alloc((*self.mat).size1) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_complex_float_get_col(tmp, self.mat as *const ffi::gsl_matrix_complex_float, x) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    pub fn set_row(&self, y: u64, v: &VectorComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_set_row(self.mat, y, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_complex_float) }
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&self, x: u64, v: &VectorComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_set_col(self.mat, x, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_complex_float) }
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&self, y1: u64, y2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_swap_rows(self.mat, y1, y2) }
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&self, x1: u64, x2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_swap_columns(self.mat, x1, x2) }
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&self, i: u64, j: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_swap_rowcol(self.mat, i, j) }
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(MatrixComplexF32, enums::Value)> {
        let dest = unsafe { ffi::gsl_matrix_complex_float_alloc((*self.mat).size1, (*self.mat).size2) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_complex_float_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix_complex_float) };

            Some((MatrixComplexF32{mat: dest}, ret))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&self) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_transpose(self.mat) }
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&self, other: &MatrixComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_add(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&self, other: &MatrixComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_sub(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&self, other: &MatrixComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&self, other: &MatrixComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_div_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&self, x: &ComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_scale(self.mat, ::std::mem::transmute(*x)) }
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&self, x: &ComplexF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_complex_float_add_constant(self.mat, ::std::mem::transmute(*x)) }
    }

    /// This function returns true if all the elements of the self matrix are stricly zero.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_float_isnull(self.mat as *const ffi::gsl_matrix_complex_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_float_ispos(self.mat as *const ffi::gsl_matrix_complex_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_float_isneg(self.mat as *const ffi::gsl_matrix_complex_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_complex_float_isnonneg(self.mat as *const ffi::gsl_matrix_complex_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all elements of the two matrix are equal.
    pub fn equal(&self, other: &MatrixComplexF32) -> bool {
        match unsafe { ffi::gsl_matrix_complex_float_equal(self.mat as *const ffi::gsl_matrix_complex_float,
            other.mat as *const ffi::gsl_matrix_complex_float) } {
            1 => true,
            _ => false
        }
    }

    pub fn clone(&self) -> Option<MatrixComplexF32> {
        unsafe {
            if self.mat.is_null() {
                None
            } else {
                match MatrixComplexF32::new((*self.mat).size1, (*self.mat).size2) {
                    Some(m) => {
                        m.copy_from(self);
                        Some(m)
                    }
                    None => None
                }
            }
        }
    }
}

impl Drop for MatrixComplexF32 {
    fn drop(&mut self) {
        unsafe { ffi::gsl_matrix_complex_float_free(self.mat) };
        self.mat = ::std::ptr::mut_null();
    }
}

impl Show for MatrixComplexF32 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        unsafe {
            for y in range(0u64, (*self.mat).size1) {
                write!(f, "[");
                for x in range(0u64, (*self.mat).size2) {
                    if x < (*self.mat).size2 - 1 {
                        write!(f, "{}, ", self.get(y, x));
                    } else {
                        write!(f, "{}", self.get(y, x));
                    }
                }
                if y < (*self.mat).size1 - 1 {
                    write!(f, "]\n");
                }
            }
        }
        write!(f, "]")
    }
}