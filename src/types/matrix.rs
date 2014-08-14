//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::fmt;
use std::fmt::{Formatter, Show};
use types::{VectorF64, VectorF32};
use ffi;
use enums;

pub struct MatrixF64 {
    mat: *mut ffi::gsl_matrix
}

impl MatrixF64 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_matrix {
        self.mat
    }

    /// Creates a new MatrixF64 with all elements set to zero
    /// 
    /// Example with n1 = 2 and n2 = 3 :
    /// 
    /// XX XX XX
    /// 
    /// XX XX XX
    pub fn new(n1: u64, n2: u64) -> Option<MatrixF64> {
        let tmp = unsafe { ffi::gsl_matrix_calloc(n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(MatrixF64 {
                mat: tmp
            })
        }
    }

    /// This function returns the (i,j)-th element of the matrix.
    /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, y: u64, x: u64) -> f64 {
        unsafe { ffi::gsl_matrix_get(self.mat as *const ffi::gsl_matrix, y, x) }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
    pub fn set(&self, y: u64, x: u64, value: f64) -> &MatrixF64 {
        unsafe { ffi::gsl_matrix_set(self.mat, y, x, value) };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    pub fn set_all(&self, x: f64) -> &MatrixF64 {
        unsafe { ffi::gsl_matrix_set_all(self.mat, x) };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    pub fn set_zero(&self) -> &MatrixF64 {
        unsafe { ffi::gsl_matrix_set_zero(self.mat) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    pub fn set_identity(&self) -> &MatrixF64 {
        unsafe { ffi::gsl_matrix_set_identity(self.mat) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
    pub fn copy_from(&self, other: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_memcpy(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_memcpy(other.mat, self.mat as *const ffi::gsl_matrix) }
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&self, other: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_swap(self.mat, other.mat) }
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: u64) -> Option<(VectorF64, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_alloc((*self.mat).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_get_row(tmp, self.mat as *const ffi::gsl_matrix, y) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: u64) -> Option<(VectorF64, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_alloc((*self.mat).size1) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_get_col(tmp, self.mat as *const ffi::gsl_matrix, x) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    pub fn set_row(&self, y: u64, v: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_set_row(self.mat, y, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&self, x: u64, v: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_set_col(self.mat, x, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&self, y1: u64, y2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_swap_rows(self.mat, y1, y2) }
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&self, x1: u64, x2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_swap_columns(self.mat, x1, x2) }
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&self, i: u64, j: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_swap_rowcol(self.mat, i, j) }
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(MatrixF64, enums::Value)> {
        let dest = unsafe { ffi::gsl_matrix_alloc((*self.mat).size1, (*self.mat).size2) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix) };

            Some((MatrixF64 {mat: dest}, ret))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&self) -> enums::Value {
        unsafe { ffi::gsl_matrix_transpose(self.mat) }
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&self, other: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_add(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&self, other: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_sub(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&self, other: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&self, other: &MatrixF64) -> enums::Value {
        unsafe { ffi::gsl_matrix_div_elements(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&self, x: f64) -> enums::Value {
        unsafe { ffi::gsl_matrix_scale(self.mat, x) }
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&self, x: f64) -> enums::Value {
        unsafe { ffi::gsl_matrix_add_constant(self.mat, x) }
    }

    /// This function returns the maximum value in the self matrix.
    pub fn max(&self) -> f64 {
        unsafe { ffi::gsl_matrix_max(self.mat as *const ffi::gsl_matrix) }
    }

    /// This function returns the minimum value in the self matrix.
    pub fn min(&self) -> f64 {
        unsafe { ffi::gsl_matrix_min(self.mat as *const ffi::gsl_matrix) }
    }

    /// This function returns the minimum and maximum values in the self matrix, storing them in min_out and max_out.
    pub fn minmax(&self, min_out: &mut f64, max_out: &mut f64) {
        unsafe { ffi::gsl_matrix_minmax(self.mat as *const ffi::gsl_matrix, min_out, max_out) }
    }

    /// This function returns the indices of the maximum value in the self matrix, storing them in imax and jmax.
    /// When there are several equal maximum elements then the first element found is returned, searching in row-major order.
    pub fn max_index(&self) -> (u64, u64) {
        let mut imax = 0u64;
        let mut jmax = 0u64;

        unsafe { ffi::gsl_matrix_max_index(self.mat as *const ffi::gsl_matrix, &mut imax, &mut jmax) };
        (imax, jmax)
    }

    /// This function returns the indices of the minimum value in the self matrix, storing them in imin and jmin.
    /// When there are several equal minimum elements then the first element found is returned, searching in row-major order.
    pub fn min_index(&self) -> (u64, u64) {
        let mut imax = 0u64;
        let mut jmax = 0u64;

        unsafe { ffi::gsl_matrix_min_index(self.mat as *const ffi::gsl_matrix, &mut imax, &mut jmax) };
        (imax, jmax)
    }

    /// This function returns the indices of the minimum and maximum values in the self matrix, storing them in (imin,jmin) and (imax,jmax).
    /// When there are several equal minimum or maximum elements then the first elements found are returned, searching in row-major order.
    pub fn minmax_index(&self) -> (u64, u64, u64, u64) {
        let mut imin = 0u64;
        let mut jmin = 0u64;
        let mut imax = 0u64;
        let mut jmax = 0u64;

        unsafe { ffi::gsl_matrix_minmax_index(self.mat as *const ffi::gsl_matrix, &mut imin, &mut jmin, &mut imax, &mut jmax) };
        (imin, jmin, imax, jmax)
    }

    /// This function returns true if all the elements of the self matrix are stricly zero.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_matrix_isnull(self.mat as *const ffi::gsl_matrix) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_matrix_ispos(self.mat as *const ffi::gsl_matrix) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_isneg(self.mat as *const ffi::gsl_matrix) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_isnonneg(self.mat as *const ffi::gsl_matrix) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all elements of the two matrix are equal.
    pub fn equal(&self, other: &MatrixF64) -> bool {
        match unsafe { ffi::gsl_matrix_equal(self.mat as *const ffi::gsl_matrix, other.mat as *const ffi::gsl_matrix) } {
            1 => true,
            _ => false
        }
    }

    pub fn clone(&self) -> Option<MatrixF64> {
        unsafe {
            if self.mat.is_null() {
                None
            } else {
                match MatrixF64::new((*self.mat).size1, (*self.mat).size2) {
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

impl Drop for MatrixF64 {
    fn drop(&mut self) {
        unsafe { ffi::gsl_matrix_free(self.mat) };
        self.mat = ::std::ptr::mut_null();
    }
}

impl Show for MatrixF64 {
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

pub struct MatrixF32 {
    mat: *mut ffi::gsl_matrix_float
}

impl MatrixF32 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_matrix_float {
        self.mat
    }

    /// Creates a new MatrixF64 with all elements set to zero
    /// 
    /// Example with n1 = 2 and n2 = 3 :
    /// 
    /// XX XX XX
    /// 
    /// XX XX XX
    pub fn new(n1: u64, n2: u64) -> Option<MatrixF32> {
        let tmp = unsafe { ffi::gsl_matrix_float_calloc(n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(MatrixF32 {
                mat: tmp
            })
        }
    }

    /// This function returns the (i,j)-th element of the matrix.
    /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, y: u64, x: u64) -> f32 {
        unsafe { ffi::gsl_matrix_float_get(self.mat as *const ffi::gsl_matrix_float, y, x) }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
    pub fn set(&self, y: u64, x: u64, value: f32) -> &MatrixF32 {
        unsafe { ffi::gsl_matrix_float_set(self.mat, y, x, value) };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    pub fn set_all(&self, x: f32) -> &MatrixF32 {
        unsafe { ffi::gsl_matrix_float_set_all(self.mat, x) };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    pub fn set_zero(&self) -> &MatrixF32 {
        unsafe { ffi::gsl_matrix_float_set_zero(self.mat) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    pub fn set_identity(&self) -> &MatrixF32 {
        unsafe { ffi::gsl_matrix_float_set_identity(self.mat) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
    pub fn copy_from(&self, other: &MatrixF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_memcpy(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &MatrixF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_memcpy(other.mat, self.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&self, other: &MatrixF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_swap(self.mat, other.mat) }
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: u64) -> Option<(VectorF32, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_float_alloc((*self.mat).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_float_get_row(tmp, self.mat as *const ffi::gsl_matrix_float, y) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: u64) -> Option<(VectorF32, enums::Value)> {
        let tmp = unsafe { ffi::gsl_vector_float_alloc((*self.mat).size1) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_float_get_col(tmp, self.mat as *const ffi::gsl_matrix_float, x) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    pub fn set_row(&self, y: u64, v: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_set_row(self.mat, y, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_float) }
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&self, x: u64, v: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_set_col(self.mat, x, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_float) }
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&self, y1: u64, y2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_swap_rows(self.mat, y1, y2) }
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&self, x1: u64, x2: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_swap_columns(self.mat, x1, x2) }
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&self, i: u64, j: u64) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_swap_rowcol(self.mat, i, j) }
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(MatrixF32, enums::Value)> {
        let dest = unsafe { ffi::gsl_matrix_float_alloc((*self.mat).size1, (*self.mat).size2) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_float_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix_float) };

            Some((MatrixF32{mat: dest}, ret))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&self) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_transpose(self.mat) }
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&self, other: &MatrixF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_add(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&self, other: &MatrixF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_sub(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&self, other: &MatrixF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&self, other: &MatrixF32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_div_elements(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&self, x: f32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_scale(self.mat, x) }
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&self, x: f32) -> enums::Value {
        unsafe { ffi::gsl_matrix_float_add_constant(self.mat, x) }
    }

    /// This function returns the maximum value in the self matrix.
    pub fn max(&self) -> f32 {
        unsafe { ffi::gsl_matrix_float_max(self.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function returns the minimum value in the self matrix.
    pub fn min(&self) -> f32 {
        unsafe { ffi::gsl_matrix_float_min(self.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function returns the minimum and maximum values in the self matrix, storing them in min_out and max_out.
    pub fn minmax(&self, min_out: &mut f32, max_out: &mut f32) {
        unsafe { ffi::gsl_matrix_float_minmax(self.mat as *const ffi::gsl_matrix_float, min_out, max_out) }
    }

    /// This function returns the indices of the maximum value in the self matrix, storing them in imax and jmax.
    /// When there are several equal maximum elements then the first element found is returned, searching in row-major order.
    pub fn max_index(&self) -> (u64, u64) {
        let mut imax = 0u64;
        let mut jmax = 0u64;

        unsafe { ffi::gsl_matrix_float_max_index(self.mat as *const ffi::gsl_matrix_float, &mut imax, &mut jmax) };
        (imax, jmax)
    }

    /// This function returns the indices of the minimum value in the self matrix, storing them in imin and jmin.
    /// When there are several equal minimum elements then the first element found is returned, searching in row-major order.
    pub fn min_index(&self) -> (u64, u64) {
        let mut imax = 0u64;
        let mut jmax = 0u64;

        unsafe { ffi::gsl_matrix_float_min_index(self.mat as *const ffi::gsl_matrix_float, &mut imax, &mut jmax) };
        (imax, jmax)
    }

    /// This function returns the indices of the minimum and maximum values in the self matrix, storing them in (imin,jmin) and (imax,jmax).
    /// When there are several equal minimum or maximum elements then the first elements found are returned, searching in row-major order.
    pub fn minmax_index(&self) -> (u64, u64, u64, u64) {
        let mut imin = 0u64;
        let mut jmin = 0u64;
        let mut imax = 0u64;
        let mut jmax = 0u64;

        unsafe { ffi::gsl_matrix_float_minmax_index(self.mat as *const ffi::gsl_matrix_float, &mut imin, &mut jmin, &mut imax, &mut jmax) };
        (imin, jmin, imax, jmax)
    }

    /// This function returns true if all the elements of the self matrix are stricly zero.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_matrix_float_isnull(self.mat as *const ffi::gsl_matrix_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_matrix_float_ispos(self.mat as *const ffi::gsl_matrix_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_float_isneg(self.mat as *const ffi::gsl_matrix_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_matrix_float_isnonneg(self.mat as *const ffi::gsl_matrix_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all elements of the two matrix are equal.
    pub fn equal(&self, other: &MatrixF32) -> bool {
        match unsafe { ffi::gsl_matrix_float_equal(self.mat as *const ffi::gsl_matrix_float, other.mat as *const ffi::gsl_matrix_float) } {
            1 => true,
            _ => false
        }
    }

    pub fn clone(&self) -> Option<MatrixF32> {
        unsafe {
            if self.mat.is_null() {
                None
            } else {
                match MatrixF32::new((*self.mat).size1, (*self.mat).size2) {
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

impl Drop for MatrixF32 {
    fn drop(&mut self) {
        unsafe { ffi::gsl_matrix_float_free(self.mat) };
        self.mat = ::std::ptr::mut_null();
    }
}

impl Show for MatrixF32 {
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