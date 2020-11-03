//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use enums;
use ffi::{self, FFI};
use std::fmt;
use std::fmt::{Debug, Formatter};
use types::{ComplexF32, ComplexF64};
use types::{VectorComplexF32, VectorComplexF64};

ffi_wrapper!(
    MatrixComplexF64,
    *mut sys::gsl_matrix_complex,
    gsl_matrix_complex_free
);

impl MatrixComplexF64 {
    /// Creates a new MatrixF64 with all elements set to zero
    ///
    /// Example with n1 = 2 and n2 = 3 :
    ///
    /// XX XX XX
    ///
    /// XX XX XX
    pub fn new(n1: usize, n2: usize) -> Option<MatrixComplexF64> {
        let tmp = unsafe { sys::gsl_matrix_complex_calloc(n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function returns the (i,j)-th element of the matrix.
    /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, y: usize, x: usize) -> ComplexF64 {
        unsafe { ::std::mem::transmute(sys::gsl_matrix_complex_get(self.unwrap_shared(), y, x)) }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
    pub fn set(&mut self, y: usize, x: usize, value: &ComplexF64) -> &MatrixComplexF64 {
        unsafe {
            sys::gsl_matrix_complex_set(self.unwrap_unique(), y, x, ::std::mem::transmute(*value))
        };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    pub fn set_all(&mut self, x: &ComplexF64) -> &MatrixComplexF64 {
        unsafe { sys::gsl_matrix_complex_set_all(self.unwrap_unique(), ::std::mem::transmute(*x)) };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    pub fn set_zero(&mut self) -> &MatrixComplexF64 {
        unsafe { sys::gsl_matrix_complex_set_zero(self.unwrap_unique()) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    pub fn set_identity(&mut self) -> &MatrixComplexF64 {
        unsafe { sys::gsl_matrix_complex_set_identity(self.unwrap_unique()) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
    pub fn copy_from(&mut self, other: &MatrixComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_memcpy(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &mut MatrixComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_memcpy(other.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&mut self, other: &mut MatrixComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_swap(self.unwrap_unique(), other.unwrap_unique())
        })
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: usize) -> Option<(enums::Value, VectorComplexF64)> {
        let tmp = unsafe { sys::gsl_vector_complex_alloc((*self.unwrap_shared()).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::gsl_matrix_complex_get_row(tmp, self.unwrap_shared(), y) };

            Some((enums::Value::from(ret), ffi::FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: usize) -> Option<(enums::Value, VectorComplexF64)> {
        let tmp = unsafe { sys::gsl_vector_complex_alloc((*self.unwrap_shared()).size1) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::gsl_matrix_complex_get_col(tmp, self.unwrap_shared(), x) };

            Some((enums::Value::from(ret), ffi::FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    pub fn set_row(&mut self, y: usize, v: &VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_set_row(self.unwrap_unique(), y, v.unwrap_shared())
        })
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&mut self, x: usize, v: &VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_set_col(self.unwrap_unique(), x, v.unwrap_shared())
        })
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&mut self, y1: usize, y2: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_swap_rows(self.unwrap_unique(), y1, y2)
        })
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&mut self, x1: usize, x2: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_swap_columns(self.unwrap_unique(), x1, x2)
        })
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&mut self, i: usize, j: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_swap_rowcol(self.unwrap_unique(), i, j)
        })
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(enums::Value, MatrixComplexF64)> {
        let ptr = self.unwrap_shared();
        let dest = unsafe { sys::gsl_matrix_complex_alloc((*ptr).size2, (*ptr).size1) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { sys::gsl_matrix_complex_transpose_memcpy(dest, ptr) };

            Some((enums::Value::from(ret), MatrixComplexF64::wrap(dest)))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&mut self) -> enums::Value {
        enums::Value::from(unsafe { sys::gsl_matrix_complex_transpose(self.unwrap_unique()) })
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&mut self, other: &MatrixComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_add(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&mut self, other: &MatrixComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_sub(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&mut self, other: &MatrixComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_mul_elements(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&mut self, other: &MatrixComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_div_elements(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&mut self, x: &ComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_scale(self.unwrap_unique(), ::std::mem::transmute(*x))
        })
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&mut self, x: &ComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_add_constant(self.unwrap_unique(), ::std::mem::transmute(*x))
        })
    }

    /// This function returns true if all the elements of the self matrix are stricly zero.
    pub fn is_null(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_isnull(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_ispos(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_isneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_isnonneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all elements of the two matrix are equal.
    pub fn equal(&self, other: &MatrixComplexF64) -> bool {
        match unsafe { sys::gsl_matrix_complex_equal(self.unwrap_shared(), other.unwrap_shared()) }
        {
            1 => true,
            _ => false,
        }
    }

    pub fn clone(&self) -> Option<MatrixComplexF64> {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            None
        } else {
            unsafe {
                match MatrixComplexF64::new((*ptr).size1, (*ptr).size2) {
                    Some(mut m) => {
                        m.copy_from(self);
                        Some(m)
                    }
                    None => None,
                }
            }
        }
    }
}

impl Debug for MatrixComplexF64 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            write!(f, "<null>")
        } else {
            unsafe {
                for y in 0..(*ptr).size1 {
                    write!(f, "[")?;
                    for x in 0..(*ptr).size2 {
                        if x < (*ptr).size2 - 1 {
                            write!(f, "{:?}, ", self.get(y, x));
                        } else {
                            write!(f, "{:?}", self.get(y, x));
                        }
                    }
                    if y < (*ptr).size1 - 1 {
                        write!(f, "]\n");
                    }
                }
            }
            write!(f, "]")
        }
    }
}

ffi_wrapper!(
    MatrixComplexF32,
    *mut sys::gsl_matrix_complex_float,
    gsl_matrix_complex_float_free
);

impl MatrixComplexF32 {
    /// Creates a new MatrixF64 with all elements set to zero
    ///
    /// Example with n1 = 2 and n2 = 3 :
    ///
    /// XX XX XX
    ///
    /// XX XX XX
    pub fn new(n1: usize, n2: usize) -> Option<MatrixComplexF32> {
        let tmp = unsafe { sys::gsl_matrix_complex_float_calloc(n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(MatrixComplexF32::wrap(tmp))
        }
    }

    /// This function returns the (i,j)-th element of the matrix.
    /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, y: usize, x: usize) -> ComplexF32 {
        unsafe {
            ::std::mem::transmute(sys::gsl_matrix_complex_float_get(
                self.unwrap_shared(),
                y,
                x,
            ))
        }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
    pub fn set(&mut self, y: usize, x: usize, value: &ComplexF32) -> &MatrixComplexF32 {
        unsafe {
            sys::gsl_matrix_complex_float_set(
                self.unwrap_unique(),
                y,
                x,
                ::std::mem::transmute(*value),
            )
        };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    pub fn set_all(&mut self, x: &ComplexF32) -> &MatrixComplexF32 {
        unsafe {
            sys::gsl_matrix_complex_float_set_all(self.unwrap_unique(), ::std::mem::transmute(*x))
        };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    pub fn set_zero(&mut self) -> &MatrixComplexF32 {
        unsafe { sys::gsl_matrix_complex_float_set_zero(self.unwrap_unique()) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    pub fn set_identity(&mut self) -> &MatrixComplexF32 {
        unsafe { sys::gsl_matrix_complex_float_set_identity(self.unwrap_unique()) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
    pub fn copy_from(&mut self, other: &MatrixComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_memcpy(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &mut MatrixComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_memcpy(other.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&mut self, other: &mut MatrixComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_swap(self.unwrap_unique(), other.unwrap_unique())
        })
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: usize) -> Option<(enums::Value, VectorComplexF32)> {
        let ptr = self.unwrap_shared();
        let tmp = unsafe { sys::gsl_vector_complex_float_alloc((*ptr).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::gsl_matrix_complex_float_get_row(tmp, ptr, y) };

            Some((enums::Value::from(ret), ffi::FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: usize) -> Option<(enums::Value, VectorComplexF32)> {
        let ptr = self.unwrap_shared();
        let tmp = unsafe { sys::gsl_vector_complex_float_alloc((*ptr).size1) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::gsl_matrix_complex_float_get_col(tmp, ptr, x) };

            Some((enums::Value::from(ret), ffi::FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    pub fn set_row(&mut self, y: usize, v: &VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_set_row(self.unwrap_unique(), y, v.unwrap_shared())
        })
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&mut self, x: usize, v: &VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_set_col(self.unwrap_unique(), x, v.unwrap_shared())
        })
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&mut self, y1: usize, y2: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_swap_rows(self.unwrap_unique(), y1, y2)
        })
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&mut self, x1: usize, x2: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_swap_columns(self.unwrap_unique(), x1, x2)
        })
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&mut self, i: usize, j: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_swap_rowcol(self.unwrap_unique(), i, j)
        })
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(enums::Value, MatrixComplexF32)> {
        let ptr = self.unwrap_shared();
        let dest = unsafe { sys::gsl_matrix_complex_float_alloc((*ptr).size2, (*ptr).size1) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { sys::gsl_matrix_complex_float_transpose_memcpy(dest, ptr) };

            Some((enums::Value::from(ret), MatrixComplexF32::wrap(dest)))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&mut self) -> enums::Value {
        enums::Value::from(unsafe { sys::gsl_matrix_complex_float_transpose(self.unwrap_unique()) })
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&mut self, other: &MatrixComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_add(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&mut self, other: &MatrixComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_sub(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&mut self, other: &MatrixComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_mul_elements(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&mut self, other: &MatrixComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_div_elements(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&mut self, x: &ComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_scale(self.unwrap_unique(), ::std::mem::transmute(*x))
        })
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&mut self, x: &ComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_matrix_complex_float_add_constant(
                self.unwrap_unique(),
                ::std::mem::transmute(*x),
            )
        })
    }

    /// This function returns true if all the elements of the self matrix are stricly zero.
    pub fn is_null(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_float_isnull(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_float_ispos(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_float_isneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self matrix are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { sys::gsl_matrix_complex_float_isnonneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all elements of the two matrix are equal.
    pub fn equal(&self, other: &MatrixComplexF32) -> bool {
        match unsafe {
            sys::gsl_matrix_complex_float_equal(self.unwrap_shared(), other.unwrap_shared())
        } {
            1 => true,
            _ => false,
        }
    }

    pub fn clone(&self) -> Option<MatrixComplexF32> {
        unsafe {
            if self.unwrap_shared().is_null() {
                None
            } else {
                let ptr = self.unwrap_shared();
                match MatrixComplexF32::new((*ptr).size1, (*ptr).size2) {
                    Some(mut m) => {
                        m.copy_from(self);
                        Some(m)
                    }
                    None => None,
                }
            }
        }
    }
}

impl Debug for MatrixComplexF32 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            write!(f, "<null>")
        } else {
            unsafe {
                for y in 0..(*ptr).size1 {
                    write!(f, "[");
                    for x in 0..(*ptr).size2 {
                        if x < (*ptr).size2 - 1 {
                            write!(f, "{:?}, ", self.get(y, x));
                        } else {
                            write!(f, "{:?}", self.get(y, x));
                        }
                    }
                    if y < (*ptr).size1 - 1 {
                        write!(f, "]\n");
                    }
                }
            }
            write!(f, "]")
        }
    }
}
