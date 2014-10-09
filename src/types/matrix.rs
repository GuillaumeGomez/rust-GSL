//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Matrices

Matrices are defined by a gsl_matrix structure which describes a generalized slice of a block. Like a vector it represents a set of 
elements in an area of memory, but uses two indices instead of one.

The gsl_matrix structure contains six components, the two dimensions of the matrix, a physical dimension, a pointer to the memory where 
the elements of the matrix are stored, data, a pointer to the block owned by the matrix block, if any, and an ownership flag, owner. The 
physical dimension determines the memory layout and can differ from the matrix dimension to allow the use of submatrices. The gsl_matrix 
structure is very simple and looks like this,

```C
typedef struct
{
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  gsl_block * block;
  int owner;
} gsl_matrix;
```

Matrices are stored in row-major order, meaning that each row of elements forms a contiguous block in memory. This is the standard 
“C-language ordering” of two-dimensional arrays. Note that FORTRAN stores arrays in column-major order. The number of rows is size1. The 
range of valid row indices runs from 0 to size1-1. Similarly size2 is the number of columns. The range of valid column indices runs from 
0 to size2-1. The physical row dimension tda, or trailing dimension, specifies the size of a row of the matrix as laid out in memory.

For example, in the following matrix size1 is 3, size2 is 4, and tda is 8. The physical memory layout of the matrix begins in the top 
left hand-corner and proceeds from left to right along each row in turn.

00 01 02 03 XX XX XX XX
10 11 12 13 XX XX XX XX
20 21 22 23 XX XX XX XX

Each unused memory location is represented by “XX”. The pointer data gives the location of the first element of the matrix in memory. The 
pointer block stores the location of the memory block in which the elements of the matrix are located (if any). If the matrix owns this 
block then the owner field is set to one and the block will be deallocated when the matrix is freed. If the matrix is only a slice of a 
block owned by another object then the owner field is zero and any underlying block will not be freed.

##References and Further Reading

The block, vector and matrix objects in GSL follow the valarray model of C++. A description of this model can be found in the following 
reference,

B. Stroustrup, The C++ Programming Language (3rd Ed), Section 22.4 Vector Arithmetic. Addison-Wesley 1997, ISBN 0-201-88954-4.
!*/

use std::fmt;
use std::fmt::{Formatter, Show};
use types::{VectorF64, VectorF32};
use ffi;
use enums;

pub struct MatrixView {
    mat: ffi::gsl_matrix
}

impl MatrixView {
    /// These functions return a matrix view of a submatrix of the matrix m. The upper-left element of the submatrix is the element (k1,k2)
    /// of the original matrix. The submatrix has n1 rows and n2 columns. The physical number of columns in memory given by tda is unchanged.
    /// Mathematically, the (i,j)-th element of the new matrix is given by,
    /// 
    /// m'(i,j) = m->data[(k1*m->tda + k2) + i*m->tda + j]
    /// 
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    /// 
    /// The data pointer of the returned matrix struct is set to null if the combined parameters (i,j,n1,n2,tda) overrun the ends of the original
    /// matrix.
    /// 
    /// The new matrix view is only a view of the block underlying the existing matrix, m. The block containing the elements of m is not
    /// owned by the new matrix view. When the view goes out of scope the original matrix m and its block will continue to exist. The original
    /// memory can only be deallocated by freeing the original matrix. Of course, the original matrix should not be deallocated while the view
    /// is still in use.
    /// 
    /// The function gsl_matrix_const_submatrix is equivalent to gsl_matrix_submatrix but can be used for matrices which are declared const.
    pub fn from_matrix(m: &MatrixF64, k1: u64, k2: u64, n1: u64, n2: u64) -> MatrixView {
        unsafe {
            MatrixView {
                mat: ffi::gsl_matrix_submatrix(m.mat, k1, k2, n1, n2).mat
            }
        }
    }

    /// These functions return a matrix view of the array base. The matrix has n1 rows and n2 columns. The physical number of columns in memory
    /// is also given by n2. Mathematically, the (i,j)-th element of the new matrix is given by,
    /// 
    /// m'(i,j) = base[i*n2 + j]
    /// 
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    /// 
    /// The new matrix is only a view of the array base. When the view goes out of scope the original array base will continue to exist. The
    /// original memory can only be deallocated by freeing the original array. Of course, the original array should not be deallocated while
    /// the view is still in use.
    /// 
    /// The function gsl_matrix_const_view_array is equivalent to gsl_matrix_view_array but can be used for matrices which are declared const.
    pub fn from_array(base: &mut [f64], n1: u64, n2: u64) -> MatrixView {
        unsafe {
            MatrixView {
                mat: ffi::gsl_matrix_view_array(base.as_mut_ptr(), n1, n2).mat
            }
        }
    }

    /// These functions return a matrix view of the array base with a physical number of columns tda which may differ from the corresponding
    /// dimension of the matrix. The matrix has n1 rows and n2 columns, and the physical number of columns in memory is given by tda.
    /// Mathematically, the (i,j)-th element of the new matrix is given by,
    /// 
    /// m'(i,j) = base[i*tda + j]
    ///
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    /// 
    /// The new matrix is only a view of the array base. When the view goes out of scope the original array base will continue to exist. The
    /// original memory can only be deallocated by freeing the original array. Of course, the original array should not be deallocated while
    /// the view is still in use.
    /// 
    /// The function gsl_matrix_const_view_array_with_tda is equivalent to gsl_matrix_view_array_with_tda but can be used for matrices which
    /// are declared const.
    pub fn from_array_with_tda(base: &mut [f64], n1: u64, n2: u64, tda: u64) -> MatrixView {
        unsafe {
            MatrixView {
                mat: ffi::gsl_matrix_view_array_with_tda(base.as_mut_ptr(), n1, n2, tda).mat
            }
        }
    }

    /// These functions return a matrix view of the vector v. The matrix has n1 rows and n2 columns. The vector must have unit stride. The
    /// physical number of columns in memory is also given by n2. Mathematically, the (i,j)-th element of the new matrix is given by,
    /// 
    /// m'(i,j) = v->data[i*n2 + j]
    /// 
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    /// 
    /// The new matrix is only a view of the vector v. When the view goes out of scope the original vector v will continue to exist. The original
    /// memory can only be deallocated by freeing the original vector. Of course, the original vector should not be deallocated while the view
    /// is still in use.
    /// 
    /// The function gsl_matrix_const_view_vector is equivalent to gsl_matrix_view_vector but can be used for matrices which are declared const.
    pub fn from_vector(v: &VectorF64, n1: u64, n2: u64) -> MatrixView {
        unsafe {
            MatrixView {
                mat: ffi::gsl_matrix_view_vector(ffi::FFI::unwrap(v), n1, n2).mat
            }
        }
    }

    /// These functions return a matrix view of the vector v with a physical number of columns tda which may differ from the corresponding
    /// matrix dimension. The vector must have unit stride. The matrix has n1 rows and n2 columns, and the physical number of columns in
    /// memory is given by tda. Mathematically, the (i,j)-th element of the new matrix is given by,
    /// 
    /// m'(i,j) = v->data[i*tda + j]
    /// 
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    /// 
    /// The new matrix is only a view of the vector v. When the view goes out of scope the original vector v will continue to exist. The original
    /// memory can only be deallocated by freeing the original vector. Of course, the original vector should not be deallocated while the view
    /// is still in use.
    /// 
    /// The function gsl_matrix_const_view_vector_with_tda is equivalent to gsl_matrix_view_vector_with_tda but can be used for matrices which
    /// are declared const.
    pub fn from_vector_with_tda(v: &VectorF64, n1: u64, n2: u64, tda: u64) -> MatrixView {
        unsafe {
            MatrixView {
                mat: ffi::gsl_matrix_view_vector_with_tda(ffi::FFI::unwrap(v), n1, n2, tda).mat
            }
        }
    }

    pub fn matrix(&mut self) -> MatrixF64 {
        unsafe {
            MatrixF64 {
                mat: ::std::mem::transmute(&mut self.mat),
                can_free: false
            }
        }
    }
}

pub struct MatrixF64 {
    mat: *mut ffi::gsl_matrix,
    can_free: bool
}

impl MatrixF64 {
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
                mat: tmp,
                can_free: true
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
    pub fn copy_from(&self, other: &MatrixF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_memcpy(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &MatrixF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_memcpy(other.mat, self.mat as *const ffi::gsl_matrix) }
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&self, other: &MatrixF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_swap(self.mat, other.mat) }
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: u64) -> Option<(VectorF64, enums::value::Value)> {
        let tmp = unsafe { ffi::gsl_vector_alloc((*self.mat).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_get_row(tmp, self.mat as *const ffi::gsl_matrix, y) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: u64) -> Option<(VectorF64, enums::value::Value)> {
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
    pub fn set_row(&self, y: u64, v: &VectorF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_set_row(self.mat, y, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&self, x: u64, v: &VectorF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_set_col(self.mat, x, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&self, y1: u64, y2: u64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_swap_rows(self.mat, y1, y2) }
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&self, x1: u64, x2: u64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_swap_columns(self.mat, x1, x2) }
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&self, i: u64, j: u64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_swap_rowcol(self.mat, i, j) }
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(MatrixF64, enums::value::Value)> {
        let dest = unsafe { ffi::gsl_matrix_alloc((*self.mat).size1, (*self.mat).size2) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix) };

            Some((MatrixF64 {mat: dest, can_free: true}, ret))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&self) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_transpose(self.mat) }
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&self, other: &MatrixF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_add(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&self, other: &MatrixF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_sub(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&self, other: &MatrixF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&self, other: &MatrixF64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_div_elements(self.mat, other.mat as *const ffi::gsl_matrix) }
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&self, x: f64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_scale(self.mat, x) }
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&self, x: f64) -> enums::value::Value {
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
        if self.can_free {
            unsafe { ffi::gsl_matrix_free(self.mat) };
            self.mat = ::std::ptr::null_mut();
        }
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

impl ffi::FFI<ffi::gsl_matrix> for MatrixF64 {
    fn wrap(r: *mut ffi::gsl_matrix) -> MatrixF64 {
        MatrixF64 {
            mat: r,
            can_free: true
        }
    }

    fn unwrap(m: &MatrixF64) -> *mut ffi::gsl_matrix {
        m.mat
    }
}

pub struct MatrixF32 {
    mat: *mut ffi::gsl_matrix_float,
    can_free: bool
}

impl MatrixF32 {
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
                mat: tmp,
                can_free: true
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
    pub fn copy_from(&self, other: &MatrixF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_memcpy(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
    pub fn copy_to(&self, other: &MatrixF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_memcpy(other.mat, self.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
    pub fn swap(&self, other: &MatrixF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_swap(self.mat, other.mat) }
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    pub fn get_row(&self, y: u64) -> Option<(VectorF32, enums::value::Value)> {
        let tmp = unsafe { ffi::gsl_vector_float_alloc((*self.mat).size2) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_float_get_row(tmp, self.mat as *const ffi::gsl_matrix_float, y) };

            Some((ffi::FFI::wrap(tmp), ret))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    pub fn get_col(&self, x: u64) -> Option<(VectorF32, enums::value::Value)> {
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
    pub fn set_row(&self, y: u64, v: &VectorF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_set_row(self.mat, y, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_float) }
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    pub fn set_col(&self, x: u64, v: &VectorF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_set_col(self.mat, x, ffi::FFI::unwrap(v) as *const ffi::gsl_vector_float) }
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    pub fn swap_rows(&self, y1: u64, y2: u64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_swap_rows(self.mat, y1, y2) }
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    pub fn swap_columns(&self, x1: u64, x2: u64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_swap_columns(self.mat, x1, x2) }
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
    pub fn swap_row_col(&self, i: u64, j: u64) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_swap_rowcol(self.mat, i, j) }
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
    pub fn transpose_memcpy(&self) -> Option<(MatrixF32, enums::value::Value)> {
        let dest = unsafe { ffi::gsl_matrix_float_alloc((*self.mat).size1, (*self.mat).size2) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { ffi::gsl_matrix_float_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix_float) };

            Some((MatrixF32{
                    mat: dest,
                    can_free: true
                }, ret))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    pub fn transpose(&self) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_transpose(self.mat) }
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn add(&self, other: &MatrixF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_add(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function subtracts the elements of the other matrix from the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn sub(&self, other: &MatrixF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_sub(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function multiplies the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn mul_elements(&self, other: &MatrixF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
    pub fn div_elements(&self, other: &MatrixF32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_div_elements(self.mat, other.mat as *const ffi::gsl_matrix_float) }
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
    pub fn scale(&self, x: f32) -> enums::value::Value {
        unsafe { ffi::gsl_matrix_float_scale(self.mat, x) }
    }

    /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
    pub fn add_constant(&self, x: f32) -> enums::value::Value {
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
        if self.can_free {
            unsafe { ffi::gsl_matrix_float_free(self.mat) };
            self.mat = ::std::ptr::null_mut();
        }
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

impl ffi::FFI<ffi::gsl_matrix_float> for MatrixF32 {
    fn wrap(r: *mut ffi::gsl_matrix_float) -> MatrixF32 {
        MatrixF32 {
            mat: r,
            can_free: true
        }
    }

    fn unwrap(m: &MatrixF32) -> *mut ffi::gsl_matrix_float {
        m.mat
    }
}