//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Matrices

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

## References and Further Reading

The block, vector and matrix objects in GSL follow the valarray model of C++. A description of this model can be found in the following
reference,

B. Stroustrup, The C++ Programming Language (3rd Ed), Section 22.4 Vector Arithmetic. Addison-Wesley 1997, ISBN 0-201-88954-4.
!*/

use crate::paste::paste;
use crate::Value;
use ffi::{self, FFI};
use std::fmt::{self, Debug, Formatter};
use std::marker::PhantomData;
use types::{VectorF32, VectorF64, VectorI32, VectorU32};
use types::{VectorF32View, VectorF64View, VectorI32View, VectorU32View};

macro_rules! gsl_matrix {
    ($rust_name:ident, $name:ident, $rust_ty:ident, $vec_name:ident, $vec_c_name:ident) => (
paste! {
pub struct $rust_name {
    mat: *mut sys::$name,
    can_free: bool,
}

// macro_rules! lol {
//     ($doc:expr, $($t:tt)*) => (
//         #[doc(alias = $doc)]
//         $($t)*
//     );
// }

impl $rust_name {
    #[doc = "Creates a new " $rust_name " with all elements set to zero"]
    #[doc(alias = $name _calloc)]
    pub fn new(n1: usize, n2: usize) -> Option<$rust_name> {
        let tmp = unsafe { sys::[<$name _calloc>](n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function returns the (i,j)-th element of the matrix.
    /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is
    /// invoked and 0 is returned.
    #[doc(alias = $name _get)]
    pub fn get(&self, y: usize, x: usize) -> $rust_ty {
        unsafe { sys::[<$name _get>](self.unwrap_shared(), y, x) }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler
    /// is invoked.
    #[doc(alias = $name _set)]
    pub fn set(&mut self, y: usize, x: usize, value: $rust_ty) -> &$rust_name {
        unsafe { sys::[<$name _set>](self.unwrap_unique(), y, x, value) };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    #[doc(alias = $name _set_all)]
    pub fn set_all(&mut self, x: $rust_ty) -> &$rust_name {
        unsafe { sys::[<$name _set_all>](self.unwrap_unique(), x) };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    #[doc(alias = $name _set_zero)]
    pub fn set_zero(&mut self) -> &$rust_name {
        unsafe { sys::[<$name _set_zero>](self.unwrap_unique()) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity
    /// matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    #[doc(alias = $name _set_identity)]
    pub fn set_identity(&mut self) -> &$rust_name {
        unsafe { sys::[<$name _set_identity>](self.unwrap_unique()) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices
    /// must have the same size.
    #[doc(alias = $name _memcpy)]
    pub fn copy_from(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe { sys::[<$name _memcpy>](self.unwrap_unique(), other.unwrap_shared()) })
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices
    /// must have the same size.
    #[doc(alias = $name _memcpy)]
    pub fn copy_to(&self, other: &mut $rust_name) -> Value {
        Value::from(unsafe { sys::[<$name _memcpy>](other.unwrap_unique(), self.unwrap_shared()) })
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two
    /// matrices must have the same size.
    #[doc(alias = $name _swap)]
    pub fn swap(&mut self, other: &mut $rust_name) -> Value {
        Value::from(unsafe { sys::[<$name _swap>](self.unwrap_unique(), other.unwrap_unique()) })
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    #[doc(alias = $name _get_row)]
    pub fn get_row(&self, y: usize) -> Option<(Value, $vec_name)> {
        let tmp = unsafe { sys::[<$vec_c_name _alloc>](self.size2()) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::[<$name _get_row>](tmp, self.unwrap_shared(), y) };

            Some((Value::from(ret), ffi::FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    #[doc(alias = $name _get_col)]
    pub fn get_col(&self, x: usize) -> Option<(Value, $vec_name)> {
        let tmp = unsafe { sys::[<$vec_c_name _alloc>](self.size1()) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::[<$name _get_col>](tmp, self.unwrap_shared(), x) };

            Some((Value::from(ret), ffi::FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    #[doc(alias = $name _set_row)]
    pub fn set_row(&mut self, y: usize, v: &$vec_name) -> Value {
        Value::from(unsafe { sys::[<$name _set_row>](self.unwrap_unique(), y, v.unwrap_shared()) })
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    #[doc(alias = $name _set_col)]
    pub fn set_col(&mut self, x: usize, v: &$vec_name) -> Value {
        Value::from(unsafe { sys::[<$name _set_col>](self.unwrap_unique(), x, v.unwrap_shared()) })
    }

    /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
    #[doc(alias = $name _swap_rows)]
    pub fn swap_rows(&mut self, y1: usize, y2: usize) -> Value {
        Value::from(unsafe { sys::[<$name _swap_rows>](self.unwrap_unique(), y1, y2) })
    }

    /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
    #[doc(alias = $name _swap_columns)]
    pub fn swap_columns(&mut self, x1: usize, x2: usize) -> Value {
        Value::from(unsafe { sys::[<$name _swap_columns>](self.unwrap_unique(), x1, x2) })
    }

    /// This function exchanges the i-th row and j-th column of the matrix in-place.
    /// The matrix must be square for this operation to be possible.
    #[doc(alias = $name _swap_row_col)]
    pub fn swap_row_col(&mut self, i: usize, j: usize) -> Value {
        Value::from(unsafe { sys::[<$name _swap_rowcol>](self.unwrap_unique(), i, j) })
    }

    /// This function returns the transpose of the matrix by copying the elements into it.
    /// This function works for all matrices provided that the dimensions of the matrix dest match
    /// the transposed dimensions of the matrix.
    #[doc(alias = $name _transpose_memcpy)]
    pub fn transpose_memcpy(&self) -> Option<(Value, $rust_name)> {
        let dest = unsafe { sys::[<$name _alloc>](self.size2(), self.size1()) };

        if dest.is_null() {
            None
        } else {
            let ret = unsafe { sys::[<$name _transpose_memcpy>](dest, self.unwrap_shared()) };

            Some((Value::from(ret), $rust_name::wrap(dest)))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix
    /// in-place. The matrix must be square for this operation to be possible.
    #[doc(alias = $name _transpose)]
    pub fn transpose(&mut self) -> Value {
        Value::from(unsafe { sys::[<$name _transpose>](self.unwrap_unique()) })
    }

    /// This function adds the elements of the other matrix to the elements of the self matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains
    /// unchanged. The two matrices must have the same dimensions.
    #[doc(alias = $name _add)]
    pub fn add(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe { sys::[<$name _add>](self.unwrap_unique(), other.unwrap_shared()) })
    }

    /// This function subtracts the elements of the other matrix from the elements of the self
    /// matrix. The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains
    /// unchanged. The two matrices must have the same dimensions.
    #[doc(alias = $name _sub)]
    pub fn sub(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe { sys::[<$name _sub>](self.unwrap_unique(), other.unwrap_shared()) })
    }

    /// This function multiplies the elements of the self matrix by the elements of the other
    /// matrix. The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains
    /// unchanged. The two matrices must have the same dimensions.
    #[doc(alias = $name _mul_elements)]
    pub fn mul_elements(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe {
            sys::[<$name _mul_elements>](self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function divides the elements of the self matrix by the elements of the other matrix.
    /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains
    /// unchanged. The two matrices must have the same dimensions.
    #[doc(alias = $name _div_elements)]
    pub fn div_elements(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe {
            sys::[<$name _div_elements>](self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self matrix by the constant factor x. The
    /// result self(i,j) <- x self(i,j) is stored in self.
    #[doc(alias = $name _scale)]
    pub fn scale(&mut self, x: f64) -> Value {
        Value::from(unsafe { sys::[<$name _scale>](self.unwrap_unique(), x) })
    }

    /// This function adds the constant value x to the elements of the self matrix. The result
    /// self(i,j) <- self(i,j) + x is stored in self.
    #[doc(alias = $name _add_constant)]
    pub fn add_constant(&mut self, x: f64) -> Value {
        Value::from(unsafe { sys::[<$name _add_constant>](self.unwrap_unique(), x) })
    }

    #[doc(alias = $name _add_constant)]
    pub fn add_diagonal(&mut self, x: f64) -> Value {
        Value::from(unsafe { sys::[<$name _add_constant>](self.unwrap_unique(), x) })
    }

    /// This function returns the maximum value in the self matrix.
    #[doc(alias = $name _max)]
    pub fn max(&self) -> $rust_ty {
        unsafe { sys::[<$name _max>](self.unwrap_shared()) }
    }

    /// This function returns the minimum value in the self matrix.
    #[doc(alias = $name _min)]
    pub fn min(&self) -> $rust_ty {
        unsafe { sys::[<$name _min>](self.unwrap_shared()) }
    }

    /// This function returns the minimum and maximum values in the self matrix.
    #[doc(alias = $name _minmax)]
    pub fn minmax(&self) -> ($rust_ty, $rust_ty) {
        let mut min_out = 0 as _;
        let mut max_out = 0 as _;
        unsafe { sys::[<$name _minmax>](self.unwrap_shared(), &mut min_out, &mut max_out) };
        (min_out, max_out)
    }

    /// This function returns the indices of the maximum value in the self matrix. When there are
    /// several equal maximum elements then the first element found is returned, searching in
    /// row-major order.
    #[doc(alias = $name _max_index)]
    pub fn max_index(&self) -> (usize, usize) {
        let mut imax = 0;
        let mut jmax = 0;

        unsafe { sys::[<$name _max_index>](self.unwrap_shared(), &mut imax, &mut jmax) };
        (imax, jmax)
    }

    /// This function returns the indices of the minimum value in the self matrix. When there are
    /// several equal minimum elements then the first element found is returned, searching in row
    /// major order.
    #[doc(alias = $name _min_index)]
    pub fn min_index(&self) -> (usize, usize) {
        let mut imax = 0;
        let mut jmax = 0;

        unsafe { sys::[<$name _min_index>](self.unwrap_shared(), &mut imax, &mut jmax) };
        (imax, jmax)
    }

    /// This function returns the indices of the minimum and maximum values in the self matrix. When
    /// there are several equal minimum or maximum elements then the first elements found are
    /// returned, searching in row-major order.
    #[doc(alias = $name _minmax_index)]
    pub fn minmax_index(&self) -> (usize, usize, usize, usize) {
        let mut imin = 0;
        let mut jmin = 0;
        let mut imax = 0;
        let mut jmax = 0;

        unsafe {
            sys::[<$name _minmax_index>](
                self.unwrap_shared(),
                &mut imin,
                &mut jmin,
                &mut imax,
                &mut jmax,
            )
        };
        (imin, jmin, imax, jmax)
    }

    /// This function returns true if all the elements of the self matrix are stricly zero.
    #[doc(alias = $name _isnull)]
    pub fn is_null(&self) -> bool {
        unsafe { sys::[<$name _isnull>](self.unwrap_shared()) == 1 }
    }

    /// This function returns true if all the elements of the self matrix are stricly positive.
    #[doc(alias = $name _ispos)]
    pub fn is_pos(&self) -> bool {
        unsafe { sys::[<$name _ispos>](self.unwrap_shared()) == 1 }
    }

    /// This function returns true if all the elements of the self matrix are stricly negative.
    #[doc(alias = $name _isneg)]
    pub fn is_neg(&self) -> bool {
        unsafe { sys::[<$name _isneg>](self.unwrap_shared()) == 1 }
    }

    /// This function returns true if all the elements of the self matrix are stricly non-negative.
    #[doc(alias = $name _isnonneg)]
    pub fn is_non_neg(&self) -> bool {
        unsafe { sys::[<$name _isnonneg>](self.unwrap_shared()) == 1 }
    }

    /// This function returns true if all elements of the two matrix are equal.
    #[doc(alias = $name _equal)]
    pub fn equal(&self, other: &$rust_name) -> bool {
        unsafe { sys::[<$name _equal>](self.unwrap_shared(), other.unwrap_shared()) == 1 }
    }

    #[doc(alias = $name _row)]
    pub fn row<F: FnOnce(Option<[<$vec_name View>]>)>(&mut self, i: usize, f: F) {
        [<$vec_name View>]::wrap(unsafe { sys::[<$name _row>](self.unwrap_unique(), i) }, f)
    }

    #[doc(alias = $name _column)]
    pub fn column<F: FnOnce(Option<[<$vec_name View>]>)>(&mut self, j: usize, f: F) {
        [<$vec_name View>]::wrap(unsafe { sys::[<$name _column>](self.unwrap_unique(), j) }, f)
    }

    #[doc(alias = $name _diagonal)]
    pub fn diagonal<F: FnOnce(Option<[<$vec_name View>]>)>(&mut self, f: F) {
        [<$vec_name View>]::wrap(unsafe { sys::[<$name _diagonal>](self.unwrap_unique()) }, f)
    }

    #[doc(alias = $name _subdiagonal)]
    pub fn subdiagonal<F: FnOnce(Option<[<$vec_name View>]>)>(&mut self, k: usize, f: F) {
        [<$vec_name View>]::wrap(unsafe { sys::[<$name _subdiagonal>](self.unwrap_unique(), k) }, f)
    }

    #[doc(alias = $name _superdiagonal)]
    pub fn superdiagonal<F: FnOnce(Option<[<$vec_name View>]>)>(&mut self, k: usize, f: F) {
        [<$vec_name View>]::wrap(unsafe { sys::[<$name _superdiagonal>](self.unwrap_unique(), k) }, f)
    }

    #[doc(alias = $name _subrow)]
    pub fn subrow<F: FnOnce(Option<[<$vec_name View>]>)>(&mut self, i: usize, offset: usize, n: usize, f: F) {
        [<$vec_name View>]::wrap(unsafe { sys::[<$name _subrow>](self.unwrap_unique(), i, offset, n) }, f)
    }

    #[doc(alias = $name _subcolumn)]
    pub fn subcolumn<F: FnOnce(Option<[<$vec_name View>]>)>(&mut self, i: usize, offset: usize, n: usize, f: F) {
        [<$vec_name View>]::wrap(unsafe { sys::[<$name _subcolumn>](self.unwrap_unique(), i, offset, n) }, f)
    }

    #[doc(alias = $name _submatrix)]
    pub fn submatrix<'a>(
        &'a mut self,
        k1: usize,
        k2: usize,
        n1: usize,
        n2: usize,
    ) -> [<$rust_name View>]<'a> {
        [<$rust_name View>]::from_matrix(self, k1, k2, n1, n2)
    }

    pub fn size1(&self) -> usize {
        if self.unwrap_shared().is_null() {
            0
        } else {
            unsafe { (*self.unwrap_shared()).size1 }
        }
    }

    pub fn size2(&self) -> usize {
        if self.unwrap_shared().is_null() {
            0
        } else {
            unsafe { (*self.unwrap_shared()).size2 }
        }
    }

    pub fn clone(&self) -> Option<Self> {
        if self.unwrap_shared().is_null() {
            None
        } else {
            match Self::new(self.size1(), self.size2()) {
                Some(mut m) => {
                    m.copy_from(self);
                    Some(m)
                }
                None => None,
            }
        }
    }

    #[doc(hidden)]
    pub fn is_ptr_null(&self) -> bool {
        self.unwrap_shared().is_null()
    }
}

impl Drop for $rust_name {
    #[doc(alias = $name _free)]
    fn drop(&mut self) {
        if self.can_free {
            unsafe { sys::[<$name _free>](self.mat) };
            self.mat = ::std::ptr::null_mut();
        }
    }
}

impl FFI<sys::$name> for $rust_name {
    fn wrap(mat: *mut sys::$name) -> Self {
        Self {
            mat,
            can_free: true,
        }
    }

    fn soft_wrap(mat: *mut sys::$name) -> Self {
        Self {
            mat,
            can_free: false,
        }
    }

    fn unwrap_shared(&self) -> *const sys::$name {
        self.mat as *const _
    }

    fn unwrap_unique(&mut self) -> *mut sys::$name {
        self.mat
    }
}

impl Debug for $rust_name {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            write!(f, "<null>")
        } else {
            let size1 = self.size1();
            let size2 = self.size2();
            for y in 0..size1 {
                write!(f, "[");
                for x in 0..size2 {
                    if x < size2 - 1 {
                        write!(f, "{}, ", self.get(y, x));
                    } else {
                        write!(f, "{}", self.get(y, x));
                    }
                }
                if y < size1 - 1 {
                    write!(f, "]\n");
                }
            }
            write!(f, "]")
        }
    }
}

pub struct [<$rust_name View>]<'a> {
    mat: sys::[<$name _view>],
    #[allow(dead_code)]
    phantom: PhantomData<&'a ()>,
}

impl<'a> [<$rust_name View>]<'a> {
    /// These functions return a matrix view of a submatrix of the matrix m. The upper-left element
    /// of the submatrix is the element (k1,k2) of the original matrix. The submatrix has n1 rows
    /// and n2 columns. The physical number of columns in memory given by tda is unchanged.
    /// Mathematically, the (i,j)-th element of the new matrix is given by,
    ///
    /// m'(i,j) = m->data[(k1*m->tda + k2) + i*m->tda + j]
    ///
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    ///
    /// The data pointer of the returned matrix struct is set to null if the combined parameters
    /// (i,j,n1,n2,tda) overrun the ends of the original
    /// matrix.
    ///
    /// The new matrix view is only a view of the block underlying the existing matrix, m. The
    /// block containing the elements of m is not
    /// owned by the new matrix view. When the view goes out of scope the original matrix m and its
    /// block will continue to exist. The original memory can only be deallocated by freeing the
    /// original matrix. Of course, the original matrix should not be deallocated while the view
    /// is still in use.
    ///
    /// The function gsl_matrix_const_submatrix is equivalent to gsl_matrix_submatrix but can be
    /// used for matrices which are declared const.
    #[doc(alias = $name _submatrix)]
    pub fn from_matrix(
        m: &'a mut $rust_name,
        k1: usize,
        k2: usize,
        n1: usize,
        n2: usize,
    ) -> Self {
        unsafe {
            Self {
                mat: sys::[<$name _submatrix>](m.mat, k1, k2, n1, n2),
                phantom: PhantomData,
            }
        }
    }

    /// These functions return a matrix view of the array base. The matrix has n1 rows and n2
    /// columns. The physical number of columns in memory is also given by n2. Mathematically, the
    /// (i,j)-th element of the new matrix is given by,
    ///
    /// m'(i,j) = base[i*n2 + j]
    ///
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    ///
    /// The new matrix is only a view of the array base. When the view goes out of scope the
    /// original array base will continue to exist. The original memory can only be deallocated by
    /// freeing the original array. Of course, the original array should not be deallocated while
    /// the view is still in use.
    ///
    /// The function gsl_matrix_const_view_array is equivalent to gsl_matrix_view_array but can be
    /// used for matrices which are declared const.
    #[doc(alias = $name _view_array)]
    pub fn from_array(base: &'a mut [$rust_ty], n1: usize, n2: usize) -> Self {
        assert!(
            n1 * n2 <= base.len() as _,
            "n1 * n2 cannot be longer than base"
        );
        unsafe {
            Self {
                mat: sys::[<$name _view_array>](base.as_mut_ptr(), n1, n2),
                phantom: PhantomData,
            }
        }
    }

    /// These functions return a matrix view of the array base with a physical number of columns tda
    /// which may differ from the corresponding dimension of the matrix. The matrix has n1 rows and
    /// n2 columns, and the physical number of columns in memory is given by tda. Mathematically,
    /// the (i,j)-th element of the new matrix is given by,
    ///
    /// m'(i,j) = base[i*tda + j]
    ///
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    ///
    /// The new matrix is only a view of the array base. When the view goes out of scope the
    /// original array base will continue to exist. The original memory can only be deallocated by
    /// freeing the original array. Of course, the original array should not be deallocated while
    /// the view is still in use.
    ///
    /// The function gsl_matrix_const_view_array_with_tda is equivalent to
    /// gsl_matrix_view_array_with_tda but can be used for matrices which are declared const.
    #[doc(alias = $name _view_array_with_tda)]
    pub fn from_array_with_tda(base: &'a mut [$rust_ty], n1: usize, n2: usize, tda: usize) -> Self {
        unsafe {
            Self {
                mat: sys::[<$name _view_array_with_tda>](base.as_mut_ptr(), n1, n2, tda),
                phantom: PhantomData,
            }
        }
    }

    /// These functions return a matrix view of the vector v. The matrix has n1 rows and n2 columns.
    /// The vector must have unit stride. The physical number of columns in memory is also given by
    /// n2. Mathematically, the (i,j)-th element of the new matrix is given by,
    ///
    /// m'(i,j) = v->data[i*n2 + j]
    ///
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    ///
    /// The new matrix is only a view of the vector v. When the view goes out of scope the original
    /// vector v will continue to exist. The original memory can only be deallocated by freeing the
    /// original vector. Of course, the original vector should not be deallocated while the view
    /// is still in use.
    ///
    /// The function gsl_matrix_const_view_vector is equivalent to gsl_matrix_view_vector but can be
    /// used for matrices which are declared const.
    #[doc(alias = $name _view_vector)]
    pub fn from_vector(v: &'a mut $vec_name, n1: usize, n2: usize) -> Self {
        unsafe {
            Self {
                mat: sys::[<$name _view_vector>](v.unwrap_unique(), n1, n2),
                phantom: PhantomData,
            }
        }
    }

    /// These functions return a matrix view of the vector v with a physical number of columns tda
    /// which may differ from the corresponding matrix dimension. The vector must have unit stride.
    /// The matrix has n1 rows and n2 columns, and the physical number of columns in memory is given
    /// by tda. Mathematically, the (i,j)-th element of the new matrix is given by,
    ///
    /// m'(i,j) = v->data[i*tda + j]
    ///
    /// where the index i runs from 0 to n1-1 and the index j runs from 0 to n2-1.
    ///
    /// The new matrix is only a view of the vector v. When the view goes out of scope the original
    /// vector v will continue to exist. The original memory can only be deallocated by freeing the
    /// original vector. Of course, the original vector should not be deallocated while the view
    /// is still in use.
    ///
    /// The function gsl_matrix_const_view_vector_with_tda is equivalent to
    /// gsl_matrix_view_vector_with_tda but can be used for matrices which are declared const.
    #[doc(alias = $name _view_vector_with_tda)]
    pub fn from_vector_with_tda(v: &'a mut $vec_name, n1: usize, n2: usize, tda: usize) -> Self {
        unsafe {
            Self {
                mat: sys::[<$name _view_vector_with_tda>](v.unwrap_unique(), n1, n2, tda),
                phantom: PhantomData,
            }
        }
    }

    pub fn matrix<F: FnOnce(Option<&$rust_name>)>(&self, f: F) {
        let tmp = &self.mat.matrix;
        let tmp_mat = $rust_name::soft_wrap(tmp as *const _ as usize as *mut _);
        if tmp_mat.is_ptr_null() {
            f(None)
        } else {
            f(Some(&tmp_mat))
        }
    }

    pub fn matrix_mut<F: FnOnce(Option<&mut $rust_name>)>(&mut self, f: F) {
        let tmp = &mut self.mat.matrix;
        let mut tmp_mat = $rust_name::soft_wrap(tmp as *mut _);
        if tmp_mat.is_ptr_null() {
            f(None)
        } else {
            f(Some(&mut tmp_mat))
        }
    }
} // end of impl block
} // end of paste! block

    ); // end of the gsl_matrix macro
}

gsl_matrix!(
    MatrixF32,
    gsl_matrix_float,
    f32,
    VectorF32,
    gsl_vector_float
);
gsl_matrix!(MatrixF64, gsl_matrix, f64, VectorF64, gsl_vector);
gsl_matrix!(MatrixI32, gsl_matrix_int, i32, VectorI32, gsl_vector_int);
gsl_matrix!(MatrixU32, gsl_matrix_uint, u32, VectorU32, gsl_vector_uint);
