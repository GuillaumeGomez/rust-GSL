//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::paste::paste;
use crate::Value;
use ffi::FFI;
use std::fmt::{self, Debug, Formatter};

macro_rules! gsl_matrix_complex {
    ($rust_name:ident, $name:ident, $complex:ident, $complex_c:ident) => (
paste! {

use types::{$complex, [<Vector $complex>], [<Vector $complex View>]};

ffi_wrapper!(
    $rust_name,
    *mut sys::$name,
    [<$name _free>]
);

impl $rust_name {
    /// Creates a new MatrixF64.
    #[doc(alias = $name _alloc)]
    pub fn new(n1: usize, n2: usize) -> Option<Self> {
        let tmp = unsafe { sys::[<$name _alloc>](n1, n2) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// Creates a new MatrixF64 with all elements set to zero.
    #[doc(alias = $name _calloc)]
    pub fn new_with_init(n1: usize, n2: usize) -> Option<Self> {
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
    pub fn get(&self, y: usize, x: usize) -> $complex {
        unsafe { ::std::mem::transmute(sys::[<$name _get>](self.unwrap_shared(), y, x)) }
    }

    /// This function sets the value of the (i,j)-th element of the matrix to value.
    /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler
    /// is invoked.
    #[doc(alias = $name _set)]
    pub fn set(&mut self, y: usize, x: usize, value: &$complex) -> &Self {
        unsafe {
            sys::[<$name _set>](self.unwrap_unique(), y, x, ::std::mem::transmute(*value))
        };
        self
    }

    /// This function sets all the elements of the matrix to the value x.
    #[doc(alias = $name _set_all)]
    pub fn set_all(&mut self, x: &$complex) -> &Self {
        unsafe { sys::[<$name _set_all>](self.unwrap_unique(), ::std::mem::transmute(*x)) };
        self
    }

    /// This function sets all the elements of the matrix to zero.
    #[doc(alias = $name _set_zero)]
    pub fn set_zero(&mut self) -> &Self {
        unsafe { sys::[<$name _set_zero>](self.unwrap_unique()) };
        self
    }

    /// This function sets the elements of the matrix to the corresponding elements of the identity
    /// matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
    /// This applies to both square and rectangular matrices.
    #[doc(alias = $name _set_identity)]
    pub fn set_identity(&mut self) -> &Self {
        unsafe { sys::[<$name _set_identity>](self.unwrap_unique()) };
        self
    }

    /// This function copies the elements of the other matrix into the self matrix. The two matrices
    /// must have the same size.
    #[doc(alias = $name _memcpy)]
    pub fn copy_from(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe {
            sys::[<$name _memcpy>](self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function copies the elements of the self matrix into the other matrix. The two matrices
    /// must have the same size.
    #[doc(alias = $name _memcpy)]
    pub fn copy_to(&self, other: &mut $rust_name) -> Value {
        Value::from(unsafe {
            sys::[<$name _memcpy>](other.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function exchanges the elements of the matrices self and other by copying. The two
    /// matrices must have the same size.
    #[doc(alias = $name _swap)]
    pub fn swap(&mut self, other: &mut $rust_name) -> Value {
        Value::from(unsafe {
            sys::[<$name _swap>](self.unwrap_unique(), other.unwrap_unique())
        })
    }

    /// This function copies the elements of the y-th row of the matrix into the returned vector.
    #[doc(alias = $name _get_row)]
    pub fn get_row(&self, y: usize) -> Option<(Value, [<Vector $complex>])> {
        let tmp = unsafe { sys::[<$complex_c _alloc>](self.size2()) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::[<$name _get_row>](tmp, self.unwrap_shared(), y) };

            Some((Value::from(ret), FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the x-th column of the matrix into the returned vector.
    #[doc(alias = $name _get_col)]
    pub fn get_col(&self, x: usize) -> Option<(Value, [<Vector $complex>])> {
        let tmp = unsafe { sys::[<$complex_c _alloc>](self.size1()) };

        if tmp.is_null() {
            None
        } else {
            let ret = unsafe { sys::[<$name _get_col>](tmp, self.unwrap_shared(), x) };

            Some((Value::from(ret), FFI::wrap(tmp)))
        }
    }

    /// This function copies the elements of the vector v into the y-th row of the matrix.
    /// The length of the vector must be the same as the length of the row.
    #[doc(alias = $name _set_row)]
    pub fn set_row(&mut self, y: usize, v: &[<Vector $complex>]) -> Value {
        Value::from(unsafe {
            sys::[<$name _set_row>](self.unwrap_unique(), y, v.unwrap_shared())
        })
    }

    /// This function copies the elements of the vector v into the x-th column of the matrix.
    /// The length of the vector must be the same as the length of the column.
    #[doc(alias = $name _set_col)]
    pub fn set_col(&mut self, x: usize, v: &[<Vector $complex>]) -> Value {
        Value::from(unsafe {
            sys::[<$name _set_col>](self.unwrap_unique(), x, v.unwrap_shared())
        })
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

    /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must
    /// be square for this operation to be possible.
    #[doc(alias = $name _swap_rowcol)]
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
            let ret =
                unsafe { sys::[<$name _transpose_memcpy>](dest, self.unwrap_shared()) };

            Some((Value::from(ret), Self::wrap(dest)))
        }
    }

    /// This function replaces the matrix m by its transpose by copying the elements of the matrix
    /// in-place. The matrix must be square for this operation to be possible.
    #[doc(alias = $name _transpose)]
    pub fn transpose(&mut self) -> Value {
        Value::from(unsafe { sys::[<$name _transpose>](self.unwrap_unique()) })
    }

    /// This function adds the elements of the other matrix to the elements of the `self` matrix.
    /// The result self(i,j) <- self(i,j) + other(i,j) is stored in `self` and other remains
    /// unchanged. The two matrices must have the same dimensions.
    #[doc(alias = $name _add)]
    pub fn add(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe {
            sys::[<$name _add>](self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function subtracts the elements of the other matrix from the elements of the `self`
    /// matrix. The result self(i,j) <- self(i,j) - other(i,j) is stored in `self` and other remains
    /// unchanged. The two matrices must have the same dimensions.
    #[doc(alias = $name _sub)]
    pub fn sub(&mut self, other: &$rust_name) -> Value {
        Value::from(unsafe {
            sys::[<$name _sub>](self.unwrap_unique(), other.unwrap_shared())
        })
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
    pub fn scale(&mut self, x: &$complex) -> Value {
        Value::from(unsafe {
            sys::[<$name _scale>](self.unwrap_unique(), ::std::mem::transmute(*x))
        })
    }

    /// This function adds the constant value x to the elements of the self matrix. The result
    /// self(i,j) <- self(i,j) + x is stored in self.
    #[doc(alias = $name _add_constant)]
    pub fn add_constant(&mut self, x: &$complex) -> Value {
        Value::from(unsafe {
            sys::[<$name _add_constant>](self.unwrap_unique(), ::std::mem::transmute(*x))
        })
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
    pub fn row<F: FnOnce(Option<[<Vector $complex View>]>)>(&mut self, i: usize, f: F) {
        [<Vector $complex View>]::wrap(unsafe { sys::[<$name _row>](self.unwrap_unique(), i) }, f)
    }

    #[doc(alias = $name _column)]
    pub fn column<F: FnOnce(Option<[<Vector $complex View>]>)>(&mut self, j: usize, f: F) {
        [<Vector $complex View>]::wrap(unsafe { sys::[<$name _column>](self.unwrap_unique(), j) }, f)
    }

    #[doc(alias = $name _diagonal)]
    pub fn diagonal<F: FnOnce(Option<[<Vector $complex View>]>)>(&mut self, f: F) {
        [<Vector $complex View>]::wrap(unsafe { sys::[<$name _diagonal>](self.unwrap_unique()) }, f)
    }

    #[doc(alias = $name _subdiagonal)]
    pub fn subdiagonal<F: FnOnce(Option<[<Vector $complex View>]>)>(&mut self, k: usize, f: F) {
        [<Vector $complex View>]::wrap(unsafe { sys::[<$name _subdiagonal>](self.unwrap_unique(), k) }, f)
    }

    #[doc(alias = $name _superdiagonal)]
    pub fn superdiagonal<F: FnOnce(Option<[<Vector $complex View>]>)>(&mut self, k: usize, f: F) {
        [<Vector $complex View>]::wrap(unsafe { sys::[<$name _superdiagonal>](self.unwrap_unique(), k) }, f)
    }

    #[doc(alias = $name _subrow)]
    pub fn subrow<F: FnOnce(Option<[<Vector $complex View>]>)>(&mut self, i: usize, offset: usize, n: usize, f: F) {
        [<Vector $complex View>]::wrap(unsafe { sys::[<$name _subrow>](self.unwrap_unique(), i, offset, n) }, f)
    }

    #[doc(alias = $name _subcolumn)]
    pub fn subcolumn<F: FnOnce(Option<[<Vector $complex View>]>)>(&mut self, i: usize, offset: usize, n: usize, f: F) {
        [<Vector $complex View>]::wrap(unsafe { sys::[<$name _subcolumn>](self.unwrap_unique(), i, offset, n) }, f)
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
}

impl Debug for $rust_name {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        if self.unwrap_shared().is_null() {
            write!(f, "<null>")
        } else {
            let size1 = self.size1();
            let size2 = self.size2();
            for y in 0..size1 {
                write!(f, "[")?;
                for x in 0..size2 {
                    if x < size2 - 1 {
                        write!(f, "{:?}, ", self.get(y, x))?;
                    } else {
                        write!(f, "{:?}", self.get(y, x))?;
                    }
                }
                if y < size1 - 1 {
                    write!(f, "]\n")?;
                }
            }
            write!(f, "]")
        }
    }
}

} // end of paste! block
); // end of macro block
}

gsl_matrix_complex!(
    MatrixComplexF64,
    gsl_matrix_complex,
    ComplexF64,
    gsl_vector_complex
);
gsl_matrix_complex!(
    MatrixComplexF32,
    gsl_matrix_complex_float,
    ComplexF32,
    gsl_vector_complex_float
);
