//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Vectors

Vectors are defined by a gsl_vector structure which describes a slice of a block. Different vectors can be created which point to the 
same block. A vector slice is a set of equally-spaced elements of an area of memory.

The gsl_vector structure contains five components, the size, the stride, a pointer to the memory where the elements are stored, data, a 
pointer to the block owned by the vector, block, if any, and an ownership flag, owner. The structure is very simple and looks like this,

```C
typedef struct
{
  size_t size;
  size_t stride;
  double * data;
  gsl_block * block;
  int owner;
} gsl_vector;
```

The size is simply the number of vector elements. The range of valid indices runs from 0 to size-1. The stride is the step-size from one 
element to the next in physical memory, measured in units of the appropriate datatype. The pointer data gives the location of the first 
element of the vector in memory. The pointer block stores the location of the memory block in which the vector elements are located (if 
any). If the vector owns this block then the owner field is set to one and the block will be deallocated when the vector is freed. If the 
vector points to a block owned by another object then the owner field is zero and any underlying block will not be deallocated with the 
vector.
!*/

use std::fmt;
use std::fmt::{Formatter, Debug};
use ffi;
use enums;

pub struct VectorView {
    v: ffi::gsl_vector_view
}

impl VectorView {
    /// These functions return a vector view of a subvector of another vector v. The start of the new vector is offset by offset elements
    /// from the start of the original vector. The new vector has n elements. Mathematically, the i-th element of the new vector v’ is given by,
    /// 
    /// v'(i) = v->data[(offset + i)*v->stride]
    /// 
    /// where the index i runs from 0 to n-1.
    /// 
    /// The data pointer of the returned vector struct is set to null if the combined parameters (offset,n) overrun the end of the original
    /// vector.
    /// 
    /// The new vector is only a view of the block underlying the original vector, v. The block containing the elements of v is not owned by
    /// the new vector. When the view goes out of scope the original vector v and its block will continue to exist. The original memory can
    /// only be deallocated by freeing the original vector. Of course, the original vector should not be deallocated while the view is still
    /// in use.
    /// 
    /// The function gsl_vector_const_subvector is equivalent to gsl_vector_subvector but can be used for vectors which are declared const.
    pub fn from_vector(v: &VectorF64, offset: usize, n: usize) -> VectorView {
        unsafe {
            VectorView {
                v: ffi::gsl_vector_subvector(v.vec, offset, n)
            }
        }
    }

    /// These functions return a vector view of a subvector of another vector v with an additional stride argument. The subvector is formed
    /// in the same way as for gsl_vector_subvector but the new vector has n elements with a step-size of stride from one element to the
    /// next in the original vector. Mathematically, the i-th element of the new vector v’ is given by,
    /// 
    /// v'(i) = v->data[(offset + i*stride)*v->stride]
    /// where the index i runs from 0 to n-1.
    /// 
    /// Note that subvector views give direct access to the underlying elements of the original vector. For example, the following code will
    /// zero the even elements of the vector v of length n, while leaving the odd elements untouched,
    /// 
    /// ```C
    /// gsl_vector_view v_even 
    ///   = gsl_vector_subvector_with_stride (v, 0, 2, n/2);
    /// gsl_vector_set_zero (&v_even.vector);
    /// ```
    /// A vector view can be passed to any subroutine which takes a vector argument just as a directly allocated vector would be, using &view.vector.
    /// For example, the following code computes the norm of the odd elements of v using the BLAS routine DNRM2,
    /// 
    /// ```C
    /// gsl_vector_view v_odd 
    ///   = gsl_vector_subvector_with_stride (v, 1, 2, n/2);
    /// double r = gsl_blas_dnrm2 (&v_odd.vector);
    /// ```
    /// The function gsl_vector_const_subvector_with_stride is equivalent to gsl_vector_subvector_with_stride but can be used for vectors which
    /// are declared const.
    pub fn from_vector_with_stride(v: &VectorF64, offset: usize, stride: usize, n: usize) -> VectorView {
        unsafe {
            VectorView {
                v: ffi::gsl_vector_subvector_with_stride(v.vec, offset, stride, n)
            }
        }
    }

    /// These functions return a vector view of an array. The start of the new vector is given by base and has n elements. Mathematically,
    /// the i-th element of the new vector v’ is given by,
    /// 
    /// v'(i) = base[i]
    /// 
    /// where the index i runs from 0 to n-1.
    /// 
    /// The array containing the elements of v is not owned by the new vector view. When the view goes out of scope the original array will
    /// continue to exist. The original memory can only be deallocated by freeing the original pointer base. Of course, the original array
    /// should not be deallocated while the view is still in use.
    /// 
    /// The function gsl_vector_const_view_array is equivalent to gsl_vector_view_array but can be used for arrays which are declared const.
    pub fn from_array(base: &mut [f64]) -> VectorView {
        unsafe {
            VectorView {
                v: ffi::gsl_vector_view_array(base.as_mut_ptr(), base.len() as usize)
            }
        }
    }

    /// These functions return a vector view of an array base with an additional stride argument. The subvector is formed in the same way as
    /// for gsl_vector_view_array but the new vector has n elements with a step-size of stride from one element to the next in the original
    /// array. Mathematically, the i-th element of the new vector v’ is given by,
    /// 
    /// v'(i) = base[i*stride]
    /// 
    /// where the index i runs from 0 to n-1.
    /// 
    /// Note that the view gives direct access to the underlying elements of the original array. A vector view can be passed to any subroutine
    /// which takes a vector argument just as a directly allocated vector would be, using &view.vector.
    /// 
    /// The function gsl_vector_const_view_array_with_stride is equivalent to gsl_vector_view_array_with_stride but can be used for arrays
    /// which are declared const.
    pub fn from_array_with_stride(base: &mut [f64], stride: usize) -> VectorView {
        unsafe {
            VectorView {
                v: ffi::gsl_vector_view_array_with_stride(base.as_mut_ptr(), stride, base.len() as usize)
            }
        }
    }

    pub fn vector(&mut self) -> VectorF64 {
        unsafe {
            VectorF64 {
                vec: ::std::mem::transmute(&mut self.v),
                can_free: false,
            }
        }
    }
}

pub struct VectorF32 {
    vec: *mut ffi::gsl_vector_float,
    can_free: bool,
}

impl VectorF32 {
    /// create a new VectorF32 with all elements set to zero
    pub fn new(size: usize) -> Option<VectorF32> {
        let tmp = unsafe { ffi::gsl_vector_float_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorF32 {
                vec: tmp,
                can_free: true,
            })
        }
    }

    pub fn from_slice(slice: &[f32]) -> Option<VectorF32> {
        let tmp = unsafe { ffi::gsl_vector_float_alloc(slice.len() as usize) };

        if tmp.is_null() {
            None
        } else {
            let v = VectorF32 {
                vec: tmp,
                can_free: true,
            };
            let mut pos = 0usize;

            for tmp in slice.iter() {
                v.set(pos, *tmp);
                pos += 1;
            }
            Some(v)
        }
    }

    pub fn len(&self) -> usize {
        if self.vec.is_null() {
            0usize
        } else {
            unsafe { (*self.vec).size }
        }
    }

    /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, i: usize) -> f32 {
        unsafe { ffi::gsl_vector_float_get(self.vec, i) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&self, i: usize, x: f32) -> &VectorF32 {
        unsafe { ffi::gsl_vector_float_set(self.vec, i, x) };
        self
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&self, x: f32) -> &VectorF32 {
        unsafe { ffi::gsl_vector_float_set_all(self.vec, x) };
        self
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&self) -> &VectorF32 {
        unsafe { ffi::gsl_vector_float_set_zero(self.vec) };
        self
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&self, i: usize) -> &VectorF32 {
        unsafe { ffi::gsl_vector_float_set_basis(self.vec, i) };
        self
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&self, other: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_memcpy(self.vec, other.vec) }
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_memcpy(other.vec, self.vec) }
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&self, other: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_swap(other.vec, self.vec) }
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&self, i: usize, j: usize) -> enums::Value {
        unsafe { ffi::gsl_vector_float_swap_elements(self.vec, i, j) }
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&self) -> enums::Value {
        unsafe { ffi::gsl_vector_float_reverse(self.vec) }
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&self, other: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_add(self.vec, other.vec) }
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&self, other: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_sub(self.vec, other.vec) }
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&self, other: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_mul(self.vec, other.vec) }
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&self, other: &VectorF32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_div(self.vec, other.vec) }
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&self, x: f32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_scale(self.vec, x) }
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&self, x: f32) -> enums::Value {
        unsafe { ffi::gsl_vector_float_add_constant(self.vec, x) }
    }

    /// This function returns the maximum value in the self vector.
    pub fn max(&self) -> f32 {
        unsafe { ffi::gsl_vector_float_max(self.vec) }
    }

    /// This function returns the minimum value in the self vector.
    pub fn min(&self) -> f32 {
        unsafe { ffi::gsl_vector_float_min(self.vec) }
    }

    /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
    pub fn minmax(&self) -> (f32, f32) {
        let mut min_out = 0.;
        let mut max_out = 0.;

        unsafe { ffi::gsl_vector_float_minmax(self.vec, &mut min_out, &mut max_out); }
        (min_out, max_out)
    }

    /// This function returns the index of the maximum value in the self vector.
    /// When there are several equal maximum elements then the lowest index is returned.
    pub fn max_index(&self) -> usize {
        unsafe { ffi::gsl_vector_float_max_index(self.vec) }
    }

    /// This function returns the index of the minimum value in the self vector.
    /// When there are several equal minimum elements then the lowest index is returned.
    pub fn min_index(&self) -> usize {
        unsafe { ffi::gsl_vector_float_min_index(self.vec) }
    }

    /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
    /// When there are several equal minimum or maximum elements then the lowest indices are returned.
    pub fn minmax_index(&self) -> (usize, usize) {
        let mut imin = 0usize;
        let mut imax = 0usize;

        unsafe { ffi::gsl_vector_float_minmax_index(self.vec, &mut imin, &mut imax) };
        (imin, imax)
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_isnull(self.vec) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_ispos(self.vec) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_isneg(self.vec) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_isnonneg(self.vec) } {
            1 => true,
            _ => false
        }
    }

    pub fn equal(&self, other: &VectorF32) -> bool {
        match unsafe { ffi::gsl_vector_float_equal(self.vec, other.vec) } {
            1 => true,
            _ => false
        }
    }

    // I'll find a way to do that later
    /*pub fn as_slice<'a>(&self) -> &'a [f32] {
        unsafe {
            if self.vec.is_null() {
                let tmp : Vec<f32> = Vec::new();

                tmp.as_ref()
            } else {
                let tmp : CSlice<f32> = CSlice::new((*self.vec).data, (*self.vec).size as usize);

                tmp.as_ref()
            }
        }
    }*/

    pub fn clone(&self) -> Option<VectorF32> {
        unsafe {
            if self.vec.is_null() {
                None
            } else {
                match VectorF32::new((*self.vec).size) {
                    Some(v) => {
                        v.copy_from(self);
                        Some(v)
                    }
                    None => None
                }
            }
        }
    }
}

impl Drop for VectorF32 {
    fn drop(&mut self) {
        if self.can_free {
            unsafe { ffi::gsl_vector_float_free(self.vec) };
            self.vec = ::std::ptr::null_mut();
        }
    }
}

impl Debug for VectorF32 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        unsafe {
            write!(f, "[");
            for x in 0usize..(*self.vec).size {
                if x < (*self.vec).size - 1 {
                    write!(f, "{}, ", self.get(x));
                } else {
                    write!(f, "{}", self.get(x));
                }
            }
        }
        write!(f, "]")
    }
}

impl ffi::FFI<ffi::gsl_vector_float> for VectorF32 {
    fn wrap(r: *mut ffi::gsl_vector_float) -> VectorF32 {
        VectorF32 {
            vec: r,
            can_free: true
        }
    }

    fn unwrap(v: &VectorF32) -> *mut ffi::gsl_vector_float {
        v.vec
    }
}

pub struct VectorF64 {
    vec: *mut ffi::gsl_vector,
    can_free: bool
}

impl VectorF64 {
    /// create a new VectorF64 with all elements set to zero
    pub fn new(size: usize) -> Option<VectorF64> {
        let tmp = unsafe { ffi::gsl_vector_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorF64 {
                vec: tmp,
                can_free: true
            })
        }
    }

    pub fn from_slice(slice: &[f64]) -> Option<VectorF64> {
        let tmp = unsafe { ffi::gsl_vector_alloc(slice.len() as usize) };

        if tmp.is_null() {
            None
        } else {
            let v = VectorF64 {
                vec: tmp,
                can_free: true
            };
            let mut pos = 0usize;

            for tmp in slice.iter() {
                v.set(pos, *tmp);
                pos += 1;
            }
            Some(v)
        }
    }

    pub fn len(&self) -> usize {
        if self.vec.is_null() {
            0usize
        } else {
            unsafe { (*self.vec).size }
        }
    }

    /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, i: usize) -> f64 {
        unsafe { ffi::gsl_vector_get(self.vec, i) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&self, i: usize, x: f64) -> &VectorF64 {
        unsafe { ffi::gsl_vector_set(self.vec, i, x) };
        self
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&self, x: f64) -> &VectorF64 {
        unsafe { ffi::gsl_vector_set_all(self.vec, x) };
        self
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&self) -> &VectorF64 {
        unsafe { ffi::gsl_vector_set_zero(self.vec) };
        self
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&self, i: usize) -> &VectorF64 {
        unsafe { ffi::gsl_vector_set_basis(self.vec, i) };
        self
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&self, other: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_vector_memcpy(self.vec, other.vec) }
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_vector_memcpy(other.vec, self.vec) }
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&self, other: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_vector_swap(other.vec, self.vec) }
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&self, i: usize, j: usize) -> enums::Value {
        unsafe { ffi::gsl_vector_swap_elements(self.vec, i, j) }
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&self) -> enums::Value {
        unsafe { ffi::gsl_vector_reverse(self.vec) }
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&self, other: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_vector_add(self.vec, other.vec) }
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&self, other: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_vector_sub(self.vec, other.vec) }
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&self, other: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_vector_mul(self.vec, other.vec) }
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&self, other: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_vector_div(self.vec, other.vec) }
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&self, x: f64) -> enums::Value {
        unsafe { ffi::gsl_vector_scale(self.vec, x) }
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&self, x: f64) -> enums::Value {
        unsafe { ffi::gsl_vector_add_constant(self.vec, x) }
    }

    /// This function returns the maximum value in the self vector.
    pub fn max(&self) -> f64 {
        unsafe { ffi::gsl_vector_max(self.vec) }
    }

    /// This function returns the minimum value in the self vector.
    pub fn min(&self) -> f64 {
        unsafe { ffi::gsl_vector_min(self.vec) }
    }

    /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
    pub fn minmax(&self) -> (f64, f64) {
        let mut min_out = 0.;
        let mut max_out = 0.;

        unsafe { ffi::gsl_vector_minmax(self.vec, &mut min_out, &mut max_out); }
        (min_out, max_out)
    }

    /// This function returns the index of the maximum value in the self vector.
    /// When there are several equal maximum elements then the lowest index is returned.
    pub fn max_index(&self) -> usize {
        unsafe { ffi::gsl_vector_max_index(self.vec) }
    }

    /// This function returns the index of the minimum value in the self vector.
    /// When there are several equal minimum elements then the lowest index is returned.
    pub fn min_index(&self) -> usize {
        unsafe { ffi::gsl_vector_min_index(self.vec) }
    }

    /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
    /// When there are several equal minimum or maximum elements then the lowest indices are returned.
    pub fn minmax_index(&self) -> (usize, usize) {
        let mut imin = 0usize;
        let mut imax = 0usize;

        unsafe { ffi::gsl_vector_minmax_index(self.vec, &mut imin, &mut imax) };
        (imin, imax)
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_vector_isnull(self.vec) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_vector_ispos(self.vec) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_isneg(self.vec) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_isnonneg(self.vec) } {
            1 => true,
            _ => false
        }
    }

    pub fn equal(&self, other: &VectorF64) -> bool {
        match unsafe { ffi::gsl_vector_equal(self.vec, other.vec) } {
            1 => true,
            _ => false
        }
    }

    // I'll find a way to do that later
    /*pub fn as_slice<'a>(&self) -> &'a [f64] {
        unsafe {
            if self.vec.is_null() {
                let tmp : Vec<f64> = Vec::new();

                tmp.as_ref()
            } else {
                let tmp : CSlice<f64> = CSlice::new((*self.vec).data, (*self.vec).size as usize);

                tmp.as_ref()
            }
        }
    }*/

    pub fn clone(&self) -> Option<VectorF64> {
        unsafe {
            if self.vec.is_null() {
                None
            } else {
                match VectorF64::new((*self.vec).size) {
                    Some(v) => {
                        v.copy_from(self);
                        Some(v)
                    }
                    None => None
                }
            }
        }
    }
}

impl Drop for VectorF64 {
    fn drop(&mut self) {
        if self.can_free {
            unsafe { ffi::gsl_vector_free(self.vec) };
            self.vec = ::std::ptr::null_mut();
        }
    }
}

impl Debug for VectorF64 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        unsafe {
            write!(f, "[");
            for x in 0usize..(*self.vec).size {
                if x < (*self.vec).size - 1 {
                    write!(f, "{}, ", self.get(x));
                } else {
                    write!(f, "{}", self.get(x));
                }
            }
        }
        write!(f, "]")
    }
}

impl ffi::FFI<ffi::gsl_vector> for VectorF64 {
    fn wrap(r: *mut ffi::gsl_vector) -> VectorF64 {
        VectorF64 {
            vec: r,
            can_free: true
        }
    }

    fn unwrap(v: &VectorF64) -> *mut ffi::gsl_vector {
        v.vec
    }
}

pub fn wrap(v: *mut ffi::gsl_vector) -> VectorF64 {
    VectorF64 {
        vec: v,
        can_free: false
    }
}