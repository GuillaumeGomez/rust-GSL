//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::fmt;
use std::fmt::{Formatter,Show};
use ffi;

pub struct VectorF32 {
    vec: *mut ffi::gsl_vector_float
}

impl VectorF32 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_vector_float {
        self.vec
    }

    /// create a new VectorF32 with all elements set to zero
    pub fn new(size: u64) -> Option<VectorF32> {
        let tmp = unsafe { ffi::gsl_vector_float_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorF32 {
                vec: tmp
            })
        }
    }

    pub fn from_slice(slice: &[f32]) -> Option<VectorF32> {
        let tmp = unsafe { ffi::gsl_vector_float_alloc(slice.len() as u64) };

        if tmp.is_null() {
            None
        } else {
            let v = VectorF32 {
                vec: tmp
            };
            let mut pos = 0u64;

            for tmp in slice.iter() {
                v.set(pos, *tmp);
                pos += 1;
            }
            Some(v)
        }
    }

    pub fn len(&self) -> u64 {
        if self.vec.is_null() {
            0u64
        } else {
            unsafe { (*self.vec).size }
        }
    }

    /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, i: u64) -> f32 {
        unsafe { ffi::gsl_vector_float_get(self.vec as *const ffi::gsl_vector_float, i) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&self, i: u64, x: f32) {
        unsafe { ffi::gsl_vector_float_set(self.vec, i, x) }
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&self, x: f32) {
        unsafe { ffi::gsl_vector_float_set_all(self.vec, x) }
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&self) {
        unsafe { ffi::gsl_vector_float_set_zero(self.vec) }
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&self, i: u64) {
        unsafe { ffi::gsl_vector_float_set_basis(self.vec, i) }
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&self, other: &VectorF32) -> i32 {
        unsafe { ffi::gsl_vector_float_memcpy(self.vec, other.vec as *const ffi::gsl_vector_float) }
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &VectorF32) -> i32 {
        unsafe { ffi::gsl_vector_float_memcpy(other.vec, self.vec as *const ffi::gsl_vector_float) }
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&self, other: &VectorF32) -> i32 {
        unsafe { ffi::gsl_vector_float_swap(other.vec, self.vec) }
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&self, i: u64, j: u64) -> i32 {
        unsafe { ffi::gsl_vector_float_swap_elements(self.vec, i, j) }
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&self) -> i32 {
        unsafe { ffi::gsl_vector_float_reverse(self.vec) }
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&self, other: &VectorF32) -> i32 {
        unsafe { ffi::gsl_vector_float_add(self.vec, other.vec as *const ffi::gsl_vector_float) }
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&self, other: &VectorF32) -> i32 {
        unsafe { ffi::gsl_vector_float_sub(self.vec, other.vec as *const ffi::gsl_vector_float) }
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&self, other: &VectorF32) -> i32 {
        unsafe { ffi::gsl_vector_float_mul(self.vec, other.vec as *const ffi::gsl_vector_float) }
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&self, other: &VectorF32) -> i32 {
        unsafe { ffi::gsl_vector_float_div(self.vec, other.vec as *const ffi::gsl_vector_float) }
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&self, x: f32) -> i32 {
        unsafe { ffi::gsl_vector_float_scale(self.vec, x) }
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&self, x: f32) -> i32 {
        unsafe { ffi::gsl_vector_float_add_constant(self.vec, x) }
    }

    /// This function returns the maximum value in the self vector.
    pub fn max(&self) -> f32 {
        unsafe { ffi::gsl_vector_float_max(self.vec as *const ffi::gsl_vector_float) }
    }

    /// This function returns the minimum value in the self vector.
    pub fn min(&self) -> f32 {
        unsafe { ffi::gsl_vector_float_min(self.vec as *const ffi::gsl_vector_float) }
    }

    /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
    pub fn minmax(&self, min_out: &mut f32, max_out: &mut f32) {
        unsafe { ffi::gsl_vector_float_minmax(self.vec as *const ffi::gsl_vector_float, min_out, max_out) }
    }

    /// This function returns the index of the maximum value in the self vector.
    /// When there are several equal maximum elements then the lowest index is returned.
    pub fn max_index(&self) -> u64 {
        unsafe { ffi::gsl_vector_float_max_index(self.vec as *const ffi::gsl_vector_float) }
    }

    /// This function returns the index of the minimum value in the self vector.
    /// When there are several equal minimum elements then the lowest index is returned.
    pub fn min_index(&self) -> u64 {
        unsafe { ffi::gsl_vector_float_min_index(self.vec as *const ffi::gsl_vector_float) }
    }

    /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
    /// When there are several equal minimum or maximum elements then the lowest indices are returned.
    pub fn minmax_index(&self) -> (u64, u64) {
        let mut imin = 0u64;
        let mut imax = 0u64;

        unsafe { ffi::gsl_vector_float_minmax_index(self.vec as *const ffi::gsl_vector_float, &mut imin, &mut imax) };
        (imin, imax)
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_isnull(self.vec as *const ffi::gsl_vector_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_ispos(self.vec as *const ffi::gsl_vector_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_isneg(self.vec as *const ffi::gsl_vector_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_float_isnonneg(self.vec as *const ffi::gsl_vector_float) } {
            1 => true,
            _ => false
        }
    }

    pub fn equal(&self, other: &VectorF32) -> bool {
        match unsafe { ffi::gsl_vector_float_equal(self.vec as *const ffi::gsl_vector_float, other.vec as *const ffi::gsl_vector_float) } {
            1 => true,
            _ => false
        }
    }

    // I'll find a way to do that later
    /*pub fn as_slice<'a>(&self) -> &'a [f32] {
        unsafe {
            if self.vec.is_null() {
                let tmp : Vec<f32> = Vec::new();

                tmp.as_slice()
            } else {
                let tmp : CVec<f32> = CVec::new((*self.vec).data, (*self.vec).size as uint);

                tmp.as_slice()
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
        unsafe { ffi::gsl_vector_float_free(self.vec) };
        self.vec = ::std::ptr::mut_null();
    }
}

impl Show for VectorF32 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        unsafe {
            write!(f, "[");
            for x in range(0u64, (*self.vec).size) {
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
            vec: r
        }
    }

    fn unwrap(v: &VectorF32) -> *mut ffi::gsl_vector_float {
        v.vec
    }
}

pub struct VectorF64 {
    vec: *mut ffi::gsl_vector
}

impl VectorF64 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_vector {
        self.vec
    }

    /// create a new VectorF64 with all elements set to zero
    pub fn new(size: u64) -> Option<VectorF64> {
        let tmp = unsafe { ffi::gsl_vector_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorF64 {
                vec: tmp
            })
        }
    }

    pub fn from_slice(slice: &[f64]) -> Option<VectorF64> {
        let tmp = unsafe { ffi::gsl_vector_alloc(slice.len() as u64) };

        if tmp.is_null() {
            None
        } else {
            let v = VectorF64 {
                vec: tmp
            };
            let mut pos = 0u64;

            for tmp in slice.iter() {
                v.set(pos, *tmp);
                pos += 1;
            }
            Some(v)
        }
    }

    pub fn len(&self) -> u64 {
        if self.vec.is_null() {
            0u64
        } else {
            unsafe { (*self.vec).size }
        }
    }

    /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, i: u64) -> f64 {
        unsafe { ffi::gsl_vector_get(self.vec as *const ffi::gsl_vector, i) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&self, i: u64, x: f64) {
        unsafe { ffi::gsl_vector_set(self.vec, i, x) }
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&self, x: f64) {
        unsafe { ffi::gsl_vector_set_all(self.vec, x) }
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&self) {
        unsafe { ffi::gsl_vector_set_zero(self.vec) }
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&self, i: u64) {
        unsafe { ffi::gsl_vector_set_basis(self.vec, i) }
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&self, other: &VectorF64) -> i32 {
        unsafe { ffi::gsl_vector_memcpy(self.vec, other.vec as *const ffi::gsl_vector) }
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &VectorF64) -> i32 {
        unsafe { ffi::gsl_vector_memcpy(other.vec, self.vec as *const ffi::gsl_vector) }
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&self, other: &VectorF64) -> i32 {
        unsafe { ffi::gsl_vector_swap(other.vec, self.vec) }
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&self, i: u64, j: u64) -> i32 {
        unsafe { ffi::gsl_vector_swap_elements(self.vec, i, j) }
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&self) -> i32 {
        unsafe { ffi::gsl_vector_reverse(self.vec) }
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&self, other: &VectorF64) -> i32 {
        unsafe { ffi::gsl_vector_add(self.vec, other.vec as *const ffi::gsl_vector) }
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&self, other: &VectorF64) -> i32 {
        unsafe { ffi::gsl_vector_sub(self.vec, other.vec as *const ffi::gsl_vector) }
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&self, other: &VectorF64) -> i32 {
        unsafe { ffi::gsl_vector_mul(self.vec, other.vec as *const ffi::gsl_vector) }
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&self, other: &VectorF64) -> i32 {
        unsafe { ffi::gsl_vector_div(self.vec, other.vec as *const ffi::gsl_vector) }
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&self, x: f64) -> i32 {
        unsafe { ffi::gsl_vector_scale(self.vec, x) }
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&self, x: f64) -> i32 {
        unsafe { ffi::gsl_vector_add_constant(self.vec, x) }
    }

    /// This function returns the maximum value in the self vector.
    pub fn max(&self) -> f64 {
        unsafe { ffi::gsl_vector_max(self.vec as *const ffi::gsl_vector) }
    }

    /// This function returns the minimum value in the self vector.
    pub fn min(&self) -> f64 {
        unsafe { ffi::gsl_vector_min(self.vec as *const ffi::gsl_vector) }
    }

    /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
    pub fn minmax(&self, min_out: &mut f64, max_out: &mut f64) {
        unsafe { ffi::gsl_vector_minmax(self.vec as *const ffi::gsl_vector, min_out, max_out) }
    }

    /// This function returns the index of the maximum value in the self vector.
    /// When there are several equal maximum elements then the lowest index is returned.
    pub fn max_index(&self) -> u64 {
        unsafe { ffi::gsl_vector_max_index(self.vec as *const ffi::gsl_vector) }
    }

    /// This function returns the index of the minimum value in the self vector.
    /// When there are several equal minimum elements then the lowest index is returned.
    pub fn min_index(&self) -> u64 {
        unsafe { ffi::gsl_vector_min_index(self.vec as *const ffi::gsl_vector) }
    }

    /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
    /// When there are several equal minimum or maximum elements then the lowest indices are returned.
    pub fn minmax_index(&self) -> (u64, u64) {
        let mut imin = 0u64;
        let mut imax = 0u64;

        unsafe { ffi::gsl_vector_minmax_index(self.vec as *const ffi::gsl_vector, &mut imin, &mut imax) };
        (imin, imax)
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_vector_isnull(self.vec as *const ffi::gsl_vector) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_vector_ispos(self.vec as *const ffi::gsl_vector) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_isneg(self.vec as *const ffi::gsl_vector) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_isnonneg(self.vec as *const ffi::gsl_vector) } {
            1 => true,
            _ => false
        }
    }

    pub fn equal(&self, other: &VectorF64) -> bool {
        match unsafe { ffi::gsl_vector_equal(self.vec as *const ffi::gsl_vector, other.vec as *const ffi::gsl_vector) } {
            1 => true,
            _ => false
        }
    }

    // I'll find a way to do that later
    /*pub fn as_slice<'a>(&self) -> &'a [f64] {
        unsafe {
            if self.vec.is_null() {
                let tmp : Vec<f64> = Vec::new();

                tmp.as_slice()
            } else {
                let tmp : CVec<f64> = CVec::new((*self.vec).data, (*self.vec).size as uint);

                tmp.as_slice()
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
        unsafe { ffi::gsl_vector_free(self.vec) };
        self.vec = ::std::ptr::mut_null();
    }
}

impl Show for VectorF64 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        unsafe {
            write!(f, "[");
            for x in range(0u64, (*self.vec).size) {
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
            vec: r
        }
    }

    fn unwrap(v: &VectorF64) -> *mut ffi::gsl_vector {
        v.vec
    }
}