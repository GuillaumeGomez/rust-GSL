//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::fmt;
use std::fmt::{Formatter,Show};
use types::{ComplexF32, ComplexF64};
use ffi;

pub struct VectorComplexF64 {
    vec: *mut ffi::gsl_vector_complex
}

impl VectorComplexF64 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_vector_complex {
        self.vec
    }

    /// create a new VectorComplexF64 with all elements set to zero
    pub fn new(size: u64) -> Option<VectorComplexF64> {
        let tmp = unsafe { ffi::gsl_vector_complex_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorComplexF64 {
                vec: tmp
            })
        }
    }

    pub fn from_slice(slice: &[ComplexF64]) -> Option<VectorComplexF64> {
        let tmp = unsafe { ffi::gsl_vector_complex_alloc(slice.len() as u64) };

        if tmp.is_null() {
            None
        } else {
            let v = VectorComplexF64 {
                vec: tmp
            };
            let mut pos = 0u64;

            for tmp in slice.iter() {
                v.set(pos, tmp);
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
    pub fn get(&self, i: u64) -> ComplexF64 {
        unsafe { ::std::mem::transmute(ffi::gsl_vector_complex_get(self.vec as *const ffi::gsl_vector_complex, i)) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&self, i: u64, x: &ComplexF64) {
        unsafe { ffi::gsl_vector_complex_set(self.vec, i, ::std::mem::transmute(*x)) }
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&self, x: &ComplexF64) {
        unsafe { ffi::gsl_vector_complex_set_all(self.vec, ::std::mem::transmute(*x)) }
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&self) {
        unsafe { ffi::gsl_vector_complex_set_zero(self.vec) }
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&self, i: u64) {
        unsafe { ffi::gsl_vector_complex_set_basis(self.vec, i) }
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&self, other: &VectorComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_memcpy(self.vec, other.vec as *const ffi::gsl_vector_complex) }
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &VectorComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_memcpy(other.vec, self.vec as *const ffi::gsl_vector_complex) }
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&self, other: &VectorComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_swap(other.vec, self.vec) }
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&self, i: u64, j: u64) -> i32 {
        unsafe { ffi::gsl_vector_complex_swap_elements(self.vec, i, j) }
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&self) -> i32 {
        unsafe { ffi::gsl_vector_complex_reverse(self.vec) }
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&self, other: &VectorComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_add(self.vec, other.vec as *const ffi::gsl_vector_complex) }
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&self, other: &VectorComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_sub(self.vec, other.vec as *const ffi::gsl_vector_complex) }
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&self, other: &VectorComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_mul(self.vec, other.vec as *const ffi::gsl_vector_complex) }
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&self, other: &VectorComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_div(self.vec, other.vec as *const ffi::gsl_vector_complex) }
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&self, x: &ComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_scale(self.vec, ::std::mem::transmute(*x)) }
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&self, x: &ComplexF64) -> i32 {
        unsafe { ffi::gsl_vector_complex_add_constant(self.vec, ::std::mem::transmute(*x)) }
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_isnull(self.vec as *const ffi::gsl_vector_complex) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_ispos(self.vec as *const ffi::gsl_vector_complex) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_isneg(self.vec as *const ffi::gsl_vector_complex) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_isnonneg(self.vec as *const ffi::gsl_vector_complex) } {
            1 => true,
            _ => false
        }
    }

    pub fn equal(&self, other: &VectorComplexF64) -> bool {
        match unsafe { ffi::gsl_vector_complex_equal(self.vec as *const ffi::gsl_vector_complex, other.vec as *const ffi::gsl_vector_complex) } {
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

    pub fn clone(&self) -> Option<VectorComplexF64> {
        unsafe {
            if self.vec.is_null() {
                None
            } else {
                match VectorComplexF64::new((*self.vec).size) {
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

impl Drop for VectorComplexF64 {
    fn drop(&mut self) {
        unsafe { ffi::gsl_vector_complex_free(self.vec) };
        self.vec = ::std::ptr::mut_null();
    }
}

impl Show for VectorComplexF64 {
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

impl ffi::FFI<ffi::gsl_vector_complex> for VectorComplexF64 {
    fn wrap(r: *mut ffi::gsl_vector_complex) -> VectorComplexF64 {
        VectorComplexF64 {
            vec: r
        }
    }

    fn unwrap(v: &VectorComplexF64) -> *mut ffi::gsl_vector_complex {
        v.vec
    }
}

pub struct VectorComplexF32 {
    vec: *mut ffi::gsl_vector_complex_float
}

impl VectorComplexF32 {
    #[doc(hidden)]
    #[allow(visible_private_types)]
    pub fn get_ffi(&self) -> *mut ffi::gsl_vector_complex_float {
        self.vec
    }

    /// create a new VectorComplexF32 with all elements set to zero
    pub fn new(size: u64) -> Option<VectorComplexF32> {
        let tmp = unsafe { ffi::gsl_vector_complex_float_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorComplexF32 {
                vec: tmp
            })
        }
    }

    pub fn from_slice(slice: &[ComplexF32]) -> Option<VectorComplexF32> {
        let tmp = unsafe { ffi::gsl_vector_complex_float_alloc(slice.len() as u64) };

        if tmp.is_null() {
            None
        } else {
            let v = VectorComplexF32 {
                vec: tmp
            };
            let mut pos = 0u64;

            for tmp in slice.iter() {
                v.set(pos, tmp);
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
    pub fn get(&self, i: u64) -> ComplexF32 {
        unsafe { ::std::mem::transmute(ffi::gsl_vector_complex_float_get(self.vec as *const ffi::gsl_vector_complex_float, i)) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&self, i: u64, x: &ComplexF32) {
        unsafe { ffi::gsl_vector_complex_float_set(self.vec, i, ::std::mem::transmute(*x)) }
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&self, x: &ComplexF32) {
        unsafe { ffi::gsl_vector_complex_float_set_all(self.vec, ::std::mem::transmute(*x)) }
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&self) {
        unsafe { ffi::gsl_vector_complex_float_set_zero(self.vec) }
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&self, i: u64) {
        unsafe { ffi::gsl_vector_complex_float_set_basis(self.vec, i) }
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&self, other: &VectorComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_memcpy(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &VectorComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_memcpy(other.vec, self.vec as *const ffi::gsl_vector_complex_float) }
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&self, other: &VectorComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_swap(other.vec, self.vec) }
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&self, i: u64, j: u64) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_swap_elements(self.vec, i, j) }
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&self) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_reverse(self.vec) }
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&self, other: &VectorComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_add(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&self, other: &VectorComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_sub(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&self, other: &VectorComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_mul(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&self, other: &VectorComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_div(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&self, x: &ComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_scale(self.vec, ::std::mem::transmute(*x)) }
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&self, x: &ComplexF32) -> i32 {
        unsafe { ffi::gsl_vector_complex_float_add_constant(self.vec, ::std::mem::transmute(*x)) }
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_float_isnull(self.vec as *const ffi::gsl_vector_complex_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_float_ispos(self.vec as *const ffi::gsl_vector_complex_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_float_isneg(self.vec as *const ffi::gsl_vector_complex_float) } {
            1 => true,
            _ => false
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { ffi::gsl_vector_complex_float_isnonneg(self.vec as *const ffi::gsl_vector_complex_float) } {
            1 => true,
            _ => false
        }
    }

    pub fn equal(&self, other: &VectorComplexF32) -> bool {
        match unsafe { ffi::gsl_vector_complex_float_equal(self.vec as *const ffi::gsl_vector_complex_float,
            other.vec as *const ffi::gsl_vector_complex_float) } {
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

    pub fn clone(&self) -> Option<VectorComplexF32> {
        unsafe {
            if self.vec.is_null() {
                None
            } else {
                match VectorComplexF32::new((*self.vec).size) {
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

impl Drop for VectorComplexF32 {
    fn drop(&mut self) {
        unsafe { ffi::gsl_vector_complex_float_free(self.vec) };
        self.vec = ::std::ptr::mut_null();
    }
}

impl Show for VectorComplexF32 {
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

impl ffi::FFI<ffi::gsl_vector_complex_float> for VectorComplexF32 {
    fn wrap(r: *mut ffi::gsl_vector_complex_float) -> VectorComplexF32 {
        VectorComplexF32 {
            vec: r
        }
    }

    fn unwrap(v: &VectorComplexF32) -> *mut ffi::gsl_vector_complex_float {
        v.vec
    }
}