//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use enums;
use ffi::FFI;
use std::fmt;
use std::fmt::{Debug, Formatter};
use types::{ComplexF32, ComplexF64};

ffi_wrapper!(
    VectorComplexF64,
    *mut sys::gsl_vector_complex,
    gsl_vector_complex_free
);

impl VectorComplexF64 {
    /// create a new VectorComplexF64 with all elements set to zero
    pub fn new(size: usize) -> Option<VectorComplexF64> {
        let tmp = unsafe { sys::gsl_vector_complex_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorComplexF64::wrap(tmp))
        }
    }

    pub fn from_slice(slice: &[ComplexF64]) -> Option<VectorComplexF64> {
        let tmp = unsafe { sys::gsl_vector_complex_alloc(slice.len() as _) };

        if tmp.is_null() {
            None
        } else {
            let mut v = VectorComplexF64::wrap(tmp);

            for (pos, tmp) in slice.iter().enumerate() {
                v.set(pos as _, tmp);
            }
            Some(v)
        }
    }

    pub fn len(&self) -> usize {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            0
        } else {
            unsafe { (*ptr).size }
        }
    }

    /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, i: usize) -> ComplexF64 {
        unsafe { ::std::mem::transmute(sys::gsl_vector_complex_get(self.unwrap_shared(), i)) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&mut self, i: usize, x: &ComplexF64) -> &VectorComplexF64 {
        unsafe { sys::gsl_vector_complex_set(self.unwrap_unique(), i, ::std::mem::transmute(*x)) };
        self
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&mut self, x: &ComplexF64) -> &VectorComplexF64 {
        unsafe { sys::gsl_vector_complex_set_all(self.unwrap_unique(), ::std::mem::transmute(*x)) };
        self
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&mut self) -> &VectorComplexF64 {
        unsafe { sys::gsl_vector_complex_set_zero(self.unwrap_unique()) };
        self
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&mut self, i: usize) -> &VectorComplexF64 {
        unsafe { sys::gsl_vector_complex_set_basis(self.unwrap_unique(), i) };
        self
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&mut self, other: &VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_memcpy(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &mut VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_memcpy(other.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&mut self, other: &mut VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_swap(other.unwrap_unique(), self.unwrap_unique())
        })
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&mut self, i: usize, j: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_swap_elements(self.unwrap_unique(), i, j)
        })
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&mut self) -> enums::Value {
        enums::Value::from(unsafe { sys::gsl_vector_complex_reverse(self.unwrap_unique()) })
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&mut self, other: &VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_add(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&mut self, other: &VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_sub(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&mut self, other: &VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_mul(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&mut self, other: &VectorComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_div(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&mut self, x: &ComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_scale(self.unwrap_unique(), ::std::mem::transmute(*x))
        })
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&mut self, x: &ComplexF64) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_add_constant(self.unwrap_unique(), ::std::mem::transmute(*x))
        })
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_isnull(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_ispos(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_isneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_isnonneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    pub fn equal(&self, other: &VectorComplexF64) -> bool {
        match unsafe { sys::gsl_vector_complex_equal(self.unwrap_shared(), other.unwrap_shared()) }
        {
            1 => true,
            _ => false,
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

    pub fn clone(&self) -> Option<VectorComplexF64> {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            None
        } else {
            unsafe {
                match VectorComplexF64::new((*ptr).size) {
                    Some(mut v) => {
                        v.copy_from(self);
                        Some(v)
                    }
                    None => None,
                }
            }
        }
    }
}

impl Debug for VectorComplexF64 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            write!(f, "<null>")
        } else {
            unsafe {
                write!(f, "[");
                for x in 0..(*ptr).size {
                    if x < (*ptr).size - 1 {
                        write!(f, "{:?}, ", self.get(x));
                    } else {
                        write!(f, "{:?}", self.get(x));
                    }
                }
            }
            write!(f, "]")
        }
    }
}

ffi_wrapper!(
    VectorComplexF32,
    *mut sys::gsl_vector_complex_float,
    gsl_vector_complex_float_free
);

impl VectorComplexF32 {
    /// create a new VectorComplexF32 with all elements set to zero
    pub fn new(size: usize) -> Option<VectorComplexF32> {
        let tmp = unsafe { sys::gsl_vector_complex_float_calloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(VectorComplexF32::wrap(tmp))
        }
    }

    pub fn from_slice(slice: &[ComplexF32]) -> Option<VectorComplexF32> {
        let tmp = unsafe { sys::gsl_vector_complex_float_alloc(slice.len() as _) };

        if tmp.is_null() {
            None
        } else {
            let mut v = VectorComplexF32::wrap(tmp);

            for (pos, tmp) in slice.iter().enumerate() {
                v.set(pos as _, tmp);
            }
            Some(v)
        }
    }

    pub fn len(&self) -> usize {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            0
        } else {
            unsafe { (*ptr).size }
        }
    }

    /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, i: usize) -> ComplexF32 {
        unsafe { ::std::mem::transmute(sys::gsl_vector_complex_float_get(self.unwrap_shared(), i)) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
    pub fn set(&mut self, i: usize, x: &ComplexF32) -> &VectorComplexF32 {
        unsafe {
            sys::gsl_vector_complex_float_set(self.unwrap_unique(), i, ::std::mem::transmute(*x))
        };
        self
    }

    /// This function sets all the elements of the vector v to the value x.
    pub fn set_all(&mut self, x: &ComplexF32) -> &VectorComplexF32 {
        unsafe {
            sys::gsl_vector_complex_float_set_all(self.unwrap_unique(), ::std::mem::transmute(*x))
        };
        self
    }

    /// This function sets all the elements of the vector v to zero.
    pub fn set_zero(&mut self) -> &VectorComplexF32 {
        unsafe { sys::gsl_vector_complex_float_set_zero(self.unwrap_unique()) };
        self
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
    pub fn set_basis(&mut self, i: usize) -> &VectorComplexF32 {
        unsafe { sys::gsl_vector_complex_float_set_basis(self.unwrap_unique(), i) };
        self
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
    pub fn copy_from(&mut self, other: &VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_memcpy(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
    pub fn copy_to(&self, other: &mut VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_memcpy(other.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
    pub fn swap(&mut self, other: &mut VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_swap(other.unwrap_unique(), self.unwrap_unique())
        })
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    pub fn swap_elements(&mut self, i: usize, j: usize) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_swap_elements(self.unwrap_unique(), i, j)
        })
    }

    /// This function reverses the order of the elements of the vector v.
    pub fn reverse(&mut self) -> enums::Value {
        enums::Value::from(unsafe { sys::gsl_vector_complex_float_reverse(self.unwrap_unique()) })
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn add(&mut self, other: &VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_add(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function subtracts the elements of the self vector from the elements of the other vector.
    /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn sub(&mut self, other: &VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_sub(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self vector a by the elements of the other vector.
    /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn mul(&mut self, other: &VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_mul(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
    pub fn div(&mut self, other: &VectorComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_div(self.unwrap_unique(), other.unwrap_shared())
        })
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
    pub fn scale(&mut self, x: &ComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_scale(self.unwrap_unique(), ::std::mem::transmute(*x))
        })
    }

    /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
    pub fn add_constant(&mut self, x: &ComplexF32) -> enums::Value {
        enums::Value::from(unsafe {
            sys::gsl_vector_complex_float_add_constant(
                self.unwrap_unique(),
                ::std::mem::transmute(*x),
            )
        })
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    pub fn is_null(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_float_isnull(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    pub fn is_pos(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_float_ispos(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    pub fn is_neg(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_float_isneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    pub fn is_non_neg(&self) -> bool {
        match unsafe { sys::gsl_vector_complex_float_isnonneg(self.unwrap_shared()) } {
            1 => true,
            _ => false,
        }
    }

    pub fn equal(&self, other: &VectorComplexF32) -> bool {
        match unsafe {
            sys::gsl_vector_complex_float_equal(self.unwrap_shared(), other.unwrap_shared())
        } {
            1 => true,
            _ => false,
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

    pub fn clone(&self) -> Option<VectorComplexF32> {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            None
        } else {
            unsafe {
                match VectorComplexF32::new((*ptr).size) {
                    Some(mut v) => {
                        v.copy_from(self);
                        Some(v)
                    }
                    None => None,
                }
            }
        }
    }
}

impl Debug for VectorComplexF32 {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            write!(f, "<null>")
        } else {
            unsafe {
                write!(f, "[");
                for x in 0..(*ptr).size {
                    if x < (*ptr).size - 1 {
                        write!(f, "{:?}, ", self.get(x));
                    } else {
                        write!(f, "{:?}", self.get(x));
                    }
                }
            }
            write!(f, "]")
        }
    }
}
