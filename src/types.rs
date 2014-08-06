//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub type gsl_mode_t = u32;
pub struct CblasIndex(pub u32);

pub mod Gsl {
    use std::default::Default;
    use ffi;
    use std::fmt::{Formatter,Show};
    use std::fmt;

    /// The error handling form of the special functions always calculate an error estimate along with the value of the result. Therefore, structures are provided for amalgamating a value and error estimate.
    pub struct Result {
        /// Contains the value.
        pub val: f64,
        /// Contains an estimate of the absolute error in the value.
        pub err: f64
    }

    impl Result {
        pub fn new() -> Result {
            Result {
                val: 0f64,
                err: 0f64
            }
        }
    }

    impl Default for Result {
        fn default() -> Result {
            Result::new()
        }
    }

    pub struct Complex {
        pub data: [f64, ..2]
    }

    impl Show for Complex {
        #[allow(unused_must_use)]
        fn fmt(&self, f: &mut Formatter) -> fmt::Result {
            write!(f, "[{}, {}]", self.data[0], self.data[1])
        }
    }

    impl Default for Complex {
        fn default() -> Complex {
            Complex {
                data: [0f64, 0f64]
            }
        }
    }

    pub struct ComplexFloat {
        pub data: [f32, ..2]
    }

    impl Show for ComplexFloat {
        #[allow(unused_must_use)]
        fn fmt(&self, f: &mut Formatter) -> fmt::Result {
            write!(f, "[{}, {}]", self.data[0], self.data[1])
        }
    }

    impl Default for ComplexFloat {
        fn default() -> ComplexFloat {
            ComplexFloat {
                data: [0f32, 0f32]
            }
        }
    }

    pub struct VectorFloat {
        vec: *mut ffi::gsl_vector_float
    }

    impl VectorFloat {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_vector_float {
            self.vec
        }

        /// create a new VectorFloat with all elements set to zero
        pub fn new(size: u64) -> Option<VectorFloat> {
            let tmp = unsafe { ffi::gsl_vector_float_calloc(size) };

            if tmp.is_null() {
                None
            } else {
                Some(VectorFloat {
                    vec: tmp
                })
            }
        }

        pub fn from_slice(slice: &[f32]) -> Option<VectorFloat> {
            let tmp = unsafe { ffi::gsl_vector_float_alloc(slice.len() as u64) };

            if tmp.is_null() {
                None
            } else {
                let v = VectorFloat {
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
        pub fn copy_from(&self, other: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_vector_float_memcpy(self.vec, other.vec as *const ffi::gsl_vector_float) }
        }

        /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
        pub fn copy_to(&self, other: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_vector_float_memcpy(other.vec, self.vec as *const ffi::gsl_vector_float) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &VectorFloat) -> i32 {
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
        pub fn add(&self, other: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_vector_float_add(self.vec, other.vec as *const ffi::gsl_vector_float) }
        }

        /// This function subtracts the elements of the self vector from the elements of the other vector.
        /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn sub(&self, other: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_vector_float_sub(self.vec, other.vec as *const ffi::gsl_vector_float) }
        }

        /// This function multiplies the elements of the self vector a by the elements of the other vector.
        /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn mul(&self, other: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_vector_float_mul(self.vec, other.vec as *const ffi::gsl_vector_float) }
        }

        /// This function divides the elements of the self vector by the elements of the other vector.
        /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn div(&self, other: &VectorFloat) -> i32 {
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

        pub fn equal(&self, other: &VectorFloat) -> bool {
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

        pub fn clone(&self) -> Option<VectorFloat> {
            unsafe {
                if self.vec.is_null() {
                    None
                } else {
                    match VectorFloat::new((*self.vec).size) {
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

    impl Drop for VectorFloat {
        fn drop(&mut self) {
            unsafe { ffi::gsl_vector_float_free(self.vec) };
            self.vec = ::std::ptr::mut_null();
        }
    }

    impl Show for VectorFloat {
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

    pub struct Vector {
        vec: *mut ffi::gsl_vector
    }

    impl Vector {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_vector {
            self.vec
        }

        /// create a new Vector with all elements set to zero
        pub fn new(size: u64) -> Option<Vector> {
            let tmp = unsafe { ffi::gsl_vector_calloc(size) };

            if tmp.is_null() {
                None
            } else {
                Some(Vector {
                    vec: tmp
                })
            }
        }

        pub fn from_slice(slice: &[f64]) -> Option<Vector> {
            let tmp = unsafe { ffi::gsl_vector_alloc(slice.len() as u64) };

            if tmp.is_null() {
                None
            } else {
                let v = Vector {
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
        pub fn copy_from(&self, other: &Vector) -> i32 {
            unsafe { ffi::gsl_vector_memcpy(self.vec, other.vec as *const ffi::gsl_vector) }
        }

        /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
        pub fn copy_to(&self, other: &Vector) -> i32 {
            unsafe { ffi::gsl_vector_memcpy(other.vec, self.vec as *const ffi::gsl_vector) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &Vector) -> i32 {
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
        pub fn add(&self, other: &Vector) -> i32 {
            unsafe { ffi::gsl_vector_add(self.vec, other.vec as *const ffi::gsl_vector) }
        }

        /// This function subtracts the elements of the self vector from the elements of the other vector.
        /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn sub(&self, other: &Vector) -> i32 {
            unsafe { ffi::gsl_vector_sub(self.vec, other.vec as *const ffi::gsl_vector) }
        }

        /// This function multiplies the elements of the self vector a by the elements of the other vector.
        /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn mul(&self, other: &Vector) -> i32 {
            unsafe { ffi::gsl_vector_mul(self.vec, other.vec as *const ffi::gsl_vector) }
        }

        /// This function divides the elements of the self vector by the elements of the other vector.
        /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn div(&self, other: &Vector) -> i32 {
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

        pub fn equal(&self, other: &Vector) -> bool {
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

        pub fn clone(&self) -> Option<Vector> {
            unsafe {
                if self.vec.is_null() {
                    None
                } else {
                    match Vector::new((*self.vec).size) {
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

    impl Drop for Vector {
        fn drop(&mut self) {
            unsafe { ffi::gsl_vector_free(self.vec) };
            self.vec = ::std::ptr::mut_null();
        }
    }

    impl Show for Vector {
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

    pub struct VectorComplex {
        vec: *mut ffi::gsl_vector_complex
    }

    impl VectorComplex {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_vector_complex {
            self.vec
        }

        /// create a new VectorComplex with all elements set to zero
        pub fn new(size: u64) -> Option<VectorComplex> {
            let tmp = unsafe { ffi::gsl_vector_complex_calloc(size) };

            if tmp.is_null() {
                None
            } else {
                Some(VectorComplex {
                    vec: tmp
                })
            }
        }

        pub fn from_slice(slice: &[Complex]) -> Option<VectorComplex> {
            let tmp = unsafe { ffi::gsl_vector_complex_alloc(slice.len() as u64) };

            if tmp.is_null() {
                None
            } else {
                let v = VectorComplex {
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
        pub fn get(&self, i: u64) -> Complex {
            unsafe { ::std::mem::transmute(ffi::gsl_vector_complex_get(self.vec as *const ffi::gsl_vector_complex, i)) }
        }

        /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
        pub fn set(&self, i: u64, x: &Complex) {
            unsafe { ffi::gsl_vector_complex_set(self.vec, i, ::std::mem::transmute(*x)) }
        }

        /// This function sets all the elements of the vector v to the value x.
        pub fn set_all(&self, x: &Complex) {
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
        pub fn copy_from(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_memcpy(self.vec, other.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
        pub fn copy_to(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_memcpy(other.vec, self.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &VectorComplex) -> i32 {
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
        pub fn add(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_add(self.vec, other.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function subtracts the elements of the self vector from the elements of the other vector.
        /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn sub(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_sub(self.vec, other.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function multiplies the elements of the self vector a by the elements of the other vector.
        /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn mul(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_mul(self.vec, other.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function divides the elements of the self vector by the elements of the other vector.
        /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn div(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_div(self.vec, other.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
        pub fn scale(&self, x: &Complex) -> i32 {
            unsafe { ffi::gsl_vector_complex_scale(self.vec, ::std::mem::transmute(*x)) }
        }

        /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
        pub fn add_constant(&self, x: &Complex) -> i32 {
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

        pub fn equal(&self, other: &VectorComplex) -> bool {
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

        pub fn clone(&self) -> Option<VectorComplex> {
            unsafe {
                if self.vec.is_null() {
                    None
                } else {
                    match VectorComplex::new((*self.vec).size) {
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

    impl Drop for VectorComplex {
        fn drop(&mut self) {
            unsafe { ffi::gsl_vector_complex_free(self.vec) };
            self.vec = ::std::ptr::mut_null();
        }
    }

    impl Show for VectorComplex {
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

    pub struct VectorComplexFloat {
        vec: *mut ffi::gsl_vector_complex_float
    }

    impl VectorComplexFloat {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_vector_complex_float {
            self.vec
        }

        /// create a new VectorComplexFloat with all elements set to zero
        pub fn new(size: u64) -> Option<VectorComplexFloat> {
            let tmp = unsafe { ffi::gsl_vector_complex_float_calloc(size) };

            if tmp.is_null() {
                None
            } else {
                Some(VectorComplexFloat {
                    vec: tmp
                })
            }
        }

        pub fn from_slice(slice: &[ComplexFloat]) -> Option<VectorComplexFloat> {
            let tmp = unsafe { ffi::gsl_vector_complex_float_alloc(slice.len() as u64) };

            if tmp.is_null() {
                None
            } else {
                let v = VectorComplexFloat {
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
        pub fn get(&self, i: u64) -> ComplexFloat {
            unsafe { ::std::mem::transmute(ffi::gsl_vector_complex_float_get(self.vec as *const ffi::gsl_vector_complex_float, i)) }
        }

        /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
        pub fn set(&self, i: u64, x: &ComplexFloat) {
            unsafe { ffi::gsl_vector_complex_float_set(self.vec, i, ::std::mem::transmute(*x)) }
        }

        /// This function sets all the elements of the vector v to the value x.
        pub fn set_all(&self, x: &ComplexFloat) {
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
        pub fn copy_from(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_memcpy(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function copies the elements of the self vector into the other vector. The two vectors must have the same length.
        pub fn copy_to(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_memcpy(other.vec, self.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &VectorComplexFloat) -> i32 {
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
        pub fn add(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_add(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function subtracts the elements of the self vector from the elements of the other vector.
        /// The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn sub(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_sub(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function multiplies the elements of the self vector a by the elements of the other vector.
        /// The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn mul(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_mul(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function divides the elements of the self vector by the elements of the other vector.
        /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors must have the same length.
        pub fn div(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_div(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function multiplies the elements of the self vector by the constant factor x. The result a_i <- a_i is stored in self.
        pub fn scale(&self, x: &ComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_scale(self.vec, ::std::mem::transmute(*x)) }
        }

        /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
        pub fn add_constant(&self, x: &ComplexFloat) -> i32 {
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

        pub fn equal(&self, other: &VectorComplexFloat) -> bool {
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

        pub fn clone(&self) -> Option<VectorComplexFloat> {
            unsafe {
                if self.vec.is_null() {
                    None
                } else {
                    match VectorComplexFloat::new((*self.vec).size) {
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

    impl Drop for VectorComplexFloat {
        fn drop(&mut self) {
            unsafe { ffi::gsl_vector_complex_float_free(self.vec) };
            self.vec = ::std::ptr::mut_null();
        }
    }

    impl Show for VectorComplexFloat {
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

    pub struct Matrix {
        mat: *mut ffi::gsl_matrix
    }

    impl Matrix {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_matrix {
            self.mat
        }

        /// Creates a new Matrix with all elements set to zero
        /// 
        /// Example with n1 = 2 and n2 = 3 :
        /// 
        /// XX XX XX
        /// 
        /// XX XX XX
        pub fn new(n1: u64, n2: u64) -> Option<Matrix> {
            let tmp = unsafe { ffi::gsl_matrix_calloc(n1, n2) };

            if tmp.is_null() {
                None
            } else {
                Some(Matrix {
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
        pub fn set(&self, y: u64, x: u64, value: f64) {
            unsafe { ffi::gsl_matrix_set(self.mat, y, x, value) }
        }

        /// This function sets all the elements of the matrix to the value x.
        pub fn set_all(&self, x: f64) {
            unsafe { ffi::gsl_matrix_set_all(self.mat, x) }
        }

        /// This function sets all the elements of the matrix to zero.
        pub fn set_zero(&self) {
            unsafe { ffi::gsl_matrix_set_zero(self.mat) }
        }

        /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
        /// This applies to both square and rectangular matrices.
        pub fn set_identity(&self) {
            unsafe { ffi::gsl_matrix_set_identity(self.mat) }
        }

        /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
        pub fn copy_from(&self, other: &Matrix) -> i32 {
            unsafe { ffi::gsl_matrix_memcpy(self.mat, other.mat as *const ffi::gsl_matrix) }
        }

        /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
        pub fn copy_to(&self, other: &Matrix) -> i32 {
            unsafe { ffi::gsl_matrix_memcpy(other.mat, self.mat as *const ffi::gsl_matrix) }
        }

        /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
        pub fn swap(&self, other: &Matrix) -> i32 {
            unsafe { ffi::gsl_matrix_swap(self.mat, other.mat) }
        }

        /// This function copies the elements of the y-th row of the matrix into the returned vector.
        pub fn get_row(&self, y: u64) -> Option<(Vector, i32)> {
            let tmp = unsafe { ffi::gsl_vector_alloc((*self.mat).size2) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_get_row(tmp, self.mat as *const ffi::gsl_matrix, y) };

                Some((Vector{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the x-th column of the matrix into the returned vector.
        pub fn get_col(&self, x: u64) -> Option<(Vector, i32)> {
            let tmp = unsafe { ffi::gsl_vector_alloc((*self.mat).size1) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_get_col(tmp, self.mat as *const ffi::gsl_matrix, x) };

                Some((Vector{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the vector v into the y-th row of the matrix.
        /// The length of the vector must be the same as the length of the row.
        pub fn set_row(&self, y: u64, v: &Vector) -> i32 {
            unsafe { ffi::gsl_matrix_set_row(self.mat, y, v.vec as *const ffi::gsl_vector) }
        }

        /// This function copies the elements of the vector v into the x-th column of the matrix.
        /// The length of the vector must be the same as the length of the column.
        pub fn set_col(&self, x: u64, v: &Vector) -> i32 {
            unsafe { ffi::gsl_matrix_set_col(self.mat, x, v.vec as *const ffi::gsl_vector) }
        }

        /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
        pub fn swap_rows(&self, y1: u64, y2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_swap_rows(self.mat, y1, y2) }
        }

        /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
        pub fn swap_columns(&self, x1: u64, x2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_swap_columns(self.mat, x1, x2) }
        }

        /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
        pub fn swap_row_col(&self, i: u64, j: u64) -> i32 {
            unsafe { ffi::gsl_matrix_swap_rowcol(self.mat, i, j) }
        }

        /// This function returns the transpose of the matrix by copying the elements into it.
        /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
        pub fn transpose_memcpy(&self) -> Option<(Matrix, i32)> {
            let dest = unsafe { ffi::gsl_matrix_alloc((*self.mat).size1, (*self.mat).size2) };

            if dest.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix) };

                Some((Matrix {mat: dest}, ret))
            }
        }

        /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
        /// The matrix must be square for this operation to be possible.
        pub fn transpose(&self) -> i32 {
            unsafe { ffi::gsl_matrix_transpose(self.mat) }
        }

        /// This function adds the elements of the other matrix to the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn add(&self, other: &Matrix) -> i32 {
            unsafe { ffi::gsl_matrix_add(self.mat, other.mat as *const ffi::gsl_matrix) }
        }

        /// This function subtracts the elements of the other matrix from the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn sub(&self, other: &Matrix) -> i32 {
            unsafe { ffi::gsl_matrix_sub(self.mat, other.mat as *const ffi::gsl_matrix) }
        }

        /// This function multiplies the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn mul_elements(&self, other: &Matrix) -> i32 {
            unsafe { ffi::gsl_matrix_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix) }
        }

        /// This function divides the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn div_elements(&self, other: &Matrix) -> i32 {
            unsafe { ffi::gsl_matrix_div_elements(self.mat, other.mat as *const ffi::gsl_matrix) }
        }

        /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
        pub fn scale(&self, x: f64) -> i32 {
            unsafe { ffi::gsl_matrix_scale(self.mat, x) }
        }

        /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
        pub fn add_constant(&self, x: f64) -> i32 {
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
        pub fn equal(&self, other: &Matrix) -> bool {
            match unsafe { ffi::gsl_matrix_equal(self.mat as *const ffi::gsl_matrix, other.mat as *const ffi::gsl_matrix) } {
                1 => true,
                _ => false
            }
        }

        pub fn clone(&self) -> Option<Matrix> {
            unsafe {
                if self.mat.is_null() {
                    None
                } else {
                    match Matrix::new((*self.mat).size1, (*self.mat).size2) {
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

    impl Drop for Matrix {
        fn drop(&mut self) {
            unsafe { ffi::gsl_matrix_free(self.mat) };
            self.mat = ::std::ptr::mut_null();
        }
    }

    impl Show for Matrix {
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

    pub struct MatrixFloat {
        mat: *mut ffi::gsl_matrix_float
    }

    impl MatrixFloat {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_matrix_float {
            self.mat
        }

        /// Creates a new Matrix with all elements set to zero
        /// 
        /// Example with n1 = 2 and n2 = 3 :
        /// 
        /// XX XX XX
        /// 
        /// XX XX XX
        pub fn new(n1: u64, n2: u64) -> Option<MatrixFloat> {
            let tmp = unsafe { ffi::gsl_matrix_float_calloc(n1, n2) };

            if tmp.is_null() {
                None
            } else {
                Some(MatrixFloat {
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
        pub fn set(&self, y: u64, x: u64, value: f32) {
            unsafe { ffi::gsl_matrix_float_set(self.mat, y, x, value) }
        }

        /// This function sets all the elements of the matrix to the value x.
        pub fn set_all(&self, x: f32) {
            unsafe { ffi::gsl_matrix_float_set_all(self.mat, x) }
        }

        /// This function sets all the elements of the matrix to zero.
        pub fn set_zero(&self) {
            unsafe { ffi::gsl_matrix_float_set_zero(self.mat) }
        }

        /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
        /// This applies to both square and rectangular matrices.
        pub fn set_identity(&self) {
            unsafe { ffi::gsl_matrix_float_set_identity(self.mat) }
        }

        /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
        pub fn copy_from(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_memcpy(self.mat, other.mat as *const ffi::gsl_matrix_float) }
        }

        /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
        pub fn copy_to(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_memcpy(other.mat, self.mat as *const ffi::gsl_matrix_float) }
        }

        /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
        pub fn swap(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_swap(self.mat, other.mat) }
        }

        /// This function copies the elements of the y-th row of the matrix into the returned vector.
        pub fn get_row(&self, y: u64) -> Option<(VectorFloat, i32)> {
            let tmp = unsafe { ffi::gsl_vector_float_alloc((*self.mat).size2) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_float_get_row(tmp, self.mat as *const ffi::gsl_matrix_float, y) };

                Some((VectorFloat{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the x-th column of the matrix into the returned vector.
        pub fn get_col(&self, x: u64) -> Option<(VectorFloat, i32)> {
            let tmp = unsafe { ffi::gsl_vector_float_alloc((*self.mat).size1) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_float_get_col(tmp, self.mat as *const ffi::gsl_matrix_float, x) };

                Some((VectorFloat{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the vector v into the y-th row of the matrix.
        /// The length of the vector must be the same as the length of the row.
        pub fn set_row(&self, y: u64, v: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_set_row(self.mat, y, v.vec as *const ffi::gsl_vector_float) }
        }

        /// This function copies the elements of the vector v into the x-th column of the matrix.
        /// The length of the vector must be the same as the length of the column.
        pub fn set_col(&self, x: u64, v: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_set_col(self.mat, x, v.vec as *const ffi::gsl_vector_float) }
        }

        /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
        pub fn swap_rows(&self, y1: u64, y2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_float_swap_rows(self.mat, y1, y2) }
        }

        /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
        pub fn swap_columns(&self, x1: u64, x2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_float_swap_columns(self.mat, x1, x2) }
        }

        /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
        pub fn swap_row_col(&self, i: u64, j: u64) -> i32 {
            unsafe { ffi::gsl_matrix_float_swap_rowcol(self.mat, i, j) }
        }

        /// This function returns the transpose of the matrix by copying the elements into it.
        /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
        pub fn transpose_memcpy(&self) -> Option<(MatrixFloat, i32)> {
            let dest = unsafe { ffi::gsl_matrix_float_alloc((*self.mat).size1, (*self.mat).size2) };

            if dest.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_float_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix_float) };

                Some((MatrixFloat{mat: dest}, ret))
            }
        }

        /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
        /// The matrix must be square for this operation to be possible.
        pub fn transpose(&self) -> i32 {
            unsafe { ffi::gsl_matrix_float_transpose(self.mat) }
        }

        /// This function adds the elements of the other matrix to the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn add(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_add(self.mat, other.mat as *const ffi::gsl_matrix_float) }
        }

        /// This function subtracts the elements of the other matrix from the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn sub(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_sub(self.mat, other.mat as *const ffi::gsl_matrix_float) }
        }

        /// This function multiplies the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn mul_elements(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix_float) }
        }

        /// This function divides the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn div_elements(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_float_div_elements(self.mat, other.mat as *const ffi::gsl_matrix_float) }
        }

        /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
        pub fn scale(&self, x: f32) -> i32 {
            unsafe { ffi::gsl_matrix_float_scale(self.mat, x) }
        }

        /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
        pub fn add_constant(&self, x: f32) -> i32 {
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
        pub fn equal(&self, other: &MatrixFloat) -> bool {
            match unsafe { ffi::gsl_matrix_float_equal(self.mat as *const ffi::gsl_matrix_float, other.mat as *const ffi::gsl_matrix_float) } {
                1 => true,
                _ => false
            }
        }

        pub fn clone(&self) -> Option<MatrixFloat> {
            unsafe {
                if self.mat.is_null() {
                    None
                } else {
                    match MatrixFloat::new((*self.mat).size1, (*self.mat).size2) {
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

    impl Drop for MatrixFloat {
        fn drop(&mut self) {
            unsafe { ffi::gsl_matrix_float_free(self.mat) };
            self.mat = ::std::ptr::mut_null();
        }
    }

    impl Show for MatrixFloat {
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

    pub struct MatrixComplex {
        mat: *mut ffi::gsl_matrix_complex
    }

    impl MatrixComplex {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_matrix_complex {
            self.mat
        }

        /// Creates a new Matrix with all elements set to zero
        /// 
        /// Example with n1 = 2 and n2 = 3 :
        /// 
        /// XX XX XX
        /// 
        /// XX XX XX
        pub fn new(n1: u64, n2: u64) -> Option<MatrixComplex> {
            let tmp = unsafe { ffi::gsl_matrix_complex_calloc(n1, n2) };

            if tmp.is_null() {
                None
            } else {
                Some(MatrixComplex {
                    mat: tmp
                })
            }
        }

        /// This function returns the (i,j)-th element of the matrix.
        /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
        pub fn get(&self, y: u64, x: u64) -> Complex {
            unsafe { ::std::mem::transmute(ffi::gsl_matrix_complex_get(self.mat as *const ffi::gsl_matrix_complex, y, x)) }
        }

        /// This function sets the value of the (i,j)-th element of the matrix to value.
        /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
        pub fn set(&self, y: u64, x: u64, value: &Complex) {
            unsafe { ffi::gsl_matrix_complex_set(self.mat, y, x, ::std::mem::transmute(*value)) }
        }

        /// This function sets all the elements of the matrix to the value x.
        pub fn set_all(&self, x: &Complex) {
            unsafe { ffi::gsl_matrix_complex_set_all(self.mat, ::std::mem::transmute(*x)) }
        }

        /// This function sets all the elements of the matrix to zero.
        pub fn set_zero(&self) {
            unsafe { ffi::gsl_matrix_complex_set_zero(self.mat) }
        }

        /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
        /// This applies to both square and rectangular matrices.
        pub fn set_identity(&self) {
            unsafe { ffi::gsl_matrix_complex_set_identity(self.mat) }
        }

        /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
        pub fn copy_from(&self, other: &MatrixComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_memcpy(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
        }

        /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
        pub fn copy_to(&self, other: &MatrixComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_memcpy(other.mat, self.mat as *const ffi::gsl_matrix_complex) }
        }

        /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
        pub fn swap(&self, other: &MatrixComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_swap(self.mat, other.mat) }
        }

        /// This function copies the elements of the y-th row of the matrix into the returned vector.
        pub fn get_row(&self, y: u64) -> Option<(VectorComplex, i32)> {
            let tmp = unsafe { ffi::gsl_vector_complex_alloc((*self.mat).size2) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_complex_get_row(tmp, self.mat as *const ffi::gsl_matrix_complex, y) };

                Some((VectorComplex{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the x-th column of the matrix into the returned vector.
        pub fn get_col(&self, x: u64) -> Option<(VectorComplex, i32)> {
            let tmp = unsafe { ffi::gsl_vector_complex_alloc((*self.mat).size1) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_complex_get_col(tmp, self.mat as *const ffi::gsl_matrix_complex, x) };

                Some((VectorComplex{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the vector v into the y-th row of the matrix.
        /// The length of the vector must be the same as the length of the row.
        pub fn set_row(&self, y: u64, v: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_set_row(self.mat, y, v.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function copies the elements of the vector v into the x-th column of the matrix.
        /// The length of the vector must be the same as the length of the column.
        pub fn set_col(&self, x: u64, v: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_set_col(self.mat, x, v.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
        pub fn swap_rows(&self, y1: u64, y2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_complex_swap_rows(self.mat, y1, y2) }
        }

        /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
        pub fn swap_columns(&self, x1: u64, x2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_complex_swap_columns(self.mat, x1, x2) }
        }

        /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
        pub fn swap_row_col(&self, i: u64, j: u64) -> i32 {
            unsafe { ffi::gsl_matrix_complex_swap_rowcol(self.mat, i, j) }
        }

        /// This function returns the transpose of the matrix by copying the elements into it.
        /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
        pub fn transpose_memcpy(&self) -> Option<(MatrixComplex, i32)> {
            let dest = unsafe { ffi::gsl_matrix_complex_alloc((*self.mat).size1, (*self.mat).size2) };

            if dest.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_complex_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix_complex) };

                Some((MatrixComplex{mat: dest}, ret))
            }
        }

        /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
        /// The matrix must be square for this operation to be possible.
        pub fn transpose(&self) -> i32 {
            unsafe { ffi::gsl_matrix_complex_transpose(self.mat) }
        }

        /// This function adds the elements of the other matrix to the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn add(&self, other: &MatrixComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_add(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
        }

        /// This function subtracts the elements of the other matrix from the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn sub(&self, other: &MatrixComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_sub(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
        }

        /// This function multiplies the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn mul_elements(&self, other: &MatrixComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
        }

        /// This function divides the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn div_elements(&self, other: &MatrixComplex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_div_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex) }
        }

        /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
        pub fn scale(&self, x: &Complex) -> i32 {
            unsafe { ffi::gsl_matrix_complex_scale(self.mat, ::std::mem::transmute(*x)) }
        }

        /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
        pub fn add_constant(&self, x: &Complex) -> i32 {
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
        pub fn equal(&self, other: &MatrixComplex) -> bool {
            match unsafe { ffi::gsl_matrix_complex_equal(self.mat as *const ffi::gsl_matrix_complex, other.mat as *const ffi::gsl_matrix_complex) } {
                1 => true,
                _ => false
            }
        }

        pub fn clone(&self) -> Option<MatrixComplex> {
            unsafe {
                if self.mat.is_null() {
                    None
                } else {
                    match MatrixComplex::new((*self.mat).size1, (*self.mat).size2) {
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

    impl Drop for MatrixComplex {
        fn drop(&mut self) {
            unsafe { ffi::gsl_matrix_complex_free(self.mat) };
            self.mat = ::std::ptr::mut_null();
        }
    }

    impl Show for MatrixComplex {
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






    pub struct MatrixComplexFloat {
        mat: *mut ffi::gsl_matrix_complex_float
    }

    impl MatrixComplexFloat {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_matrix_complex_float {
            self.mat
        }

        /// Creates a new Matrix with all elements set to zero
        /// 
        /// Example with n1 = 2 and n2 = 3 :
        /// 
        /// XX XX XX
        /// 
        /// XX XX XX
        pub fn new(n1: u64, n2: u64) -> Option<MatrixComplexFloat> {
            let tmp = unsafe { ffi::gsl_matrix_complex_float_calloc(n1, n2) };

            if tmp.is_null() {
                None
            } else {
                Some(MatrixComplexFloat {
                    mat: tmp
                })
            }
        }

        /// This function returns the (i,j)-th element of the matrix.
        /// If y or x lie outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked and 0 is returned.
        pub fn get(&self, y: u64, x: u64) -> ComplexFloat {
            unsafe { ::std::mem::transmute(ffi::gsl_matrix_complex_float_get(self.mat as *const ffi::gsl_matrix_complex_float, y, x)) }
        }

        /// This function sets the value of the (i,j)-th element of the matrix to value.
        /// If y or x lies outside the allowed range of 0 to n1-1 and 0 to n2-1 then the error handler is invoked.
        pub fn set(&self, y: u64, x: u64, value: &ComplexFloat) {
            unsafe { ffi::gsl_matrix_complex_float_set(self.mat, y, x, ::std::mem::transmute(*value)) }
        }

        /// This function sets all the elements of the matrix to the value x.
        pub fn set_all(&self, x: &ComplexFloat) {
            unsafe { ffi::gsl_matrix_complex_float_set_all(self.mat, ::std::mem::transmute(*x)) }
        }

        /// This function sets all the elements of the matrix to zero.
        pub fn set_zero(&self) {
            unsafe { ffi::gsl_matrix_complex_float_set_zero(self.mat) }
        }

        /// This function sets the elements of the matrix to the corresponding elements of the identity matrix, m(i,j) = \delta(i,j), i.e. a unit diagonal with all off-diagonal elements zero.
        /// This applies to both square and rectangular matrices.
        pub fn set_identity(&self) {
            unsafe { ffi::gsl_matrix_complex_float_set_identity(self.mat) }
        }

        /// This function copies the elements of the other matrix into the self matrix. The two matrices must have the same size.
        pub fn copy_from(&self, other: &MatrixComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_memcpy(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
        }

        /// This function copies the elements of the self matrix into the other matrix. The two matrices must have the same size.
        pub fn copy_to(&self, other: &MatrixComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_memcpy(other.mat, self.mat as *const ffi::gsl_matrix_complex_float) }
        }

        /// This function exchanges the elements of the matrices self and other by copying. The two matrices must have the same size.
        pub fn swap(&self, other: &MatrixComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_swap(self.mat, other.mat) }
        }

        /// This function copies the elements of the y-th row of the matrix into the returned vector.
        pub fn get_row(&self, y: u64) -> Option<(VectorComplexFloat, i32)> {
            let tmp = unsafe { ffi::gsl_vector_complex_float_alloc((*self.mat).size2) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_complex_float_get_row(tmp, self.mat as *const ffi::gsl_matrix_complex_float, y) };

                Some((VectorComplexFloat{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the x-th column of the matrix into the returned vector.
        pub fn get_col(&self, x: u64) -> Option<(VectorComplexFloat, i32)> {
            let tmp = unsafe { ffi::gsl_vector_complex_float_alloc((*self.mat).size1) };

            if tmp.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_complex_float_get_col(tmp, self.mat as *const ffi::gsl_matrix_complex_float, x) };

                Some((VectorComplexFloat{vec: tmp}, ret))
            }
        }

        /// This function copies the elements of the vector v into the y-th row of the matrix.
        /// The length of the vector must be the same as the length of the row.
        pub fn set_row(&self, y: u64, v: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_set_row(self.mat, y, v.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function copies the elements of the vector v into the x-th column of the matrix.
        /// The length of the vector must be the same as the length of the column.
        pub fn set_col(&self, x: u64, v: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_set_col(self.mat, x, v.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function exchanges the y1-th and y2-th rows of the matrix in-place.
        pub fn swap_rows(&self, y1: u64, y2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_swap_rows(self.mat, y1, y2) }
        }

        /// This function exchanges the x1-th and x2-th columns of the matrix in-place.
        pub fn swap_columns(&self, x1: u64, x2: u64) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_swap_columns(self.mat, x1, x2) }
        }

        /// This function exchanges the i-th row and j-th column of the matrix in-place. The matrix must be square for this operation to be possible.
        pub fn swap_row_col(&self, i: u64, j: u64) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_swap_rowcol(self.mat, i, j) }
        }

        /// This function returns the transpose of the matrix by copying the elements into it.
        /// This function works for all matrices provided that the dimensions of the matrix dest match the transposed dimensions of the matrix.
        pub fn transpose_memcpy(&self) -> Option<(MatrixComplexFloat, i32)> {
            let dest = unsafe { ffi::gsl_matrix_complex_float_alloc((*self.mat).size1, (*self.mat).size2) };

            if dest.is_null() {
                None
            } else {
                let ret = unsafe { ffi::gsl_matrix_complex_float_transpose_memcpy(dest, self.mat as *const ffi::gsl_matrix_complex_float) };

                Some((MatrixComplexFloat{mat: dest}, ret))
            }
        }

        /// This function replaces the matrix m by its transpose by copying the elements of the matrix in-place.
        /// The matrix must be square for this operation to be possible.
        pub fn transpose(&self) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_transpose(self.mat) }
        }

        /// This function adds the elements of the other matrix to the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) + other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn add(&self, other: &MatrixComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_add(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
        }

        /// This function subtracts the elements of the other matrix from the elements of the self matrix.
        /// The result self(i,j) <- self(i,j) - other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn sub(&self, other: &MatrixComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_sub(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
        }

        /// This function multiplies the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) * other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn mul_elements(&self, other: &MatrixComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_mul_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
        }

        /// This function divides the elements of the self matrix by the elements of the other matrix.
        /// The result self(i,j) <- self(i,j) / other(i,j) is stored in self and other remains unchanged. The two matrices must have the same dimensions.
        pub fn div_elements(&self, other: &MatrixFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_div_elements(self.mat, other.mat as *const ffi::gsl_matrix_complex_float) }
        }

        /// This function multiplies the elements of the self matrix by the constant factor x. The result self(i,j) <- x self(i,j) is stored in self.
        pub fn scale(&self, x: &ComplexFloat) -> i32 {
            unsafe { ffi::gsl_matrix_complex_float_scale(self.mat, ::std::mem::transmute(*x)) }
        }

        /// This function adds the constant value x to the elements of the self matrix. The result self(i,j) <- self(i,j) + x is stored in self.
        pub fn add_constant(&self, x: &ComplexFloat) -> i32 {
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
        pub fn equal(&self, other: &MatrixComplexFloat) -> bool {
            match unsafe { ffi::gsl_matrix_complex_float_equal(self.mat as *const ffi::gsl_matrix_complex_float,
                other.mat as *const ffi::gsl_matrix_complex_float) } {
                1 => true,
                _ => false
            }
        }

        pub fn clone(&self) -> Option<MatrixComplexFloat> {
            unsafe {
                if self.mat.is_null() {
                    None
                } else {
                    match MatrixComplexFloat::new((*self.mat).size1, (*self.mat).size2) {
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

    impl Drop for MatrixComplexFloat {
        fn drop(&mut self) {
            unsafe { ffi::gsl_matrix_complex_float_free(self.mat) };
            self.mat = ::std::ptr::mut_null();
        }
    }

    impl Show for MatrixComplexFloat {
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
}