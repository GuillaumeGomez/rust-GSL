//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub type gsl_mode_t = u32;
pub struct CblasIndex(pub u32);

pub mod Gsl {
    use std::default::Default;
    use ffi;

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
        pub fn new(size: u32) -> Option<VectorFloat> {
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
            let tmp = unsafe { ffi::gsl_vector_float_alloc(slice.len() as u32) };

            if tmp.is_null() {
                None
            } else {
                let v = VectorFloat {
                    vec: tmp
                };
                let mut pos = 0u32;

                for tmp in slice.iter() {
                    v.set(pos, *tmp);
                    pos += 1;
                }
                Some(v)
            }
        }

        pub fn len(&self) -> u32 {
            if self.vec.is_null() {
                0u32
            } else {
                unsafe { (*self.vec).size }
            }
        }

        /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
        pub fn get(&self, i: u32) -> f32 {
            unsafe { ffi::gsl_vector_float_get(self.vec, i) }
        }

        /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
        pub fn set(&self, i: u32, x: f32) {
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
        pub fn set_basis(&self, i: u32) {
            unsafe { ffi::gsl_vector_float_set_basis(self.vec, i) }
        }

        /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
        pub fn copy_from(&self, other: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_vector_float_memcpy(self.vec, other.vec as *const ffi::gsl_vector_float) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &VectorFloat) -> i32 {
            unsafe { ffi::gsl_vector_float_swap(other.vec, self.vec) }
        }

        /// This function exchanges the i-th and j-th elements of the vector v in-place.
        pub fn swap_elements(&self, i: u32, j: u32) -> i32 {
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
            unsafe { ffi::gsl_vector_float_max(self.vec) }
        }

        /// This function returns the minimum value in the self vector.
        pub fn min(&self) -> f32 {
            unsafe { ffi::gsl_vector_float_min(self.vec) }
        }

        /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
        pub fn minmax(&self, min_out: &mut f32, max_out: &mut f32) {
            unsafe { ffi::gsl_vector_float_minmax(self.vec, min_out, max_out) }
        }

        /// This function returns the index of the maximum value in the self vector.
        /// When there are several equal maximum elements then the lowest index is returned.
        pub fn max_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_float_max_index(self.vec) }
        }

        /// This function returns the index of the minimum value in the self vector.
        /// When there are several equal minimum elements then the lowest index is returned.
        pub fn min_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_float_min_index(self.vec) }
        }

        /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
        /// When there are several equal minimum or maximum elements then the lowest indices are returned.
        pub fn minmax_index(&self) -> (u32, u32) {
            let mut imin = 0u32;
            let mut imax = 0u32;

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

    pub struct Vector {
        vec: *mut ffi::gsl_vector
    }

    impl Vector {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_vector {
            self.vec
        }

        /// create a new VectorFloat with all elements set to zero
        pub fn new(size: u32) -> Option<Vector> {
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
            let tmp = unsafe { ffi::gsl_vector_alloc(slice.len() as u32) };

            if tmp.is_null() {
                None
            } else {
                let v = Vector {
                    vec: tmp
                };
                let mut pos = 0u32;

                for tmp in slice.iter() {
                    v.set(pos, *tmp);
                    pos += 1;
                }
                Some(v)
            }
        }

        pub fn len(&self) -> u32 {
            if self.vec.is_null() {
                0u32
            } else {
                unsafe { (*self.vec).size }
            }
        }

        /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
        pub fn get(&self, i: u32) -> f64 {
            unsafe { ffi::gsl_vector_get(self.vec, i) }
        }

        /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
        pub fn set(&self, i: u32, x: f64) {
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
        pub fn set_basis(&self, i: u32) {
            unsafe { ffi::gsl_vector_set_basis(self.vec, i) }
        }

        /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
        pub fn copy_from(&self, other: &Vector) -> i32 {
            unsafe { ffi::gsl_vector_memcpy(self.vec, other.vec as *const ffi::gsl_vector) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &Vector) -> i32 {
            unsafe { ffi::gsl_vector_swap(other.vec, self.vec) }
        }

        /// This function exchanges the i-th and j-th elements of the vector v in-place.
        pub fn swap_elements(&self, i: u32, j: u32) -> i32 {
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
            unsafe { ffi::gsl_vector_max(self.vec) }
        }

        /// This function returns the minimum value in the self vector.
        pub fn min(&self) -> f64 {
            unsafe { ffi::gsl_vector_min(self.vec) }
        }

        /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
        pub fn minmax(&self, min_out: &mut f64, max_out: &mut f64) {
            unsafe { ffi::gsl_vector_minmax(self.vec, min_out, max_out) }
        }

        /// This function returns the index of the maximum value in the self vector.
        /// When there are several equal maximum elements then the lowest index is returned.
        pub fn max_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_max_index(self.vec) }
        }

        /// This function returns the index of the minimum value in the self vector.
        /// When there are several equal minimum elements then the lowest index is returned.
        pub fn min_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_min_index(self.vec) }
        }

        /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
        /// When there are several equal minimum or maximum elements then the lowest indices are returned.
        pub fn minmax_index(&self) -> (u32, u32) {
            let mut imin = 0u32;
            let mut imax = 0u32;

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


    pub struct VectorComplex {
        vec: *mut ffi::gsl_vector_complex
    }

    impl VectorComplex {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_vector_complex {
            self.vec
        }

        /// create a new VectorFloat with all elements set to zero
        pub fn new(size: u32) -> Option<VectorComplex> {
            let tmp = unsafe { ffi::gsl_vector_complex_calloc(size) };

            if tmp.is_null() {
                None
            } else {
                Some(VectorComplex {
                    vec: tmp
                })
            }
        }

        pub fn from_slice(slice: &[f64]) -> Option<VectorComplex> {
            let tmp = unsafe { ffi::gsl_vector_complex_alloc(slice.len() as u32) };

            if tmp.is_null() {
                None
            } else {
                let v = VectorComplex {
                    vec: tmp
                };
                let mut pos = 0u32;

                for tmp in slice.iter() {
                    v.set(pos, *tmp);
                    pos += 1;
                }
                Some(v)
            }
        }

        pub fn len(&self) -> u32 {
            if self.vec.is_null() {
                0u32
            } else {
                unsafe { (*self.vec).size }
            }
        }

        /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
        pub fn get(&self, i: u32) -> f64 {
            unsafe { ffi::gsl_vector_complex_get(self.vec, i) }
        }

        /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
        pub fn set(&self, i: u32, x: f64) {
            unsafe { ffi::gsl_vector_complex_set(self.vec, i, x) }
        }

        /// This function sets all the elements of the vector v to the value x.
        pub fn set_all(&self, x: f64) {
            unsafe { ffi::gsl_vector_complex_set_all(self.vec, x) }
        }

        /// This function sets all the elements of the vector v to zero.
        pub fn set_zero(&self) {
            unsafe { ffi::gsl_vector_complex_set_zero(self.vec) }
        }

        /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
        pub fn set_basis(&self, i: u32) {
            unsafe { ffi::gsl_vector_complex_set_basis(self.vec, i) }
        }

        /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
        pub fn copy_from(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_memcpy(self.vec, other.vec as *const ffi::gsl_vector_complex) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &VectorComplex) -> i32 {
            unsafe { ffi::gsl_vector_complex_swap(other.vec, self.vec) }
        }

        /// This function exchanges the i-th and j-th elements of the vector v in-place.
        pub fn swap_elements(&self, i: u32, j: u32) -> i32 {
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
        pub fn scale(&self, x: f64) -> i32 {
            unsafe { ffi::gsl_vector_complex_scale(self.vec, x) }
        }

        /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
        pub fn add_constant(&self, x: f64) -> i32 {
            unsafe { ffi::gsl_vector_complex_add_constant(self.vec, x) }
        }

        /// This function returns the maximum value in the self vector.
        pub fn max(&self) -> f64 {
            unsafe { ffi::gsl_vector_complex_max(self.vec) }
        }

        /// This function returns the minimum value in the self vector.
        pub fn min(&self) -> f64 {
            unsafe { ffi::gsl_vector_complex_min(self.vec) }
        }

        /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
        pub fn minmax(&self, min_out: &mut f64, max_out: &mut f64) {
            unsafe { ffi::gsl_vector_complex_minmax(self.vec, min_out, max_out) }
        }

        /// This function returns the index of the maximum value in the self vector.
        /// When there are several equal maximum elements then the lowest index is returned.
        pub fn max_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_complex_max_index(self.vec) }
        }

        /// This function returns the index of the minimum value in the self vector.
        /// When there are several equal minimum elements then the lowest index is returned.
        pub fn min_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_complex_min_index(self.vec) }
        }

        /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
        /// When there are several equal minimum or maximum elements then the lowest indices are returned.
        pub fn minmax_index(&self) -> (u32, u32) {
            let mut imin = 0u32;
            let mut imax = 0u32;

            unsafe { ffi::gsl_vector_complex_minmax_index(self.vec, &mut imin, &mut imax) };
            (imin, imax)
        }

        /// This function returns true if all the elements of the self vector are equal to 0.
        pub fn is_null(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_isnull(self.vec) } {
                1 => true,
                _ => false
            }
        }

        /// This function returns true if all the elements of the self vector are stricly positive.
        pub fn is_pos(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_ispos(self.vec) } {
                1 => true,
                _ => false
            }
        }

        /// This function returns true if all the elements of the self vector are stricly negative.
        pub fn is_neg(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_isneg(self.vec) } {
                1 => true,
                _ => false
            }
        }

        /// This function returns true if all the elements of the self vector are stricly non-negative.
        pub fn is_non_neg(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_isnonneg(self.vec) } {
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

    pub struct VectorComplexFloat {
        vec: *mut ffi::gsl_vector_complex_float
    }

    impl VectorComplexFloat {
        #[doc(hidden)]
        #[allow(visible_private_types)]
        pub fn get_ffi(&self) -> *mut ffi::gsl_vector_complex_float {
            self.vec
        }

        /// create a new VectorFloat with all elements set to zero
        pub fn new(size: u32) -> Option<VectorComplexFloat> {
            let tmp = unsafe { ffi::gsl_vector_complex_float_calloc(size) };

            if tmp.is_null() {
                None
            } else {
                Some(VectorComplexFloat {
                    vec: tmp
                })
            }
        }

        pub fn from_slice(slice: &[f32]) -> Option<VectorComplexFloat> {
            let tmp = unsafe { ffi::gsl_vector_complex_float_alloc(slice.len() as u32) };

            if tmp.is_null() {
                None
            } else {
                let v = VectorComplexFloat {
                    vec: tmp
                };
                let mut pos = 0u32;

                for tmp in slice.iter() {
                    v.set(pos, *tmp);
                    pos += 1;
                }
                Some(v)
            }
        }

        pub fn len(&self) -> u32 {
            if self.vec.is_null() {
                0u32
            } else {
                unsafe { (*self.vec).size }
            }
        }

        /// This function returns the i-th element of a vector v. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.
        pub fn get(&self, i: u32) -> f32 {
            unsafe { ffi::gsl_vector_complex_float_get(self.vec, i) }
        }

        /// This function sets the value of the i-th element of a vector v to x. If i lies outside the allowed range of 0 to n-1 then the error handler is invoked.
        pub fn set(&self, i: u32, x: f32) {
            unsafe { ffi::gsl_vector_complex_float_set(self.vec, i, x) }
        }

        /// This function sets all the elements of the vector v to the value x.
        pub fn set_all(&self, x: f32) {
            unsafe { ffi::gsl_vector_complex_float_set_all(self.vec, x) }
        }

        /// This function sets all the elements of the vector v to zero.
        pub fn set_zero(&self) {
            unsafe { ffi::gsl_vector_complex_float_set_zero(self.vec) }
        }

        /// This function makes a basis vector by setting all the elements of the vector v to zero except for the i-th element which is set to one.
        pub fn set_basis(&self, i: u32) {
            unsafe { ffi::gsl_vector_complex_float_set_basis(self.vec, i) }
        }

        /// This function copies the elements of the other vector into the self vector. The two vectors must have the same length.
        pub fn copy_from(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_memcpy(self.vec, other.vec as *const ffi::gsl_vector_complex_float) }
        }

        /// This function exchanges the elements of the vectors by copying. The two vectors must have the same length.
        pub fn swap(&self, other: &VectorComplexFloat) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_swap(other.vec, self.vec) }
        }

        /// This function exchanges the i-th and j-th elements of the vector v in-place.
        pub fn swap_elements(&self, i: u32, j: u32) -> i32 {
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
        pub fn scale(&self, x: f32) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_scale(self.vec, x) }
        }

        /// This function adds the constant value x to the elements of the self vector. The result a_i <- a_i + x is stored in self.
        pub fn add_constant(&self, x: f32) -> i32 {
            unsafe { ffi::gsl_vector_complex_float_add_constant(self.vec, x) }
        }

        /// This function returns the maximum value in the self vector.
        pub fn max(&self) -> f32 {
            unsafe { ffi::gsl_vector_complex_float_max(self.vec) }
        }

        /// This function returns the minimum value in the self vector.
        pub fn min(&self) -> f32 {
            unsafe { ffi::gsl_vector_complex_float_min(self.vec) }
        }

        /// This function returns the minimum and maximum values in the self vector, storing them in min_out and max_out.
        pub fn minmax(&self, min_out: &mut f32, max_out: &mut f32) {
            unsafe { ffi::gsl_vector_complex_float_minmax(self.vec, min_out, max_out) }
        }

        /// This function returns the index of the maximum value in the self vector.
        /// When there are several equal maximum elements then the lowest index is returned.
        pub fn max_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_complex_float_max_index(self.vec) }
        }

        /// This function returns the index of the minimum value in the self vector.
        /// When there are several equal minimum elements then the lowest index is returned.
        pub fn min_index(&self) -> u32 {
            unsafe { ffi::gsl_vector_complex_float_min_index(self.vec) }
        }

        /// This function returns the indices of the minimum and maximum values in the self vector, storing them in imin and imax.
        /// When there are several equal minimum or maximum elements then the lowest indices are returned.
        pub fn minmax_index(&self) -> (u32, u32) {
            let mut imin = 0u32;
            let mut imax = 0u32;

            unsafe { ffi::gsl_vector_complex_float_minmax_index(self.vec, &mut imin, &mut imax) };
            (imin, imax)
        }

        /// This function returns true if all the elements of the self vector are equal to 0.
        pub fn is_null(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_float_isnull(self.vec) } {
                1 => true,
                _ => false
            }
        }

        /// This function returns true if all the elements of the self vector are stricly positive.
        pub fn is_pos(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_float_ispos(self.vec) } {
                1 => true,
                _ => false
            }
        }

        /// This function returns true if all the elements of the self vector are stricly negative.
        pub fn is_neg(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_float_isneg(self.vec) } {
                1 => true,
                _ => false
            }
        }

        /// This function returns true if all the elements of the self vector are stricly non-negative.
        pub fn is_non_neg(&self) -> bool {
            match unsafe { ffi::gsl_vector_complex_float_isnonneg(self.vec) } {
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
}