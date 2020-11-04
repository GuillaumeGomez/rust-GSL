//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::paste::paste;
use crate::Value;
use ffi::FFI;
use std::fmt;
use std::fmt::{Debug, Formatter};
use types::{ComplexF32, ComplexF64};

macro_rules! gsl_vec_complex {
    ($rust_name:ident, $name:ident, $complex:ident) => {
        paste! {

        ffi_wrapper!(
            $rust_name,
            *mut sys::$name,
            [<$name _free>]
        );

        impl $rust_name {
            doc! {
                concat!("Create a new ", stringify!($rust_name), "with all elements set to zero"),
                pub fn new(size: usize) -> Option<Self> {
                    let tmp = unsafe { sys::[<$name _calloc>](size) };

                    if tmp.is_null() {
                        None
                    } else {
                        Some(Self::wrap(tmp))
                    }
                }
            }

            pub fn from_slice(slice: &[$complex]) -> Option<Self> {
                let tmp = unsafe { sys::[<$name _alloc>](slice.len() as _) };

                if tmp.is_null() {
                    None
                } else {
                    let mut v = Self::wrap(tmp);

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

            /// This function returns the i-th element of a vector v. If i lies outside the allowed range of
            /// 0 to n-1 then the error handler is invoked and 0 is returned.
            pub fn get(&self, i: usize) -> $complex {
                unsafe { ::std::mem::transmute(sys::[<$name _get>](self.unwrap_shared(), i)) }
            }

            /// This function sets the value of the i-th element of a vector v to x. If i lies outside the
            /// allowed range of 0 to n-1 then the error handler is invoked.
            pub fn set(&mut self, i: usize, x: &$complex) -> &Self {
                unsafe {
                    sys::[<$name _set>](self.unwrap_unique(), i, ::std::mem::transmute(*x))
                };
                self
            }

            /// This function sets all the elements of the vector v to the value x.
            pub fn set_all(&mut self, x: &$complex) -> &Self {
                unsafe {
                    sys::[<$name _set_all>](self.unwrap_unique(), ::std::mem::transmute(*x))
                };
                self
            }

            /// This function sets all the elements of the vector v to zero.
            pub fn set_zero(&mut self) -> &$rust_name {
                unsafe { sys::[<$name _set_zero>](self.unwrap_unique()) };
                self
            }

            /// This function makes a basis vector by setting all the elements of the vector v to zero
            /// except for the i-th element which is set to one.
            pub fn set_basis(&mut self, i: usize) -> &Self {
                unsafe { sys::[<$name _set_basis>](self.unwrap_unique(), i) };
                self
            }

            /// This function copies the elements of the other vector into the self vector. The two vectors
            /// must have the same length.
            pub fn copy_from(&mut self, other: &$rust_name) -> Value {
                Value::from(unsafe {
                    sys::[<$name _memcpy>](self.unwrap_unique(), other.unwrap_shared())
                })
            }

            /// This function copies the elements of the self vector into the other vector. The two vectors
            /// must have the same length.
            pub fn copy_to(&self, other: &mut $rust_name) -> Value {
                Value::from(unsafe {
                    sys::[<$name _memcpy>](other.unwrap_unique(), self.unwrap_shared())
                })
            }

            /// This function exchanges the elements of the vectors by copying. The two vectors must have
            /// the same length.
            pub fn swap(&mut self, other: &mut $rust_name) -> Value {
                Value::from(unsafe {
                    sys::[<$name _swap>](other.unwrap_unique(), self.unwrap_unique())
                })
            }

            /// This function exchanges the i-th and j-th elements of the vector v in-place.
            pub fn swap_elements(&mut self, i: usize, j: usize) -> Value {
                Value::from(unsafe {
                    sys::[<$name _swap_elements>](self.unwrap_unique(), i, j)
                })
            }

            /// This function reverses the order of the elements of the vector v.
            pub fn reverse(&mut self) -> Value {
                Value::from(unsafe { sys::[<$name _reverse>](self.unwrap_unique()) })
            }

            /// This function adds the elements of the other vector to the elements of the `self` vector.
            /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors
            /// must have the same length.
            pub fn add(&mut self, other: &$rust_name) -> Value {
                Value::from(unsafe {
                    sys::[<$name _add>](self.unwrap_unique(), other.unwrap_shared())
                })
            }

            /// This function subtracts the elements of the self vector from the elements of the other
            /// vector. The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two
            /// vectors must have the same length.
            pub fn sub(&mut self, other: &$rust_name) -> Value {
                Value::from(unsafe {
                    sys::[<$name _sub>](self.unwrap_unique(), other.unwrap_shared())
                })
            }

            /// This function multiplies the elements of the self vector a by the elements of the other
            /// vector. The result a_i <- a_i * b_i is stored in self and other remains unchanged. The two
            /// vectors must have the same length.
            pub fn mul(&mut self, other: &$rust_name) -> Value {
                Value::from(unsafe {
                    sys::[<$name _mul>](self.unwrap_unique(), other.unwrap_shared())
                })
            }

            /// This function divides the elements of the self vector by the elements of the other vector.
            /// The result a_i <- a_i / b_i is stored in self and other remains unchanged. The two vectors
            /// must have the same length.
            pub fn div(&mut self, other: &$rust_name) -> Value {
                Value::from(unsafe {
                    sys::[<$name _div>](self.unwrap_unique(), other.unwrap_shared())
                })
            }

            /// This function multiplies the elements of the self vector by the constant factor x. The
            /// result a_i <- a_i is stored in self.
            pub fn scale(&mut self, x: &$complex) -> Value {
                Value::from(unsafe {
                    sys::[<$name _scale>](self.unwrap_unique(), ::std::mem::transmute(*x))
                })
            }

            /// This function adds the constant value x to the elements of the self vector. The result
            /// a_i <- a_i + x is stored in self.
            pub fn add_constant(&mut self, x: &$complex) -> Value {
                Value::from(unsafe {
                    sys::[<$name _add_constant>](
                        self.unwrap_unique(),
                        ::std::mem::transmute(*x),
                    )
                })
            }

            /// This function returns true if all the elements of the self vector are equal to 0.
            pub fn is_null(&self) -> bool {
                unsafe { sys::[<$name _isnull>](self.unwrap_shared()) == 1 }
            }

            /// This function returns true if all the elements of the self vector are stricly positive.
            pub fn is_pos(&self) -> bool {
                unsafe { sys::[<$name _ispos>](self.unwrap_shared()) == 1 }
            }

            /// This function returns true if all the elements of the self vector are stricly negative.
            pub fn is_neg(&self) -> bool {
                unsafe { sys::[<$name _isneg>](self.unwrap_shared()) == 1 }
            }

            /// This function returns true if all the elements of the self vector are stricly non-negative.
            pub fn is_non_neg(&self) -> bool {
                unsafe { sys::[<$name _isnonneg>](self.unwrap_shared()) == 1 }
            }

            pub fn equal(&self, other: &$rust_name) -> bool {
                unsafe {
                    sys::[<$name _equal>](self.unwrap_shared(), other.unwrap_shared()) == 1
                }
            }

            pub fn clone(&self) -> Option<Self> {
                if self.unwrap_shared().is_null() {
                    None
                } else {
                    unsafe {
                        match Self::new(self.len()) {
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

        impl Debug for $rust_name {
            fn fmt(&self, f: &mut Formatter) -> fmt::Result {
                let ptr = self.unwrap_shared();
                if ptr.is_null() {
                    write!(f, "<null>")
                } else {
                    unsafe {
                        write!(f, "[").unwrap();
                        let size = self.len();
                        for x in 0..size {
                            if x < size - 1 {
                                write!(f, "{:?}, ", self.get(x)).unwrap();
                            } else {
                                write!(f, "{:?}", self.get(x)).unwrap();
                            }
                        }
                    }
                    write!(f, "]")
                }
            }
        }

        } // end of paste! block
    }; // end of macro block
}

gsl_vec_complex!(VectorComplexF32, gsl_vector_complex_float, ComplexF32);
gsl_vec_complex!(VectorComplexF64, gsl_vector_complex, ComplexF64);
