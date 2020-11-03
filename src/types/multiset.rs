//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Multisets

This chapter describes functions for creating and manipulating multisets. A multiset c is represented by an array of k integers in the
range 0 to n-1, where each value c_i may occur more than once. The multiset c corresponds to indices of k elements chosen from an n element
vector with replacement. In mathematical terms, n is the cardinality of the multiset while k is the maximum multiplicity of any value.
Multisets are useful, for example, when iterating over the indices of a k-th order symmetric tensor in n-space.
!*/

use crate::Value;
use c_vec::CSlice;
use ffi;
use std::io::Result as IoResult;
use std::io::Write;

pub struct MultiSet {
    c: *mut sys::gsl_multiset,
    data: CSlice<usize>,
}

impl MultiSet {
    /// This function allocates memory for a new multiset with parameters n, k. The multiset is not initialized and its elements are
    /// undefined. Use the function gsl_multiset_calloc if you want to create a multiset which is initialized to the lexicographically
    /// first multiset element. A null pointer is returned if insufficient memory is available to create the multiset.
    pub fn new(n: usize, k: usize) -> Option<MultiSet> {
        let tmp = unsafe { sys::gsl_multiset_alloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                if (*tmp).data.is_null() {
                    Some(MultiSet {
                        c: tmp,
                        // dirty trick to avoid a failure
                        data: CSlice::new(tmp as *mut _, 0),
                    })
                } else {
                    Some(MultiSet {
                        c: tmp,
                        data: CSlice::new((*tmp).data, (*tmp).k as _),
                    })
                }
            }
        }
    }

    /// This function allocates memory for a new multiset with parameters n, k and initializes it to the lexicographically first multiset
    /// element. A null pointer is returned if insufficient memory is available to create the multiset.
    pub fn new_init(n: usize, k: usize) -> Option<MultiSet> {
        let tmp = unsafe { sys::gsl_multiset_calloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                if (*tmp).data.is_null() {
                    Some(MultiSet {
                        c: tmp,
                        // dirty trick to avoid a failure
                        data: CSlice::new(tmp as *mut _, 0),
                    })
                } else {
                    Some(MultiSet {
                        c: tmp,
                        data: CSlice::new((*tmp).data, (*tmp).k as _),
                    })
                }
            }
        }
    }

    /// This function initializes the multiset c to the lexicographically first multiset element, i.e. 0 repeated k times.
    pub fn init_first(&mut self) {
        unsafe { sys::gsl_multiset_init_first(self.c) }
    }

    /// This function initializes the multiset c to the lexicographically last multiset element, i.e. n-1 repeated k times.
    pub fn init_last(&mut self) {
        unsafe { sys::gsl_multiset_init_last(self.c) }
    }

    /// This function copies the elements of the multiset self into the multiset dest. The two multisets must have the same size.
    pub fn copy(&self, dest: &mut MultiSet) -> Value {
        Value::from(unsafe { sys::gsl_multiset_memcpy(dest.c, self.c) })
    }

    /// This function returns the value of the i-th element of the multiset c. If i lies outside the allowed range of 0 to k-1 then the
    /// error handler is invoked and 0 is returned.
    pub fn get(&self, i: usize) -> usize {
        unsafe { sys::gsl_multiset_get(self.c, i) }
    }

    /// This function returns the range (n) of the multiset self.
    pub fn n(&self) -> usize {
        unsafe { sys::gsl_multiset_n(self.c) }
    }

    /// This function returns the number of elements (k) in the multiset self.
    pub fn k(&self) -> usize {
        unsafe { sys::gsl_multiset_k(self.c) }
    }

    /// This function returns a pointer to the array of elements in the multiset self.
    pub fn data<'r>(&'r mut self) -> &'r mut [usize] {
        self.data.as_mut()
    }

    /// This function checks that the multiset self is valid. The k elements should lie in the range 0 to n-1, with each value occurring in
    /// nondecreasing order.
    pub fn valid(&self) -> Value {
        Value::from(unsafe { sys::gsl_multiset_valid(self.c) })
    }

    /// This function advances the multiset self to the next multiset element in lexicographic order and returns ::Value::Success. If no
    /// further multisets elements are available it returns Value::Failure and leaves self unmodified. Starting with the first multiset and
    /// repeatedly applying this function will iterate through all possible multisets of a given order.
    pub fn next(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_multiset_next(self.c) })
    }

    /// This function steps backwards from the multiset self to the previous multiset element in lexicographic order, returning ::Value::Success.
    /// If no previous multiset is available it returns Value::Failure and leaves self unmodified.
    pub fn prev(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_multiset_prev(self.c) })
    }

    pub fn print<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        for value in self.data.as_ref().iter() {
            write!(writer, " {}", *value)?;
        }
        Ok(())
    }
}

impl Drop for MultiSet {
    fn drop(&mut self) {
        unsafe { sys::gsl_multiset_free(self.c) };
        self.c = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<sys::gsl_multiset> for MultiSet {
    fn wrap(c: *mut sys::gsl_multiset) -> Self {
        unsafe {
            if (*c).data.is_null() {
                Self {
                    c,
                    // dirty trick to avoid a failure
                    data: CSlice::new(c as *mut _, 0),
                }
            } else {
                Self {
                    c,
                    data: CSlice::new((*c).data, (*c).k as _),
                }
            }
        }
    }

    fn soft_wrap(c: *mut sys::gsl_multiset) -> Self {
        Self::wrap(c)
    }

    fn unwrap_shared(&self) -> *const sys::gsl_multiset {
        self.c as *const _
    }

    fn unwrap_unique(&mut self) -> *mut sys::gsl_multiset {
        self.c
    }
}
