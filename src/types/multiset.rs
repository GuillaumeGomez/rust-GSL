//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Multisets

This chapter describes functions for creating and manipulating multisets. A multiset c is represented by an array of k integers in the
range 0 to n-1, where each value c_i may occur more than once. The multiset c corresponds to indices of k elements chosen from an n element
vector with replacement. In mathematical terms, n is the cardinality of the multiset while k is the maximum multiplicity of any value.
Multisets are useful, for example, when iterating over the indices of a k-th order symmetric tensor in n-space.
!*/

use crate::Value;
use ffi::FFI;
use std::io;

ffi_wrapper!(MultiSet, *mut sys::gsl_multiset, gsl_multiset_free);

impl MultiSet {
    /// This function allocates memory for a new multiset with parameters n, k. The multiset is not
    /// initialized and its elements are undefined. Use the function [`Self::new_with_init`] if you
    /// want to create a multiset which is initialized to the lexicographically first multiset
    /// element. A null pointer is returned if insufficient memory is available to create the
    /// multiset.
    #[doc(alias = "gsl_multiset_alloc")]
    pub fn new(n: usize, k: usize) -> Option<Self> {
        let tmp = unsafe { sys::gsl_multiset_alloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function allocates memory for a new multiset with parameters n, k and initializes it to
    /// the lexicographically first multiset element. A null pointer is returned if insufficient
    /// memory is available to create the multiset.
    #[doc(alias = "gsl_multiset_calloc")]
    pub fn new_with_init(n: usize, k: usize) -> Option<Self> {
        let tmp = unsafe { sys::gsl_multiset_calloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function initializes the multiset c to the lexicographically first multiset element,
    /// i.e. 0 repeated k times.
    #[doc(alias = "gsl_multiset_init_first")]
    pub fn init_first(&mut self) {
        unsafe { sys::gsl_multiset_init_first(self.unwrap_unique()) }
    }

    /// This function initializes the multiset c to the lexicographically last multiset element,
    /// i.e. n-1 repeated k times.
    #[doc(alias = "gsl_multiset_init_last")]
    pub fn init_last(&mut self) {
        unsafe { sys::gsl_multiset_init_last(self.unwrap_unique()) }
    }

    /// This function copies the elements of the multiset `self` into the multiset dest. The two
    /// multisets must have the same size.
    #[doc(alias = "gsl_multiset_memcpy")]
    pub fn copy(&self, dest: &mut MultiSet) -> Value {
        Value::from(unsafe { sys::gsl_multiset_memcpy(dest.unwrap_unique(), self.unwrap_shared()) })
    }

    /// This function returns the value of the i-th element of the multiset c. If i lies outside the
    /// allowed range of 0 to k-1 then the error handler is invoked and 0 is returned.
    #[doc(alias = "gsl_multiset_get")]
    pub fn get(&self, i: usize) -> usize {
        unsafe { sys::gsl_multiset_get(self.unwrap_shared(), i) }
    }

    /// This function returns the range (n) of the multiset `self`.
    #[doc(alias = "gsl_multiset_n")]
    pub fn n(&self) -> usize {
        unsafe { sys::gsl_multiset_n(self.unwrap_shared()) }
    }

    /// This function returns the number of elements (k) in the multiset `self`.
    #[doc(alias = "gsl_multiset_k")]
    pub fn k(&self) -> usize {
        unsafe { sys::gsl_multiset_k(self.unwrap_shared()) }
    }

    /// This function returns a pointer to the array of elements in the multiset `self`.
    #[doc(alias = "gsl_multiset_data")]
    pub fn data(&self) -> &[usize] {
        unsafe {
            let ptr = sys::gsl_multiset_data(self.unwrap_shared());
            ::std::slice::from_raw_parts(ptr, self.k())
        }
    }

    /// This function returns a pointer to the array of elements in the multiset `self`.
    #[doc(alias = "gsl_multiset_data")]
    pub fn data_mut(&mut self) -> &mut [usize] {
        unsafe {
            let ptr = sys::gsl_multiset_data(self.unwrap_shared());
            ::std::slice::from_raw_parts_mut(ptr, self.k())
        }
    }

    /// This function checks that the multiset self is valid. The k elements should lie in the range
    /// 0 to n-1, with each value occurring in non-decreasing order.
    #[doc(alias = "gsl_multiset_valid")]
    pub fn valid(&self) -> Value {
        // Little trick here: the function is expecting a mutable pointer whereas it doesn't need
        // to be...
        Value::from(unsafe { sys::gsl_multiset_valid(self.inner) })
    }

    /// This function advances the multiset self to the next multiset element in lexicographic order
    /// and returns [`Value::Success`]. If no further multisets elements are available it returns
    /// [`Value::Failure`] and leaves self unmodified. Starting with the first multiset and
    /// repeatedly applying this function will iterate through all possible multisets of a given
    /// order.
    #[doc(alias = "gsl_multiset_next")]
    pub fn next(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_multiset_next(self.unwrap_unique()) })
    }

    /// This function steps backwards from the multiset self to the previous multiset element in
    /// lexicographic order, returning [`Value::Success`]. If no previous multiset is available it
    /// returns [`Value::Failure`] and leaves self unmodified.
    #[doc(alias = "gsl_multiset_prev")]
    pub fn prev(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_multiset_prev(self.unwrap_unique()) })
    }

    pub fn print<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        write!(writer, "{:?}", self.data())
    }
}
