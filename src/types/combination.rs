//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Combinations

This chapter describes functions for creating and manipulating combinations. A combination c is
represented by an array of k integers in the  range 0 to n-1, where each value c_i occurs at most
once. The combination c corresponds to indices of k elements chosen from an n element  vector.
Combinations are useful for iterating over all k-element subsets of a set.

## References and Further Reading

Further information on combinations can be found in,

Donald L. Kreher, Douglas R. Stinson, Combinatorial Algorithms: Generation, Enumeration and Search,
1998, CRC Press LLC, ISBN 084933988X
!*/

use crate::Value;
use ffi::FFI;
use std::fmt::{self, Debug, Formatter};

ffi_wrapper!(Combination, *mut sys::gsl_combination, gsl_combination_free);

impl Combination {
    /// This function allocates memory for a new combination with parameters n, k. The combination
    /// is not initialized and its elements are undefined. Use the function
    /// `Combination::new_init_first` if you want to create a combination which is initialized to
    /// the lexicographically first combination. A null pointer is returned if insufficient memory
    /// is available to create the combination.
    #[doc(alias = "gsl_combination_alloc")]
    pub fn new(n: usize, k: usize) -> Option<Self> {
        let tmp = unsafe { sys::gsl_combination_alloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function allocates memory for a new combination with parameters n, k and initializes it
    /// to the lexicographically first combination. A null pointer is returned if insufficient
    /// memory is available to create the combination.
    #[doc(alias = "gsl_combination_calloc")]
    pub fn new_with_init(n: usize, k: usize) -> Option<Self> {
        let tmp = unsafe { sys::gsl_combination_calloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function initializes the combination c to the lexicographically first combination, i.e.
    /// (0,1,2,...,k-1).
    #[doc(alias = "gsl_combination_init_first")]
    pub fn init_first(&mut self) {
        unsafe { sys::gsl_combination_init_first(self.unwrap_unique()) }
    }

    /// This function initializes the combination c to the lexicographically last combination, i.e.
    /// (n-k,n-k+1,…,n-1).
    #[doc(alias = "gsl_combination_init_last")]
    pub fn init_last(&mut self) {
        unsafe { sys::gsl_combination_init_last(self.unwrap_unique()) }
    }

    /// This function copies the elements of the combination self into the combination dest. The two
    /// combinations must have the same size.
    #[doc(alias = "gsl_combination_memcpy")]
    pub fn copy(&self, dest: &mut Combination) -> Value {
        Value::from(unsafe {
            sys::gsl_combination_memcpy(dest.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function returns the value of the i-th element of the combination self. If i lies
    /// outside the allowed range of 0 to k-1 then the error handler is invoked and 0 is returned.
    #[doc(alias = "gsl_combination_get")]
    pub fn get(&self, i: usize) -> usize {
        unsafe { sys::gsl_combination_get(self.unwrap_shared(), i) }
    }

    /// This function returns the range (n) of the combination self.
    #[doc(alias = "gsl_combination_n")]
    pub fn n(&self) -> usize {
        unsafe { sys::gsl_combination_n(self.unwrap_shared()) }
    }

    /// This function returns the number of elements (k) in the combination self.
    #[doc(alias = "gsl_combination_k")]
    pub fn k(&self) -> usize {
        unsafe { sys::gsl_combination_k(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_combination_data")]
    pub fn as_slice(&self) -> &[usize] {
        unsafe {
            let data = sys::gsl_combination_data(self.unwrap_shared());
            ::std::slice::from_raw_parts(data, self.k())
        }
    }

    #[doc(alias = "gsl_combination_data")]
    pub fn as_mut_slice(&mut self) -> &mut [usize] {
        unsafe {
            let data = sys::gsl_combination_data(self.unwrap_shared());
            ::std::slice::from_raw_parts_mut(data, self.k())
        }
    }

    /// This function checks that the combination self is valid. The k elements should lie in the
    /// range 0 to n-1, with each value occurring once at most and in increasing order.
    #[doc(alias = "gsl_combination_valid")]
    pub fn is_valid(&self) -> Value {
        // Little hack because `gsl_combination_valid` doesn't in fact need a mutable object...
        Value::from(unsafe { sys::gsl_combination_valid(self.inner) })
    }

    /// This function advances the combination self to the next combination in lexicographic order
    /// and returns `Success`. If no further combinations are available it returns Failure and
    /// leaves self unmodified. Starting with the first combination and repeatedly applying this
    /// function will iterate through all possible combinations of a given order.
    #[doc(alias = "gsl_combination_next")]
    pub fn next(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_combination_next(self.unwrap_unique()) })
    }

    /// This function steps backwards from the combination self to the previous combination in
    /// lexicographic order, returning `Success`. If no previous combination is available it returns
    /// `Failure` and leaves self unmodified.
    #[doc(alias = "gsl_combination_prev")]
    pub fn prev(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_combination_prev(self.unwrap_unique()) })
    }
}

impl Debug for Combination {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "[");
        let s = self.as_slice();
        for (pos, v) in s.iter().enumerate() {
            if pos == 0 {
                write!(f, "{}", v).unwrap();
            } else {
                write!(f, " {}", v).unwrap();
            }
        }
        write!(f, "]")
    }
}
