//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Combinations

This chapter describes functions for creating and manipulating combinations. A combination c is
represented by an array of k integers in the  range 0 to n-1, where each value c_i occurs at most
once. The combination c corresponds to indices of k elements chosen from an n element  vector.
Combinations are useful for iterating over all k-element subsets of a set.

##References and Further Reading

Further information on combinations can be found in,

Donald L. Kreher, Douglas R. Stinson, Combinatorial Algorithms: Generation, Enumeration and Search,
1998, CRC Press LLC, ISBN 084933988X
!*/

use ffi;
use enums;
use std::fmt;
use std::fmt::{Formatter, Debug};
use c_vec::CSlice;

pub struct Combination {
    c: *mut ffi::gsl_combination,
    data: CSlice<usize>
}

impl Combination {
    /// This function allocates memory for a new combination with parameters n, k. The combination
    /// is not initialized and its elements are undefined. Use the function
    /// `Combination::new_init_first` if you want to create a combination which is initialized to
    /// the lexicographically first combination. A null pointer is returned if insufficient memory
    /// is available to create the combination.
    pub fn new(n: usize, k: usize) -> Option<Combination> {
        let tmp = unsafe { ffi::gsl_combination_alloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                if !(*tmp).data.is_null() {
                    Some(Combination {
                        c: tmp,
                        data: CSlice::new((*tmp).data, (*tmp).k as usize)
                    })
                } else {
                    Some(Combination {
                        c: tmp,
                        data: CSlice::new(tmp as *mut usize, 0usize)
                    })
                }
            }
        }
    }

    /// This function allocates memory for a new combination with parameters n, k and initializes it
    /// to the lexicographically first combination. A null pointer is returned if insufficient
    /// memory is available to create the combination.
    pub fn new_init_first(n: usize, k: usize) -> Option<Combination> {
        let tmp = unsafe { ffi::gsl_combination_calloc(n, k) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                if !(*tmp).data.is_null() {
                    Some(Combination {
                        c: tmp,
                        data: CSlice::new((*tmp).data, (*tmp).k as usize)
                    })
                } else {
                    Some(Combination {
                        c: tmp,
                        data: CSlice::new(tmp as *mut usize, 0usize)
                    })
                }
            }
        }
    }

    /// This function initializes the combination c to the lexicographically first combination, i.e.
    /// (0,1,2,...,k-1).
    pub fn init_first(&mut self) {
        unsafe { ffi::gsl_combination_init_first(self.c) }
    }

    /// This function initializes the combination c to the lexicographically last combination, i.e.
    /// (n-k,n-k+1,…,n-1).
    pub fn init_last(&mut self) {
        unsafe { ffi::gsl_combination_init_last(self.c) }
    }

    /// This function copies the elements of the combination self into the combination dest. The two
    /// combinations must have the same size.
    pub fn copy(&self, dest: &mut Combination) -> enums::Value {
        unsafe { ffi::gsl_combination_memcpy(dest.c, self.c) }
    }

    /// This function returns the value of the i-th element of the combination self. If i lies
    /// outside the allowed range of 0 to k-1 then the error handler is invoked and 0 is returned.
    pub fn get(&self, i: usize) -> usize {
        unsafe { ffi::gsl_combination_get(self.c, i) }
    }

    /// This function returns the range (n) of the combination self.
    pub fn n(&self) -> usize {
        unsafe { ffi::gsl_combination_n(self.c) }
    }

    /// This function returns the number of elements (k) in the combination self.
    pub fn k(&self) -> usize {
        unsafe { ffi::gsl_combination_k(self.c) }
    }

    /// This function returns a pointer to the array of elements in the combination self.
    pub fn as_slice<'r>(&'r self) -> &'r [usize] {
        self.data.as_ref()
    }

    /// This function returns a pointer to the array of elements in the combination self.
    pub fn as_mut_slice<'r>(&'r mut self) -> &'r mut [usize] {
        self.data.as_mut()
    }

    /// This function checks that the combination self is valid. The k elements should lie in the
    /// range 0 to n-1, with each value occurring once at most and in increasing order.
    pub fn is_valid(&self) -> enums::Value {
        unsafe { ffi::gsl_combination_valid(self.c) }
    }

    /// This function advances the combination self to the next combination in lexicographic order
    /// and returns `Success`. If no further combinations are available it returns Failure and
    /// leaves self unmodified. Starting with the first combination and repeatedly applying this
    /// function will iterate through all possible combinations of a given order.
    pub fn next(&mut self) -> enums::Value {
        unsafe { ffi::gsl_combination_next(self.c) }
    }

    /// This function steps backwards from the combination self to the previous combination in
    /// lexicographic order, returning `Success`. If no previous combination is available it returns
    /// `Failure` and leaves self unmodified.
    pub fn prev(&mut self) -> enums::Value {
        unsafe { ffi::gsl_combination_prev(self.c) }
    }
}

impl Drop for Combination {
    fn drop(&mut self) {
        unsafe { ffi::gsl_combination_free(self.c) };
        self.c = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_combination> for Combination {
    fn wrap(c: *mut ffi::gsl_combination) -> Combination {
        unsafe {
            Combination {
                c: c,
                data: CSlice::new((*c).data, (*c).k as usize)
            }
        }
    }

    fn soft_wrap(r: *mut ffi::gsl_combination) -> Combination {
        Self::wrap(r)
    }

    fn unwrap_shared(c: &Combination) -> *const ffi::gsl_combination {
        c.c as *const _
    }

    fn unwrap_unique(c: &mut Combination) -> *mut ffi::gsl_combination {
        c.c
    }
}

impl Debug for Combination {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "[");
        for tmp in 0..self.data.len() {
            if tmp == 0 {
                write!(f, "{}", self.data.get(tmp).unwrap());
            } else {
                write!(f, ", {}", self.data.get(tmp).unwrap());
            }
        }
        write!(f, "]")
    }
}