//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use types::{VectorF64};
use ffi;
use enums;
use std::fmt;
use std::fmt::{Formatter, Debug};
use c_vec::CSlice;

pub struct Permutation {
    p: *mut ffi::gsl_permutation,
    d: CSlice<u64>
}

///##Permutations in cyclic form
/// 
/// A permutation can be represented in both linear and cyclic notations. The functions described in this section convert between the two forms.
/// The linear notation is an index mapping, and has already been described above. The cyclic notation expresses a permutation as a series of
/// circular rearrangements of groups of elements, or cycles.
/// 
/// For example, under the cycle (1 2 3), 1 is replaced by 2, 2 is replaced by 3 and 3 is replaced by 1 in a circular fashion. Cycles of different
/// sets of elements can be combined independently, for example (1 2 3) (4 5) combines the cycle (1 2 3) with the cycle (4 5), which is an exchange
/// of elements 4 and 5. A cycle of length one represents an element which is unchanged by the permutation and is referred to as a singleton.
/// 
/// It can be shown that every permutation can be decomposed into combinations of cycles. The decomposition is not unique, but can always be
/// rearranged into a standard canonical form by a reordering of elements. The library uses the canonical form defined in Knuth’s Art of Computer
/// Programming (Vol 1, 3rd Ed, 1997) Section 1.3.3, p.178.
/// 
/// The procedure for obtaining the canonical form given by Knuth is,
/// 
/// Write all singleton cycles explicitly
/// Within each cycle, put the smallest number first
/// Order the cycles in decreasing order of the first number in the cycle.
/// For example, the linear representation (2 4 3 0 1) is represented as (1 4) (0 2 3) in canonical form. The permutation corresponds to an exchange
/// of elements 1 and 4, and rotation of elements 0, 2 and 3.
/// 
/// The important property of the canonical form is that it can be reconstructed from the contents of each cycle without the brackets. In addition,
/// by removing the brackets it can be considered as a linear representation of a different permutation. In the example given above the permutation
/// (2 4 3 0 1) would become (1 4 0 2 3). This mapping has many applications in the theory of permutations.
impl Permutation {
    /// This function allocates memory for a new permutation of size n. The permutation is not initialized and its elements are undefined.
    /// Use the function gsl_permutation_calloc if you want to create a permutation which is initialized to the identity. A null pointer is
    /// returned if insufficient memory is available to create the permutation.
    pub fn new(n: u64) -> Option<Permutation> {
        let tmp = unsafe { ffi::gsl_permutation_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                Some(Permutation {
                    p: tmp,
                    d: CSlice::new((*tmp).data, (*tmp).size as usize)
                })
            }
        }
    }

    /// This function allocates memory for a new permutation of size n and initializes it to the identity. A null pointer is returned if
    /// insufficient memory is available to create the permutation.
    pub fn new_with_init(n: u64) -> Option<Permutation> {
        let tmp = unsafe { ffi::gsl_permutation_calloc(n) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                Some(Permutation {
                    p: tmp,
                    d: CSlice::new((*tmp).data, (*tmp).size as usize)
                })
            }
        }
    }

    /// This function initializes the permutation p to the identity, i.e. (0,1,2,…,n-1).
    pub fn init(&self) {
        unsafe { ffi::gsl_permutation_init(self.p) }
    }

    /// This function copies the elements of the permutation src into the permutation dest. The two permutations must have the same size.
    pub fn copy(&self, dest: &Permutation) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_memcpy(dest.p, self.p as *const ffi::gsl_permutation) }
    }

    /// This function returns the value of the i-th element of the permutation p. If i lies outside the allowed range of 0 to n-1 then
    /// the error handler is invoked and 0 is returned.
    pub fn get(&self, i: u64) -> u64 {
        unsafe { ffi::gsl_permutation_get(self.p as *const ffi::gsl_permutation, i) }
    }

    /// This function exchanges the i-th and j-th elements of the permutation p.
    pub fn swap(&self, i: u64, j: u64) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_swap(self.p, i, j) }
    }

    /// This function returns the size of the permutation p.
    pub fn size(&self) -> u64 {
        unsafe { ffi::gsl_permutation_size(self.p as *const ffi::gsl_permutation) }
    }

    /// This function returns a pointer to the array of elements in the permutation p.
    pub fn data<'r>(&'r mut self) -> &'r mut [u64] {
        self.d.as_mut_slice()
    }

    /// This function checks that the permutation p is valid. The n elements should contain each of the numbers 0 to n-1 once and only once.
    pub fn is_valid(&self) -> bool {
        match unsafe { ffi::gsl_permutation_valid(self.p as *const ffi::gsl_permutation) } {
            ::Value::Success => true,
            _ => false
        }
    }

    /// This function reverses the elements of the permutation p.
    pub fn reverse(&self) {
        unsafe { ffi::gsl_permutation_reverse(self.p) }
    }

    /// This function computes the inverse of the permutation p, storing the result in inv.
    pub fn inverse(&self, inv: &Permutation) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_inverse(inv.p, self.p as *const ffi::gsl_permutation) }
    }

    /// This function advances the permutation p to the next permutation in lexicographic order and returns GSL_SUCCESS. If no further
    /// permutations are available it returns GSL_FAILURE and leaves p unmodified. Starting with the identity permutation and repeatedly
    /// applying this function will iterate through all possible permutations of a given order.
    pub fn next(&self) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_next(self.p) }
    }

    /// This function steps backwards from the permutation p to the previous permutation in lexicographic order, returning GSL_SUCCESS.
    /// If no previous permutation is available it returns GSL_FAILURE and leaves p unmodified.
    pub fn prev(&self) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_prev(self.p) }
    }

    /// This function applies the permutation to the array data of size n with stride stride.
    pub fn permute(&self, data: &mut [f64], stride: u64) -> enums::value::Value {
        unsafe { ffi::gsl_permute((*self.p).data as *const u64, data.as_mut_ptr(), stride, data.len() as u64) }
    }

    /// This function applies the inverse of the permutation p to the array data of size n with stride stride.
    pub fn permute_inverse(&self, data: &mut [f64], stride: u64) -> enums::value::Value {
        unsafe { ffi::gsl_permute_inverse((*self.p).data as *const u64, data.as_mut_ptr(), stride, data.len() as u64) }
    }

    /// This function applies the permutation p to the elements of the vector v, considered as a row-vector acted on by a permutation
    /// matrix from the right, v' = v P. The j-th column of the permutation matrix P is given by the p_j-th column of the identity matrix.
    /// The permutation p and the vector v must have the same length.
    pub fn permute_vector(&self, v: &VectorF64) -> enums::value::Value {
        unsafe { ffi::gsl_permute_vector(self.p as *const ffi::gsl_permutation, ffi::FFI::unwrap(v)) }
    }

    /// This function applies the inverse of the permutation p to the elements of the vector v, considered as a row-vector acted on by an inverse permutation
    /// matrix from the right, v' = v P^T. Note that for permutation matrices the inverse is the same as the transpose. The j-th column of the permutation
    /// matrix P is given by the p_j-th column of the identity matrix. The permutation p and the vector v must have the same length.
    pub fn permute_vector_inverse(&self, v: &VectorF64) -> enums::value::Value {
        unsafe { ffi::gsl_permute_vector_inverse(self.p as *const ffi::gsl_permutation, ffi::FFI::unwrap(v)) }
    }

    /// This function combines the two permutations pa and pb into a single permutation p, where p = pa * pb. The permutation p is equivalent to applying pb
    /// first and then pa.
    pub fn mul(&self, pa: &Permutation, pb: &Permutation) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_mul(self.p, pa.p as *const ffi::gsl_permutation, pb.p as *const ffi::gsl_permutation) }
    }

    /// This function computes the canonical form of the permutation self and stores it in the output argument q.
    pub fn linear_to_canonical(&self, q: &Permutation) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_linear_to_canonical(q.p, self.p as *const ffi::gsl_permutation) }
    }

    /// This function converts the self permutation in canonical form back into linear form storing it in the output argument p.
    pub fn canonical_to_linear(&self, p: &Permutation) -> enums::value::Value {
        unsafe { ffi::gsl_permutation_canonical_to_linear(p.p, self.p as *const ffi::gsl_permutation) }
    }

    /// This function counts the number of inversions in the self permutation. An inversion is any pair of elements that are not in order. For example, the
    /// permutation 2031 has three inversions, corresponding to the pairs (2,0) (2,1) and (3,1). The identity permutation has no inversions.
    pub fn inversions(&self) -> u64 {
        unsafe { ffi::gsl_permutation_inversions(self.p as *const ffi::gsl_permutation) }
    }

    /// This function counts the number of cycles in the self permutation, given in linear form.
    pub fn linear_cycles(&self) -> u64 {
        unsafe { ffi::gsl_permutation_linear_cycles(self.p as *const ffi::gsl_permutation) }
    }

    /// This function counts the number of cycles in the self permutation, given in canonical form.
    pub fn canonical_cycles(&self) -> u64 {
        unsafe { ffi::gsl_permutation_canonical_cycles(self.p as *const ffi::gsl_permutation) }
    }
}

impl Drop for Permutation {
    fn drop(&mut self) {
        unsafe { ffi::gsl_permutation_free(self.p) };
        self.p = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_permutation> for Permutation {
    fn wrap(p: *mut ffi::gsl_permutation) -> Permutation {
        unsafe {
            Permutation {
                p: p,
                d: CSlice::new((*p).data, (*p).size as usize)
            }
        }
    }

    fn unwrap(p: &Permutation) -> *mut ffi::gsl_permutation {
        p.p
    }
}

impl Debug for Permutation {
    #[allow(unused_must_use)]
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "[");
        unsafe {
            for x in range(0u64, (*self.p).size) {
                let tmp = (*self.p).data.offset(x as isize);

                write!(f, " {}", *tmp);
            }
        }
        write!(f, "]")
    }
}