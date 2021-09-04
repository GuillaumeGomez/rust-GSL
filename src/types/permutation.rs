//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[cfg(feature = "v2_2")]
use crate::MatrixF64;
use crate::Value;
use crate::{MatrixComplexF32, MatrixComplexF64, MatrixF32, VectorF64};
use ffi::FFI;
use std::fmt;
use std::fmt::{Debug, Formatter};

ffi_wrapper!(Permutation, *mut sys::gsl_permutation, gsl_permutation_free);

/// ## Permutations in cyclic form
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
    #[doc(alias = "gsl_permutation_alloc")]
    pub fn new(n: usize) -> Option<Permutation> {
        let tmp = unsafe { sys::gsl_permutation_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function allocates memory for a new permutation of size n and initializes it to the identity. A null pointer is returned if
    /// insufficient memory is available to create the permutation.
    #[doc(alias = "gsl_permutation_calloc")]
    pub fn new_with_init(n: usize) -> Option<Permutation> {
        let tmp = unsafe { sys::gsl_permutation_calloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function initializes the permutation p to the identity, i.e. (0,1,2,…,n-1).
    #[doc(alias = "gsl_permutation_init")]
    pub fn init(&mut self) {
        unsafe { sys::gsl_permutation_init(self.unwrap_unique()) }
    }

    /// This function copies the elements of the permutation src into the permutation dest. The two permutations must have the same size.
    #[doc(alias = "gsl_permutation_memcpy")]
    pub fn copy(&self, dest: &mut Permutation) -> Value {
        Value::from(unsafe {
            sys::gsl_permutation_memcpy(dest.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function returns the value of the i-th element of the permutation p. If i lies outside the allowed range of 0 to n-1 then
    /// the error handler is invoked and 0 is returned.
    #[doc(alias = "gsl_permutation_get")]
    pub fn get(&self, i: usize) -> usize {
        unsafe { sys::gsl_permutation_get(self.unwrap_shared(), i) }
    }

    /// This function exchanges the i-th and j-th elements of the permutation p.
    #[doc(alias = "gsl_permutation_swap")]
    pub fn swap(&mut self, i: usize, j: usize) -> Value {
        Value::from(unsafe { sys::gsl_permutation_swap(self.unwrap_unique(), i, j) })
    }

    /// This function returns the size of the permutation p.
    #[doc(alias = "gsl_permutation_size")]
    pub fn size(&self) -> usize {
        unsafe { sys::gsl_permutation_size(self.unwrap_shared()) }
    }

    #[doc(alias = "gsl_permutation_data")]
    pub fn as_slice(&self) -> &[usize] {
        unsafe {
            let data = sys::gsl_permutation_data(self.unwrap_shared());
            ::std::slice::from_raw_parts(data, self.size())
        }
    }

    #[doc(alias = "gsl_permutation_data")]
    pub fn as_mut_slice(&mut self) -> &mut [usize] {
        unsafe {
            let data = sys::gsl_permutation_data(self.unwrap_shared());
            ::std::slice::from_raw_parts_mut(data, self.size())
        }
    }

    /// This function checks that the permutation p is valid. The n elements should contain each of the numbers 0 to n-1 once and only once.
    #[doc(alias = "gsl_permutation_valid")]
    pub fn is_valid(&self) -> bool {
        Value::from(unsafe { sys::gsl_permutation_valid(self.unwrap_shared()) }) == ::Value::Success
    }

    /// This function reverses the elements of the permutation p.
    #[doc(alias = "gsl_permutation_reverse")]
    pub fn reverse(&mut self) {
        unsafe { sys::gsl_permutation_reverse(self.unwrap_unique()) }
    }

    /// This function computes the inverse of the permutation p, storing the result in inv.
    #[doc(alias = "gsl_permutation_inverse")]
    pub fn inverse(&self, inv: &mut Permutation) -> Value {
        Value::from(unsafe {
            sys::gsl_permutation_inverse(inv.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function advances the permutation p to the next permutation in lexicographic order and returns GSL_SUCCESS. If no further
    /// permutations are available it returns GSL_FAILURE and leaves p unmodified. Starting with the identity permutation and repeatedly
    /// applying this function will iterate through all possible permutations of a given order.
    #[doc(alias = "gsl_permutation_next")]
    pub fn next(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_permutation_next(self.unwrap_unique()) })
    }

    /// This function steps backwards from the permutation p to the previous permutation in lexicographic order, returning GSL_SUCCESS.
    /// If no previous permutation is available it returns GSL_FAILURE and leaves p unmodified.
    #[doc(alias = "gsl_permutation_prev")]
    pub fn prev(&mut self) -> Value {
        Value::from(unsafe { sys::gsl_permutation_prev(self.unwrap_unique()) })
    }

    /// This function applies the permutation to the array data of size n with stride stride.
    #[doc(alias = "gsl_permutation_data")]
    pub fn permute(&mut self, data: &mut [f64], stride: usize) -> Value {
        Value::from(unsafe {
            let data_ptr = sys::gsl_permutation_data(self.unwrap_shared());
            sys::gsl_permute(data_ptr, data.as_mut_ptr(), stride, data.len() as _)
        })
    }

    /// This function applies the inverse of the permutation p to the array data of size n with stride stride.
    #[doc(alias = "gsl_permutation_data")]
    pub fn permute_inverse(&mut self, data: &mut [f64], stride: usize) -> Value {
        Value::from(unsafe {
            let data_ptr = sys::gsl_permutation_data(self.unwrap_shared());
            sys::gsl_permute_inverse(data_ptr, data.as_mut_ptr(), stride, data.len() as _)
        })
    }

    /// This function applies the permutation p to the elements of the vector v, considered as a row-vector acted on by a permutation
    /// matrix from the right, v' = v P. The j-th column of the permutation matrix P is given by the p_j-th column of the identity matrix.
    /// The permutation p and the vector v must have the same length.
    #[doc(alias = "gsl_permute_vector")]
    pub fn permute_vector(&mut self, v: &mut VectorF64) -> Value {
        Value::from(unsafe { sys::gsl_permute_vector(self.unwrap_unique(), v.unwrap_unique()) })
    }

    /// This function applies the inverse of the permutation p to the elements of the vector v, considered as a row-vector acted on by an inverse permutation
    /// matrix from the right, v' = v P^T. Note that for permutation matrices the inverse is the same as the transpose. The j-th column of the permutation
    /// matrix P is given by the p_j-th column of the identity matrix. The permutation p and the vector v must have the same length.
    #[doc(alias = "gsl_permute_vector_inverse")]
    pub fn permute_vector_inverse(&self, v: &mut VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_permute_vector_inverse(self.unwrap_shared(), v.unwrap_unique())
        })
    }

    #[cfg(feature = "v2_2")]
    #[cfg_attr(feature = "dox", doc(cfg(feature = "v2_2")))]
    #[doc(alias = "gsl_permute_matrix")]
    pub fn permute_matrix(&self, A: &mut MatrixF64) -> Value {
        Value::from(unsafe { sys::gsl_permute_matrix(self.unwrap_shared(), A.unwrap_unique()) })
    }

    #[doc(alias = "gsl_permute_matrix_float")]
    pub fn permute_matrix_float(&self, A: &mut MatrixF32) -> Value {
        Value::from(unsafe {
            sys::gsl_permute_matrix_float(self.unwrap_shared(), A.unwrap_unique())
        })
    }

    #[doc(alias = "gsl_permute_matrix_complex")]
    pub fn permute_matrix_complex(&self, A: &mut MatrixComplexF64) -> Value {
        Value::from(unsafe {
            sys::gsl_permute_matrix_complex(self.unwrap_shared(), A.unwrap_unique())
        })
    }

    #[doc(alias = "gsl_permute_matrix_complex_float")]
    pub fn permute_matrix_complex_float(&self, A: &mut MatrixComplexF32) -> Value {
        Value::from(unsafe {
            sys::gsl_permute_matrix_complex_float(self.unwrap_shared(), A.unwrap_unique())
        })
    }

    /// This function combines the two permutations pa and pb into a single permutation p, where p = pa * pb. The permutation p is equivalent to applying pb
    /// first and then pa.
    #[doc(alias = "gsl_permutation_mul")]
    pub fn mul(&mut self, pa: &Permutation, pb: &Permutation) -> Value {
        Value::from(unsafe {
            sys::gsl_permutation_mul(self.unwrap_unique(), pa.unwrap_shared(), pb.unwrap_shared())
        })
    }

    /// This function computes the canonical form of the permutation self and stores it in the output argument q.
    #[doc(alias = "gsl_permutation_linear_to_canonical")]
    pub fn linear_to_canonical(&self, q: &mut Permutation) -> Value {
        Value::from(unsafe {
            sys::gsl_permutation_linear_to_canonical(q.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function converts the self permutation in canonical form back into linear form storing it in the output argument p.
    #[doc(alias = "gsl_permutation_canonical_to_linear")]
    pub fn canonical_to_linear(&self, p: &mut Permutation) -> Value {
        Value::from(unsafe {
            sys::gsl_permutation_canonical_to_linear(p.unwrap_unique(), self.unwrap_shared())
        })
    }

    /// This function counts the number of inversions in the self permutation. An inversion is any pair of elements that are not in order. For example, the
    /// permutation 2031 has three inversions, corresponding to the pairs (2,0) (2,1) and (3,1). The identity permutation has no inversions.
    #[doc(alias = "gsl_permutation_inversions")]
    pub fn inversions(&self) -> usize {
        unsafe { sys::gsl_permutation_inversions(self.unwrap_shared()) }
    }

    /// This function counts the number of cycles in the self permutation, given in linear form.
    #[doc(alias = "gsl_permutation_linear_cycles")]
    pub fn linear_cycles(&self) -> usize {
        unsafe { sys::gsl_permutation_linear_cycles(self.unwrap_shared()) }
    }

    /// This function counts the number of cycles in the self permutation, given in canonical form.
    #[doc(alias = "gsl_permutation_canonical_cycles")]
    pub fn canonical_cycles(&self) -> usize {
        unsafe { sys::gsl_permutation_canonical_cycles(self.unwrap_shared()) }
    }
}

impl Debug for Permutation {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        if self.unwrap_shared().is_null() {
            write!(f, "<null>")
        } else {
            write!(f, "{:?}", self.as_slice())
        }
    }
}
