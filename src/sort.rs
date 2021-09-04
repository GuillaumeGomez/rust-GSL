//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Sorting

This chapter describes functions for sorting data, both directly and indirectly (using an index). All the functions use the heapsort algorithm.
Heapsort is an O(N \log N) algorithm which operates in-place and does not require any additional storage. It also provides consistent performance,
the running time for its worst-case (ordered data) being not significantly longer than the average and best cases. Note that the heapsort algorithm
does not preserve the relative ordering of equal elements—it is an unstable sort. However the resulting order of equal elements will be consistent
across different platforms when using these functions.

## References and Further Reading

The subject of sorting is covered extensively in Knuth’s Sorting and Searching,

Donald E. Knuth, The Art of Computer Programming: Sorting and Searching (Vol 3, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.
The Heapsort algorithm is described in the following book,

Robert Sedgewick, Algorithms in C, Addison-Wesley, ISBN 0201514257.
!*/

/// The following functions will sort the elements of an array or vector, either directly or indirectly. They are defined for all real and
/// integer types using the normal suffix rules. For example, the float versions of the array functions are gsl_sort_float and gsl_sort_float_index.
/// The corresponding vector functions are gsl_sort_vector_float and gsl_sort_vector_float_index. The prototypes are available in the header files
/// gsl_sort_float.h gsl_sort_vector_float.h. The complete set of prototypes can be included using the header files gsl_sort.h and gsl_sort_vector.h.
///
/// There are no functions for sorting complex arrays or vectors, since the ordering of complex numbers is not uniquely defined. To sort a complex
/// vector by magnitude compute a real vector containing the magnitudes of the complex elements, and sort this vector indirectly. The resulting index
/// gives the appropriate ordering of the original complex vector.
pub mod vectors {
    use crate::Value;
    use ffi::FFI;
    use types::{Permutation, VectorF64};

    /// This function sorts the n elements of the array data with stride stride into ascending numerical order.
    #[doc(alias = "gsl_sort")]
    pub fn sort(data: &mut [f64], stride: usize, n: usize) {
        unsafe { sys::gsl_sort(data.as_mut_ptr(), stride, n) }
    }

    /// This function sorts the n elements of the array data1 with stride stride1 into ascending numerical order, while making the same rearrangement
    /// of the array data2 with stride stride2, also of size n.
    #[doc(alias = "gsl_sort2")]
    pub fn sort2(data1: &mut [f64], stride1: usize, data2: &mut [f64], stride2: usize, n: usize) {
        unsafe { sys::gsl_sort2(data1.as_mut_ptr(), stride1, data2.as_mut_ptr(), stride2, n) }
    }

    /// This function sorts the elements of the vector v into ascending numerical order.
    #[doc(alias = "gsl_sort_vector")]
    pub fn sort_vector(v: &mut VectorF64) {
        unsafe { sys::gsl_sort_vector(v.unwrap_unique()) }
    }

    /// This function sorts the elements of the vector v1 into ascending numerical order, while making the same rearrangement of the vector v2.
    #[doc(alias = "gsl_sort_vector2")]
    pub fn sort_vector2(v1: &mut VectorF64, v2: &mut VectorF64) {
        unsafe { sys::gsl_sort_vector2(v1.unwrap_unique(), v2.unwrap_unique()) }
    }

    /// This function indirectly sorts the n elements of the array data with stride stride into ascending order, storing the resulting
    /// permutation in p. The array p must be allocated with a sufficient length to store the n elements of the permutation. The elements of p
    /// give the index of the array element which would have been stored in that position if the array had been sorted in place. The array data is not changed.
    #[doc(alias = "gsl_sort_index")]
    pub fn sort_index(p: &mut [usize], data: &[f64], stride: usize, n: usize) {
        unsafe { sys::gsl_sort_index(p.as_mut_ptr(), data.as_ptr(), stride, n) }
    }

    /// This function indirectly sorts the elements of the vector v into ascending order, storing the resulting permutation in p. The elements of p give the
    /// index of the vector element which would have been stored in that position if the vector had been sorted in place. The first element of p gives the index
    /// of the least element in v, and the last element of p gives the index of the greatest element in v. The vector v is not changed.
    #[doc(alias = "gsl_sort_vector_index")]
    pub fn sort_vector_index(p: &mut Permutation, v: &VectorF64) -> Value {
        Value::from(unsafe { sys::gsl_sort_vector_index(p.unwrap_unique(), v.unwrap_shared()) })
    }
}

/// The functions described in this section select the k smallest or largest elements of a data set of size N. The routines use an O(kN) direct insertion
/// algorithm which is suited to subsets that are small compared with the total size of the dataset. For example, the routines are useful for selecting the
/// 10 largest values from one million data points, but not for selecting the largest 100,000 values. If the subset is a significant part of the total dataset
/// it may be faster to sort all the elements of the dataset directly with an O(N \log N) algorithm and obtain the smallest or largest values that way.
pub mod select {
    use crate::Value;
    use ffi::FFI;
    use types::VectorF64;

    /// This function copies the k smallest elements of the array src, of size n and stride stride, in ascending numerical order into the array dest. The size
    /// k of the subset must be less than or equal to n. The data src is not modified by this operation.
    #[doc(alias = "gsl_sort_smallest")]
    pub fn sort_smallest(dest: &mut [f64], k: usize, src: &[f64], stride: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_smallest(dest.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as _)
        })
    }

    /// This function copies the k largest elements of the array src, of size n and stride stride, in descending numerical order into the array dest. k must
    /// be less than or equal to n. The data src is not modified by this operation.
    #[doc(alias = "gsl_sort_largest")]
    pub fn sort_largest(dest: &mut [f64], k: usize, src: &[f64], stride: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_largest(dest.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as _)
        })
    }

    /// This function copies the k smallest or largest elements of the vector v into the array dest. k must be less than or equal to the length of the vector v.
    #[doc(alias = "gsl_sort_vector_smallest")]
    pub fn sort_vector_smallest(dest: &mut [f64], k: usize, v: &VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_vector_smallest(dest.as_mut_ptr(), k, v.unwrap_shared())
        })
    }

    /// This function copies the k smallest or largest elements of the vector v into the array dest. k must be less than or equal to the length of the vector v.
    #[doc(alias = "gsl_sort_vector_largest")]
    pub fn sort_vector_largest(dest: &mut [f64], k: usize, v: &VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_vector_largest(dest.as_mut_ptr(), k, v.unwrap_shared())
        })
    }

    /// This function stores the indices of the k smallest elements of the array src, of size n and stride stride, in the array p. The indices are chosen so that
    /// the corresponding data is in ascending numerical order. k must be less than or equal to n. The data src is not modified by this operation.
    #[doc(alias = "gsl_sort_smallest_index")]
    pub fn sort_smallest_index(p: &mut [usize], k: usize, src: &[f64], stride: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_smallest_index(p.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as _)
        })
    }

    /// This function stores the indices of the k largest elements of the array src, of size n and stride stride, in the array p. The indices are chosen so that
    /// the corresponding data is in descending numerical order. k must be less than or equal to n. The data src is not modified by this operation.
    #[doc(alias = "gsl_sort_largest_index")]
    pub fn sort_largest_index(p: &mut [usize], k: usize, src: &[f64], stride: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_largest_index(p.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as _)
        })
    }

    /// This function stores the indices of the k smallest or largest elements of the vector v in the array p. k must be less than or equal to the length of
    /// the vector v.
    #[doc(alias = "gsl_sort_vector_smallest_index")]
    pub fn sort_vector_smallest_index(p: &mut [usize], k: usize, v: &VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_vector_smallest_index(p.as_mut_ptr(), k, v.unwrap_shared())
        })
    }

    /// This function stores the indices of the k smallest or largest elements of the vector v in the array p. k must be less than or equal to the length of
    /// the vector v.
    #[doc(alias = "gsl_sort_vector_largest_index")]
    pub fn sort_vector_largest_index(p: &mut [usize], k: usize, v: &VectorF64) -> Value {
        Value::from(unsafe {
            sys::gsl_sort_vector_largest_index(p.as_mut_ptr(), k, v.unwrap_shared())
        })
    }
}
