//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
##Sorting

This chapter describes functions for sorting data, both directly and indirectly (using an index). All the functions use the heapsort algorithm. 
Heapsort is an O(N \log N) algorithm which operates in-place and does not require any additional storage. It also provides consistent performance, 
the running time for its worst-case (ordered data) being not significantly longer than the average and best cases. Note that the heapsort algorithm 
does not preserve the relative ordering of equal elements—it is an unstable sort. However the resulting order of equal elements will be consistent 
across different platforms when using these functions.
!*/

/// The following function provides a simple alternative to the standard library function qsort. It is intended for systems lacking qsort, not
/// as a replacement for it. The function qsort should be used whenever possible, as it will be faster and can provide stable ordering of equal
/// elements. Documentation for qsort is available in the GNU C Library Reference Manual.
/// 
/// The functions described in this section are defined in the header file gsl_heapsort.h.
pub mod objects {
    use ffi;

    /// This function sorts the count elements of the array array, each of size size, into ascending order using the comparison function compare.
    /// The type of the comparison function is defined by,
    /// 
    /// ```Rust
    /// i32 (*comparison_fn<T>) (a: &[T], b: &[T])
    /// ```
    /// 
    /// A comparison function should return a negative integer if the first argument is less than the second argument, 0 if the two arguments
    /// are equal and a positive integer if the first argumentis greater than the second argument.
    /// 
    /// For example, the following function can be used to sort doubles into ascending numerical order.
    /// 
    /// ```Rust
    /// fn compare_doubles(a: &[f64], b: &[f64]) -> i32 {
    ///    if (a[0] > b[0])
    ///       return 1;
    ///    else if (a[0] < b[0])
    ///       return -1;
    ///    else
    ///       return 0;
    /// }
    /// ```
    /// 
    /// The appropriate function call to perform the sort is,
    /// 
    /// ```Rust
    /// heapsort(array, compare_doubles);
    /// ```
    /// 
    /// Note that unlike qsort the heapsort algorithm cannot be made into a stable sort by pointer arithmetic. The trick of comparing pointers
    /// for equal elements in the comparison function does not work for the heapsort algorithm. The heapsort algorithm performs an internal
    /// rearrangement of the data which destroys its initial ordering.
    pub fn heapsort<T>(array: &mut[T], compare: compare_fn) {
        unsafe { ffi::gsl_heapsort(array.as_mut_ptr() as *mut c_void, array.len() as u64, ::std::mem::size_of::<T>() as u64, compare) }
    }

    /// This function indirectly sorts the count elements of the array array, each of size size, into ascending order using the comparison
    /// function compare. The resulting permutation is stored in p, an array of length n. The elements of p give the index of the array element
    /// which would have been stored in that position if the array had been sorted in place. The first element of p gives the index of the
    /// least element in array, and the last element of p gives the index of the greatest element in array. The array itself is not changed.
    pub fn heapsort_index<T>(p: &mut[u64], array: &[T], compare: compare_fn) -> enums::Value {
        unsafe { ffi::gsl_heapsort_index(p.as_mut_ptr(), array.as_ptr() as *const c_void, array.len() as u64, ::std::mem::size_of::<T>() as u64, compare) }
    }
}

/// The following functions will sort the elements of an array or vector, either directly or indirectly. They are defined for all real and
/// integer types using the normal suffix rules. For example, the float versions of the array functions are gsl_sort_float and gsl_sort_float_index.
/// The corresponding vector functions are gsl_sort_vector_float and gsl_sort_vector_float_index. The prototypes are available in the header files
/// gsl_sort_float.h gsl_sort_vector_float.h. The complete set of prototypes can be included using the header files gsl_sort.h and gsl_sort_vector.h.
/// 
/// There are no functions for sorting complex arrays or vectors, since the ordering of complex numbers is not uniquely defined. To sort a complex
/// vector by magnitude compute a real vector containing the magnitudes of the complex elements, and sort this vector indirectly. The resulting index
/// gives the appropriate ordering of the original complex vector.
pub mod vectors {
    use ffi;
    use types::VectorF64;

    /// This function sorts the n elements of the array data with stride stride into ascending numerical order.
    pub fn sort(data: &mut [f64], stride: u64) {
        unsafe { ffi::gsl_sort(data.as_mut_ptr(), stride, data.len() as u64) }
    }

    /// This function sorts the n elements of the array data1 with stride stride1 into ascending numerical order, while making the same rearrangement
    /// of the array data2 with stride stride2, also of size n.
    pub fn sort2(data1: &mut [f64], stride1: u64, data2: &mut [f64], stride2: u64) {
        unsafe { ffi::gsl_sort2(data1.as_mut_ptr(), stride1, data2.as_mut_ptr(), stride2, data1.len() as u64) }
    }

    /// This function sorts the elements of the vector v into ascending numerical order.
    pub fn sort_vector(v: &VectorF64) {
        unsafe { ffi::gsl_sort_vector(ffi::FFI::unwrap(v)) }
    }

    /// This function sorts the elements of the vector v1 into ascending numerical order, while making the same rearrangement of the vector v2.
    pub fn sort_vector2(v1: &VectorF64, v2: &VectorF64) {
        unsafe { ffi::gsl_sort_vector2(ffi::FFI::unwrap(v1), ffi::FFI::unwrap(v2)) }
    }

    /// This function indirectly sorts the n elements of the array data with stride stride into ascending order, storing the resulting
    /// permutation in p. The array p must be allocated with a sufficient length to store the n elements of the permutation. The elements of p
    /// give the index of the array element which would have been stored in that position if the array had been sorted in place. The array data is not changed.
    pub fn sort_index(p: &mut [u64], data: &[f64], stride: u64) {
        unsafe { ffi::gsl_sort_index(p.as_mut_ptr(), data.as_ptr(), stride, data.len() as u64) }
    }
}