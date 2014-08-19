//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Sorting

This chapter describes functions for sorting data, both directly and indirectly (using an index). All the functions use the heapsort algorithm. 
Heapsort is an O(N \log N) algorithm which operates in-place and does not require any additional storage. It also provides consistent performance, 
the running time for its worst-case (ordered data) being not significantly longer than the average and best cases. Note that the heapsort algorithm 
does not preserve the relative ordering of equal elements—it is an unstable sort. However the resulting order of equal elements will be consistent 
across different platforms when using these functions.

##References and Further Reading

The subject of sorting is covered extensively in Knuth’s Sorting and Searching,

Donald E. Knuth, The Art of Computer Programming: Sorting and Searching (Vol 3, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896850.
The Heapsort algorithm is described in the following book,

Robert Sedgewick, Algorithms in C, Addison-Wesley, ISBN 0201514257.
!*/

/// The following function provides a simple alternative to the standard library function qsort. It is intended for systems lacking qsort, not
/// as a replacement for it. The function qsort should be used whenever possible, as it will be faster and can provide stable ordering of equal
/// elements. Documentation for qsort is available in the GNU C Library Reference Manual.
/// 
/// The functions described in this section are defined in the header file gsl_heapsort.h.
pub mod objects {
    use enums;

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
    pub fn heapsort<T>(array: &mut[T], compare: ::comparison_fn<T>) {
        if array.len() == 0 {
            return;
        }
        let mut n = array.len() as u64 - 1;
        let mut k = n / 2;
        k += 1;
        while k > 0 {
            k -= 1;
            downheap(array, n, k, compare);
        }
        while n > 0 {
            array.swap(0u, n as uint);
            n -= 1;
            downheap(array, n, 0, compare);
        }
    }

    fn downheap<T>(array: &mut[T], n: u64, t_k: u64, compare: ::comparison_fn<T>) {
        let mut k = t_k;

        while k <= n / 2 {
            let mut j = 2 * k;

            if j < n && compare(&array[j as uint], &array[j as uint + 1]) < 0 {
                j += 1;
            }
            if compare(&array[k as uint], &array[j as uint]) < 0 {
                array.swap(j as uint, k as uint);
            } else {
                break;
            }
            k = j;
        }
    }

    /// This function indirectly sorts the count elements of the array array, each of size size, into ascending order using the comparison
    /// function compare. The resulting permutation is stored in p, an array of length n. The elements of p give the index of the array element
    /// which would have been stored in that position if the array had been sorted in place. The first element of p gives the index of the
    /// least element in array, and the last element of p gives the index of the greatest element in array. The array itself is not changed.
    pub fn heapsort_index<T>(p: &mut[u64], array: &[T], compare: ::comparison_fn<T>) -> enums::Value {
        if array.len() == 0 {
            return enums::Success;
        }
        for tmp in range(0u64, array.len() as u64) {
            p[tmp as uint] = tmp;
        }
        let mut n = array.len() as u64 - 1;
        let mut k = n / 2;
        k += 1;
        while k > 0 {
            k -= 1;
            downheap_index(p, array, n, k, compare);
        }
        while n > 0 {
            p.swap(0u, n as uint);
            n -= 1;
            downheap_index(p, array, n, 0, compare);
        }
        enums::Success
    }

    fn downheap_index<T>(p: &mut[u64], array: &[T], n: u64, t_k: u64, compare: ::comparison_fn<T>) {
        let pki = p[t_k as uint];
        let mut k = t_k;

        while k <= n / 2 {
            let mut j = 2 * k;

            if j < n && compare(&array[p[j as uint] as uint], &array[p[j as uint + 1] as uint]) < 0 {
                j += 1;
            }
            if compare(&array[pki as uint], &array[p[j as uint] as uint]) >= 0 {
                break;
            }
            p[k as uint] = p[j as uint];
            k = j;
        }
        p[k as uint] = pki;
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
    use types::{VectorF64, Permutation};
    use enums;

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

    /// This function indirectly sorts the elements of the vector v into ascending order, storing the resulting permutation in p. The elements of p give the
    /// index of the vector element which would have been stored in that position if the vector had been sorted in place. The first element of p gives the index
    /// of the least element in v, and the last element of p gives the index of the greatest element in v. The vector v is not changed.
    pub fn sort_vector_index(p: &Permutation, v: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_sort_vector_index(ffi::FFI::unwrap(p), ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }
}

/// The functions described in this section select the k smallest or largest elements of a data set of size N. The routines use an O(kN) direct insertion
/// algorithm which is suited to subsets that are small compared with the total size of the dataset. For example, the routines are useful for selecting the
/// 10 largest values from one million data points, but not for selecting the largest 100,000 values. If the subset is a significant part of the total dataset
/// it may be faster to sort all the elements of the dataset directly with an O(N \log N) algorithm and obtain the smallest or largest values that way.
pub mod select {
    use ffi;
    use types::VectorF64;
    use enums;

    /// This function copies the k smallest elements of the array src, of size n and stride stride, in ascending numerical order into the array dest. The size
    /// k of the subset must be less than or equal to n. The data src is not modified by this operation.
    pub fn sort_smallest(dest: &mut [f64], k: u64, src: &[f64], stride: u64) -> enums::Value {
        unsafe { ffi::gsl_sort_smallest(dest.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as u64) }
    }

    /// This function copies the k largest elements of the array src, of size n and stride stride, in descending numerical order into the array dest. k must
    /// be less than or equal to n. The data src is not modified by this operation.
    pub fn sort_largest(dest: &mut [f64], k: u64, src: &[f64], stride: u64) -> enums::Value {
        unsafe { ffi::gsl_sort_largest(dest.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as u64) }
    }

    /// This function copies the k smallest or largest elements of the vector v into the array dest. k must be less than or equal to the length of the vector v.
    pub fn sort_vector_smallest(dest: &mut [f64], k: u64, v: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_sort_vector_smallest(dest.as_mut_ptr(), k, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }

    /// This function copies the k smallest or largest elements of the vector v into the array dest. k must be less than or equal to the length of the vector v.
    pub fn sort_vector_largest(dest: &mut [f64], k: u64, v: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_sort_vector_largest(dest.as_mut_ptr(), k, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }

    /// This function stores the indices of the k smallest elements of the array src, of size n and stride stride, in the array p. The indices are chosen so that
    /// the corresponding data is in ascending numerical order. k must be less than or equal to n. The data src is not modified by this operation.
    pub fn sort_smallest_index(p: &mut [u64], k: u64, src: &[f64], stride: u64) -> enums::Value {
        unsafe { ffi::gsl_sort_smallest_index(p.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as u64) }
    }

    /// This function stores the indices of the k largest elements of the array src, of size n and stride stride, in the array p. The indices are chosen so that
    /// the corresponding data is in descending numerical order. k must be less than or equal to n. The data src is not modified by this operation.
    pub fn sort_largest_index(p: &mut [u64], k: u64, src: &[f64], stride: u64) -> enums::Value {
        unsafe { ffi::gsl_sort_largest_index(p.as_mut_ptr(), k, src.as_ptr(), stride, src.len() as u64) }
    }

    /// This function stores the indices of the k smallest or largest elements of the vector v in the array p. k must be less than or equal to the length of
    /// the vector v.
    pub fn sort_vector_smallest_index(p: &mut [u64], k: u64, v: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_sort_vector_smallest_index(p.as_mut_ptr(), k, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }

    /// This function stores the indices of the k smallest or largest elements of the vector v in the array p. k must be less than or equal to the length of
    /// the vector v.
    pub fn sort_vector_largest_index(p: &mut [u64], k: u64, v: &VectorF64) -> enums::Value {
        unsafe { ffi::gsl_sort_vector_largest_index(p.as_mut_ptr(), k, ffi::FFI::unwrap(v) as *const ffi::gsl_vector) }
    }
}