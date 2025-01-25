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
    use crate::ffi::FFI;
    use crate::types::{Permutation, VectorF64};
    use crate::vector::{self, check_equal_len, Vector, VectorMut};
    use crate::Error;

    /// This function sorts the elements of the array `data` into
    /// ascending numerical order.
    ///
    /// # Examples
    ///
    /// ```
    /// use rgsl::sort::vectors::sort;
    /// let mut data = [4., 1., 3., 2.];
    /// sort(&mut data);
    /// assert_eq!(data, [1., 2., 3., 4.]);
    /// ```
    ///
    /// The same function can also be used with GSL vectors:
    ///
    /// ```
    /// use rgsl::{vector::VectorF64, sort::vectors::sort};
    /// let mut data = VectorF64::from_slice(&[4., 1., 3., 2.]).unwrap();
    /// sort(&mut data);
    /// assert_eq!(data.as_slice().unwrap(), [1., 2., 3., 4.]);
    /// ```
    #[doc(alias = "gsl_sort")]
    pub fn sort<T>(data: &mut T)
    where
        T: VectorMut<f64> + ?Sized,
    {
        unsafe { sys::gsl_sort(vector::as_mut_ptr(data), T::stride(data), T::len(data)) }
    }

    /// This function sorts the elements of the array `data1` into
    /// ascending numerical order, while making the same rearrangement
    /// of the array `data2`.  Panic if `data1` and `data2` do not
    /// have the same length.
    ///
    /// # Example
    ///
    /// ```
    /// use rgsl::sort::vectors::sort2;
    /// let mut data1 = [4., 1., 3., 2.];
    /// let mut data2 = [10., 20., 30., 40.];
    /// sort2(&mut data1, &mut data2);
    /// assert_eq!(data1, [1., 2., 3., 4.]);
    /// assert_eq!(data2, [20., 40., 30., 10.]);
    /// ```
    ///
    /// The same function can also be used with GSL vectors:
    ///
    /// ```
    /// use rgsl::{vector::VectorF64, sort::vectors::sort2};
    /// let mut data1 = VectorF64::from_slice(&[4., 1., 3., 2.]).unwrap();
    /// let mut data2 = VectorF64::from_slice(&[10., 20., 30., 40.]).unwrap();
    /// sort2(&mut data1, &mut data2);
    /// assert_eq!(data1.as_slice().unwrap(), [1., 2., 3., 4.]);
    /// assert_eq!(data2.as_slice().unwrap(), [20., 40., 30., 10.]);
    /// ```
    #[doc(alias = "gsl_sort2")]
    pub fn sort2<T1, T2>(data1: &mut T1, data2: &mut T2)
    where
        T1: VectorMut<f64> + ?Sized,
        T2: VectorMut<f64> + ?Sized,
    {
        check_equal_len(data1, data2)
            .expect("rgsl::sort::sort2: the vectors must have the same length");
        unsafe {
            sys::gsl_sort2(
                vector::as_mut_ptr(data1),
                T1::stride(data1),
                vector::as_mut_ptr(data2),
                T2::stride(data2),
                T1::len(data1),
            )
        }
    }

    /// This function sorts the elements of the vector v into ascending numerical order.
    #[doc(alias = "gsl_sort_vector")]
    #[deprecated(since = "8.0.0", note = "Please use `sort` instead")]
    pub fn sort_vector(v: &mut VectorF64) {
        unsafe { sys::gsl_sort_vector(v.unwrap_unique()) }
    }

    /// This function sorts the elements of the vector v1 into ascending numerical order, while making the same rearrangement of the vector v2.
    #[doc(alias = "gsl_sort_vector2")]
    #[deprecated(since = "8.0.0", note = "Please use `sort2` instead")]
    pub fn sort_vector2(v1: &mut VectorF64, v2: &mut VectorF64) {
        unsafe { sys::gsl_sort_vector2(v1.unwrap_unique(), v2.unwrap_unique()) }
    }

    /// This function indirectly sorts the elements of the array
    /// `data` into ascending order, storing the resulting permutation
    /// in `p`.  The slice `p` must have the same length as `data`.
    /// The elements of `p` give the index of the array element which
    /// would have been stored in that position if the array had been
    /// sorted in place.
    #[doc(alias = "gsl_sort_index")]
    pub fn sort_index<T>(p: &mut [usize], data: &T)
    where
        T: Vector<f64> + ?Sized,
    {
        if p.len() != T::len(data) {
            panic!("rgsl::sort::vectors::sort_index: `p` and `data` must have the same length");
        }
        unsafe {
            sys::gsl_sort_index(
                p.as_mut_ptr(),
                vector::as_ptr(data),
                T::stride(data),
                T::len(data),
            )
        }
    }

    /// This function indirectly sorts the elements of the vector v into ascending order, storing the resulting permutation in p. The elements of p give the
    /// index of the vector element which would have been stored in that position if the vector had been sorted in place. The first element of p gives the index
    /// of the least element in v, and the last element of p gives the index of the greatest element in v. The vector v is not changed.
    #[doc(alias = "gsl_sort_vector_index")]
    pub fn sort_vector_index(p: &mut Permutation, v: &VectorF64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_sort_vector_index(p.unwrap_unique(), v.unwrap_shared()) };
        Error::handle(ret, ())
    }
}

/// The functions described in this section select the k smallest or largest elements of a data set of size N. The routines use an O(kN) direct insertion
/// algorithm which is suited to subsets that are small compared with the total size of the dataset. For example, the routines are useful for selecting the
/// 10 largest values from one million data points, but not for selecting the largest 100,000 values. If the subset is a significant part of the total dataset
/// it may be faster to sort all the elements of the dataset directly with an O(N \log N) algorithm and obtain the smallest or largest values that way.
pub mod select {
    use crate::ffi::FFI;
    use crate::types::VectorF64;
    use crate::vector::{self, Vector};
    use crate::Error;

    /// This function copies the `dest.len()` smallest elements of the
    /// array `src`, in ascending numerical order into the array
    /// `dest`.  Panic if `dest.len()` is larger than the size of `src`.
    #[doc(alias = "gsl_sort_smallest")]
    pub fn sort_smallest<T>(dest: &mut [f64], src: &T) -> Result<(), Error>
    where
        T: Vector<f64> + ?Sized,
    {
        if dest.len() > T::len(src) {
            panic!("rgsl::sort::select::sort_smallest: `dest.len() > src.len()`");
        }
        let ret = unsafe {
            sys::gsl_sort_smallest(
                dest.as_mut_ptr(),
                dest.len(),
                vector::as_ptr(src),
                T::stride(src),
                T::len(src),
            )
        };
        Error::handle(ret, ())
    }

    /// This function copies the `dest.len()` largest elements of the
    /// array `src` in descending numerical order into the array
    /// `dest`.  Panic if `dest.len()` is larger than the size of `src`.
    #[doc(alias = "gsl_sort_largest")]
    pub fn sort_largest<T>(dest: &mut [f64], src: &T) -> Result<(), Error>
    where
        T: Vector<f64> + ?Sized,
    {
        if dest.len() > T::len(src) {
            panic!("rgsl::sort::select::sort_largest: `dest.len() > src.len()`");
        }
        let ret = unsafe {
            sys::gsl_sort_largest(
                dest.as_mut_ptr(),
                dest.len(),
                vector::as_ptr(src),
                T::stride(src),
                T::len(src),
            )
        };
        Error::handle(ret, ())
    }

    /// This function copies the `dest.len()` smallest elements of the
    /// vector `v` into the slice `dest`.  Panic if `dest.len()` is
    /// larger than the size of `src`.
    #[doc(alias = "gsl_sort_vector_smallest")]
    #[deprecated(since = "8.0.0", note = "Please use `sort_smallest` instead")]
    pub fn sort_vector_smallest(dest: &mut [f64], v: &VectorF64) -> Result<(), Error> {
        if dest.len() > v.len() {
            panic!("rgsl::sort::select::sort_vector_smallest: `dest.len() > v.len()`");
        }
        let ret = unsafe {
            sys::gsl_sort_vector_smallest(dest.as_mut_ptr(), dest.len(), v.unwrap_shared())
        };
        Error::handle(ret, ())
    }

    /// This function copies the `dest.len()` largest elements of the
    /// vector `v` into the array dest.  Panic if `dest.len()` is
    /// larger than the size of `src`.
    #[doc(alias = "gsl_sort_vector_largest")]
    #[deprecated(since = "8.0.0", note = "Please use `sort_largest` instead")]
    pub fn sort_vector_largest(dest: &mut [f64], v: &VectorF64) -> Result<(), Error> {
        if dest.len() > v.len() {
            panic!("rgsl::sort::select::sort_vector_largest: `dest.len() > v.len()`");
        }
        let ret = unsafe {
            sys::gsl_sort_vector_largest(dest.as_mut_ptr(), dest.len(), v.unwrap_shared())
        };
        Error::handle(ret, ())
    }

    /// This function stores the indices of the `p.len()` smallest
    /// elements of the vector `src` in the slice `p`.  The indices are
    /// chosen so that the corresponding data is in ascending
    /// numerical order.  Panic if `p.len()` is larger than the size
    /// of `src`.
    #[doc(alias = "gsl_sort_smallest_index")]
    pub fn sort_smallest_index<T>(p: &mut [usize], src: &T) -> Result<(), Error>
    where
        T: Vector<f64> + ?Sized,
    {
        if p.len() > T::len(src) {
            panic!("rgsl::sort::select::sort_smallest_index: `p.len() > src.len()`");
        }
        let ret = unsafe {
            sys::gsl_sort_smallest_index(
                p.as_mut_ptr(),
                p.len(),
                vector::as_ptr(src),
                T::stride(src),
                T::len(src),
            )
        };
        Error::handle(ret, ())
    }

    /// This function stores the indices of the `p.len()` largest
    /// elements of the vector `src` in the slice `p`.  The indices
    /// are chosen so that the corresponding data is in descending
    /// numerical order.  Panic if `p.len()` is larger than the size
    /// of `src`.
    #[doc(alias = "gsl_sort_largest_index")]
    pub fn sort_largest_index<T>(p: &mut [usize], src: &T) -> Result<(), Error>
    where
        T: Vector<f64> + ?Sized,
    {
        let ret = unsafe {
            sys::gsl_sort_largest_index(
                p.as_mut_ptr(),
                p.len(),
                vector::as_ptr(src),
                T::stride(src),
                T::len(src),
            )
        };
        Error::handle(ret, ())
    }

    /// This function stores the indices of the `p.len()` smallest
    /// elements of the vector `v` in the slice `p`.  Panic if
    /// `p.len()` is larger than the size of `src`.
    #[doc(alias = "gsl_sort_vector_smallest_index")]
    #[deprecated(since = "8.0.0", note = "Please use `sort_smallest_index` instead")]
    pub fn sort_vector_smallest_index(p: &mut [usize], v: &VectorF64) -> Result<(), Error> {
        if p.len() > v.len() {
            panic!("rgsl::sort::select::sort_vector_smallest_index: `p.len() > v.len()`");
        }
        let ret = unsafe {
            sys::gsl_sort_vector_smallest_index(p.as_mut_ptr(), p.len(), v.unwrap_shared())
        };
        Error::handle(ret, ())
    }

    /// This function stores the indices of the `p.len()` largest
    /// elements of the vector `v` in the slice `p`.  Panic if
    /// `p.len()` is larger than the size of `src`.
    #[doc(alias = "gsl_sort_vector_largest_index")]
    #[deprecated(since = "8.0.0", note = "Please use `sort_largest_index` instead")]
    pub fn sort_vector_largest_index(p: &mut [usize], v: &VectorF64) -> Result<(), Error> {
        if p.len() > v.len() {
            panic!("rgsl::sort::select::sort_vector_largest_index: `p.len() > v.len()`");
        }
        let ret = unsafe {
            sys::gsl_sort_vector_largest_index(p.as_mut_ptr(), p.len(), v.unwrap_shared())
        };
        Error::handle(ret, ())
    }
}
