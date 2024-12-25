//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Vectors

Vectors are defined by a gsl_vector structure which describes a slice of a block. Different vectors can be created which point to the
same block. A vector slice is a set of equally-spaced elements of an area of memory.

The gsl_vector structure contains five components, the size, the stride, a pointer to the memory where the elements are stored, data, a
pointer to the block owned by the vector, block, if any, and an ownership flag, owner. The structure is very simple and looks like this,

```C
typedef struct
{
  size_t size;
  size_t stride;
  double * data;
  gsl_block * block;
  int owner;
} gsl_vector;
```

The size is simply the number of vector elements. The range of valid indices runs from 0 to size-1. The stride is the step-size from one
element to the next in physical memory, measured in units of the appropriate datatype. The pointer data gives the location of the first
element of the vector in memory. The pointer block stores the location of the memory block in which the vector elements are located (if
any). If the vector owns this block then the owner field is set to one and the block will be deallocated when the vector is freed. If the
vector points to a block owned by another object then the owner field is zero and any underlying block will not be deallocated with the
vector.
!*/

use crate::ffi::FFI;
use crate::Value;
use std::{
    fmt::{self, Debug, Formatter},
    marker::PhantomData,
    ops::Range,
};

use paste::paste;

#[cfg(feature = "complex")]
extern crate num_complex;
#[cfg(feature = "complex")]
use self::num_complex::Complex;

#[allow(clippy::len_without_is_empty)]
/// Trait implemented by types that are considered vectors by this crate.
/// Elements of the vector are of type `F` (`f32` or `f64`).
///
/// Bring this trait into scope in order to add methods to specify
/// strides to the types implementing `Vector`.
///
/// # Safety
/// One must make sore that `(len - 1) * stride` does not exceed the
/// length of the underlying slice.
pub unsafe trait Vector<F> {
    /// Return the number of elements in the vector.
    ///
    /// This is an associated function rather than a method to avoid
    /// conflicting with similarly named methods.
    fn len(x: &Self) -> usize;

    /// The distance in the slice between two consecutive elements of
    /// the vector in [`Vector::as_slice`] and [`VectorMut::as_mut_slice`].
    fn stride(x: &Self) -> usize;

    /// Return a reference to the underlying slice.  Note that the
    /// `i`th element of the vector, `0 <= i < len(x)`, is the
    /// `i * stride` element in the slice.
    ///
    /// This is an associated function rather than a method to avoid
    /// conflicting with similarly named methods.
    fn as_slice(x: &Self) -> &[F];

    fn slice(&self, r: Range<usize>) -> Option<Slice<'_, F>> {
        // The fields of `Slice` are not public, hence there is no way
        // to implement this function outside this module.  This
        // guarantee that there is no out-of-bounds access as soon as
        // the above methods are correctly implemented.
        let stride = Self::stride(self);
        let slice = Self::as_slice(self);
        // FIXME: use std::slice::SliceIndex methods when stable.
        if r.end == 0 {
            return Some(Slice {
                vec: &slice[0..0],
                len: 0,
                stride: 1,
            });
        }
        let end = (r.end - 1) * stride + 1;
        if r.start > r.end || end > slice.len() {
            None
        } else {
            let start = r.start * stride;
            Some(Slice {
                vec: &slice[start..end],
                len: r.end - r.start,
                stride,
            })
        }
    }
}

pub unsafe trait VectorMut<F>: Vector<F> {
    /// Same as [`Vector::as_slice`] but mutable.
    fn as_mut_slice(x: &mut Self) -> &mut [F];

    /// Same as [`Vector::slice`] but mutable.
    fn slice_mut(&mut self, r: Range<usize>) -> Option<SliceMut<'_, F>> {
        let stride = Self::stride(self);
        let slice = Self::as_mut_slice(self);
        // FIXME: use std::slice::SliceIndex methods when stable.
        if r.end == 0 {
            return Some(SliceMut {
                vec: &mut slice[0..0],
                len: 0,
                stride: 1,
            });
        }
        let end = (r.end - 1) * stride + 1;
        if r.start > r.end || end > slice.len() {
            None
        } else {
            let start = r.start * stride;
            Some(SliceMut {
                vec: &mut slice[start..end],
                len: r.end - r.start,
                stride,
            })
        }
    }
}

pub struct Slice<'a, F> {
    vec: &'a [F],
    len: usize,
    stride: usize,
}

pub struct SliceMut<'a, F> {
    vec: &'a mut [F],
    len: usize,
    stride: usize,
}

unsafe impl<'a, F> Vector<F> for Slice<'a, F> {
    #[inline]
    fn len(x: &Self) -> usize {
        x.len
    }
    #[inline]
    fn stride(x: &Self) -> usize {
        x.stride
    }
    #[inline]
    fn as_slice(x: &Self) -> &[F] {
        x.vec
    }
}

unsafe impl<'a, F> Vector<F> for SliceMut<'a, F> {
    #[inline]
    fn len(x: &Self) -> usize {
        x.len
    }
    #[inline]
    fn stride(x: &Self) -> usize {
        x.stride
    }
    #[inline]
    fn as_slice(x: &Self) -> &[F] {
        x.vec
    }
}

unsafe impl<'a, F> VectorMut<F> for SliceMut<'a, F> {
    #[inline]
    fn as_mut_slice(x: &mut Self) -> &mut [F] {
        x.vec
    }
}

/// Return the length of `x` as a `i32` value (to use in CBLAS calls).
#[inline]
pub(crate) fn len<F, T: Vector<F> + ?Sized>(x: &T) -> i32 {
    T::len(x).try_into().expect("Length must fit in `i32`")
}

#[inline]
pub(crate) fn as_ptr<F, T: Vector<F> + ?Sized>(x: &T) -> *const F {
    T::as_slice(x).as_ptr()
}

#[inline]
pub(crate) fn as_mut_ptr<F, T: VectorMut<F> + ?Sized>(x: &mut T) -> *mut F {
    T::as_mut_slice(x).as_mut_ptr()
}

/// Return the stride of `x` as a `i32` value (to use in CBLAS calls).
#[inline]
pub(crate) fn stride<F, T: Vector<F> + ?Sized>(x: &T) -> i32 {
    T::stride(x).try_into().expect("Stride must fit in `i32`")
}

#[inline]
pub(crate) fn check_equal_len<T1, T2, F>(x: &T1, y: &T2) -> Result<(), Value>
where
    T1: Vector<F> + ?Sized,
    T2: Vector<F> + ?Sized,
{
    if T1::len(x) != T2::len(y) {
        return Err(Value::Invalid);
    }
    Ok(())
}

macro_rules! gsl_vec {
    ($rust_name:ident, $name:ident, $rust_ty:ident) => (
paste! {

pub struct $rust_name {
    vec: *mut sys::$name,
    can_free: bool,
}

impl Drop for $rust_name {
    #[doc(alias = $name _free)]
    fn drop(&mut self) {
        if self.can_free {
            unsafe { sys::[<$name _free>](self.vec) };
            self.vec = std::ptr::null_mut();
        }
    }
}

impl Debug for $rust_name {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            write!(f, "<null>")
        } else {
            write!(f, "{:?}", self.as_slice().expect("conversion to slice failed"))
        }
    }
}

impl FFI<sys::$name> for $rust_name {
    fn wrap(vec: *mut sys::$name) -> Self {
        Self {
            vec,
            can_free: true,
        }
    }

    fn soft_wrap(vec: *mut sys::$name) -> Self {
        Self {
            vec,
            can_free: false,
        }
    }

    fn unwrap_shared(&self) -> *const sys::$name {
        self.vec as *const _
    }

    fn unwrap_unique(&mut self) -> *mut sys::$name {
        self.vec
    }
}

impl $rust_name {
    #[doc = "create a new " $rust_name " with all elements set to zero"]
    #[doc(alias = $name _calloc)]
    pub fn new(size: usize) -> Option<$rust_name> {
        let tmp = unsafe { sys::[<$name _calloc>](size) };

        if tmp.is_null() {
            None
        } else {
            Some($rust_name::wrap(tmp))
        }
    }

    #[doc(alias = $name _alloc)]
    pub fn from_slice(slice: &[$rust_ty]) -> Option<$rust_name> {
        let tmp = unsafe { sys::[<$name _alloc>](slice.len() as _) };

        if tmp.is_null() {
            None
        } else {
            let mut v = Self::wrap(tmp);

            for (pos, tmp) in slice.iter().enumerate() {
                v.set(pos as _, *tmp);
            }
            Some(v)
        }
    }

    pub fn len(&self) -> usize {
        let ptr = self.unwrap_shared();
        if ptr.is_null() {
            0
        } else {
            unsafe { (*ptr).size }
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn as_slice(&self) -> Option<&[$rust_ty]> {
        let ptr = unsafe { (*self.unwrap_shared()).data };
        if ptr.is_null() {
            None
        } else {
            Some(unsafe { std::slice::from_raw_parts(ptr, self.len()) })
        }
    }

    pub fn as_slice_mut(&mut self) -> Option<&mut [$rust_ty]> {
        let ptr = unsafe { (*self.unwrap_shared()).data };
        if ptr.is_null() {
            None
        } else {
            Some(unsafe { std::slice::from_raw_parts_mut(ptr, self.len()) })
        }
    }

    /// This function returns the i-th element of a vector v. If i lies outside the allowed range
    /// of 0 to n-1 then the error handler is invoked and 0 is returned.
    #[doc(alias = $name _get)]
    pub fn get(&self, i: usize) -> $rust_ty {
        unsafe { sys::[<$name _get>](self.unwrap_shared(), i) }
    }

    /// This function sets the value of the i-th element of a vector v to x. If i lies outside the
    /// allowed range of 0 to n-1 then the error handler is invoked.
    #[doc(alias = $name _set)]
    pub fn set(&mut self, i: usize, x: $rust_ty) -> &mut $rust_name {
        unsafe { sys::[<$name _set>](self.unwrap_unique(), i, x) };
        self
    }

    /// This function sets all the elements of the vector v to the value x.
    #[doc(alias = $name _set_all)]
    pub fn set_all(&mut self, x: $rust_ty) -> &mut $rust_name {
        unsafe { sys::[<$name _set_all>](self.unwrap_unique(), x) };
        self
    }

    /// This function sets all the elements of the vector v to zero.
    #[doc(alias = $name _set_zero)]
    pub fn set_zero(&mut self) -> &mut $rust_name {
        unsafe { sys::[<$name _set_zero>](self.unwrap_unique()) };
        self
    }

    /// This function makes a basis vector by setting all the elements of the vector v to zero
    /// except for the i-th element which is set to one.
    #[doc(alias = $name _set_basis)]
    pub fn set_basis(&mut self, i: usize) -> &mut $rust_name {
        unsafe { sys::[<$name _set_basis>](self.unwrap_unique(), i) };
        self
    }

    /// This function copies the elements of the other vector into the self vector. The two vectors
    /// must have the same length.
    #[doc(alias = $name _memcpy)]
    pub fn copy_from(&mut self, other: &$rust_name) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _memcpy>](
                self.unwrap_unique(),
                other.unwrap_shared()) };
        result_handler!(ret, ())
    }

    /// This function copies the elements of the self vector into the other vector. The two vectors
    /// must have the same length.
    #[doc(alias = $name _memcpy)]
    pub fn copy_to(&self, other: &mut $rust_name) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _memcpy>](other.unwrap_unique(), self.unwrap_shared()) };
        result_handler!(ret, ())
    }

    /// This function exchanges the elements of the vectors by copying. The two vectors must have
    /// the same length.
    #[doc(alias = $name _swap)]
    pub fn swap(&mut self, other: &mut $rust_name) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _swap>](other.unwrap_unique(), self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function exchanges the i-th and j-th elements of the vector v in-place.
    #[doc(alias = $name _swap_elements)]
    pub fn swap_elements(&mut self, i: usize, j: usize) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _swap_elements>](self.unwrap_unique(), i, j) };
        result_handler!(ret, ())
    }

    /// This function reverses the order of the elements of the vector v.
    #[doc(alias = $name _reverse)]
    pub fn reverse(&mut self) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _reverse>](self.unwrap_unique()) };
        result_handler!(ret, ())
    }

    /// This function adds the elements of the other vector to the elements of the self vector.
    /// The result a_i <- a_i + b_i is stored in self and other remains unchanged. The two vectors
    /// must have the same length.
    #[doc(alias = $name _add)]
    pub fn add(&mut self, other: &$rust_name) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _add>](self.unwrap_unique(), other.unwrap_shared()) };
        result_handler!(ret, ())
    }

    /// This function subtracts the elements of the self vector from the elements of the other
    /// vector. The result a_i <- a_i - b_i is stored in self and other remains unchanged. The two
    /// vectors must have the same length.
    #[doc(alias = $name _sub)]
    pub fn sub(&mut self, other: &$rust_name) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _sub>](self.unwrap_unique(), other.unwrap_shared()) };
        result_handler!(ret, ())
    }

    /// This function multiplies the elements of the self vector a by the elements of the other
    /// vector. The result `a_i <- a_i * b_i` is stored in self and other remains unchanged. The two
    /// vectors must have the same length.
    #[doc(alias = $name _mul)]
    pub fn mul(&mut self, other: &$rust_name) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _mul>](self.unwrap_unique(), other.unwrap_shared()) };
        result_handler!(ret, ())
    }

    /// This function divides the elements of the self vector by the elements of the other vector.
    /// The result `a_i <- a_i / b_i` is stored in self and other remains unchanged. The two vectors
    /// must have the same length.
    #[doc(alias = $name _div)]
    pub fn div(&mut self, other: &$rust_name) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _div>](self.unwrap_unique(), other.unwrap_shared()) };
        result_handler!(ret, ())
    }

    /// This function multiplies the elements of the self vector by the constant factor x. The
    /// result `a_i <- a_i` is stored in `self`.
    #[doc(alias = $name _scale)]
    pub fn scale(&mut self, x: $rust_ty) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _scale>](self.unwrap_unique(), x) };
        result_handler!(ret, ())
    }

    /// This function adds the constant value x to the elements of the self vector. The result
    /// `a_i <- a_i + x` is stored in `self`.
    #[doc(alias = $name _add_constant)]
    pub fn add_constant(&mut self, x: $rust_ty) -> Result<(), Value> {
        let ret = unsafe { sys::[<$name _add_constant>](self.unwrap_unique(), x) };
        result_handler!(ret, ())
    }

    /// This function returns the maximum value in the self vector.
    #[doc(alias = $name _max)]
    pub fn max(&self) -> $rust_ty {
        unsafe { sys::[<$name _max>](self.unwrap_shared()) }
    }

    /// This function returns the minimum value in the self vector.
    #[doc(alias = $name _min)]
    pub fn min(&self) -> $rust_ty {
        unsafe { sys::[<$name _min>](self.unwrap_shared()) }
    }

    /// This function returns the minimum and maximum values in the self vector.
    #[doc(alias = $name _minmax)]
    pub fn minmax(&self) -> ($rust_ty, $rust_ty) {
        let mut min_out = 0 as _;
        let mut max_out = 0 as _;

        unsafe {
            sys::[<$name _minmax>](self.unwrap_shared(), &mut min_out, &mut max_out);
        }
        (min_out, max_out)
    }

    /// This function returns the index of the maximum value in the self vector.
    /// When there are several equal maximum elements then the lowest index is returned.
    #[doc(alias = $name _max_index)]
    pub fn max_index(&self) -> usize {
        unsafe { sys::[<$name _max_index>](self.unwrap_shared()) }
    }

    /// This function returns the index of the minimum value in the self vector.
    /// When there are several equal minimum elements then the lowest index is returned.
    #[doc(alias = $name _min_index)]
    pub fn min_index(&self) -> usize {
        unsafe { sys::[<$name _min_index>](self.unwrap_shared()) }
    }

    /// This function returns the indices of the minimum and maximum values in the self vector.
    /// When there are several equal minimum or maximum elements then the lowest indices are
    /// returned.
    #[doc(alias = $name _minmax_index)]
    pub fn minmax_index(&self) -> (usize, usize) {
        let mut imin = 0;
        let mut imax = 0;

        unsafe { sys::[<$name _minmax_index>](self.unwrap_shared(), &mut imin, &mut imax) };
        (imin, imax)
    }

    /// This function returns true if all the elements of the self vector are equal to 0.
    #[doc(alias = $name _isnull)]
    pub fn is_null(&self) -> bool {
        unsafe { sys::[<$name _isnull>](self.unwrap_shared()) == 1 }
    }

    /// This function returns true if all the elements of the self vector are stricly positive.
    #[doc(alias = $name _ispos)]
    pub fn is_pos(&self) -> bool {
        unsafe { sys::[<$name _ispos>](self.unwrap_shared()) == 1 }
    }

    /// This function returns true if all the elements of the self vector are stricly negative.
    #[doc(alias = $name _isneg)]
    pub fn is_neg(&self) -> bool {
        unsafe { sys::[<$name _isneg>](self.unwrap_shared()) == 1 }
    }

    /// This function returns true if all the elements of the self vector are stricly non-negative.
    #[doc(alias = $name _isnonneg)]
    pub fn is_non_neg(&self) -> bool {
        unsafe { sys::[<$name _isnonneg>](self.unwrap_shared()) == 1 }
    }

    #[doc(alias = $name _equal)]
    pub fn equal(&self, other: &$rust_name) -> bool {
        unsafe { sys::[<$name _equal>](self.unwrap_shared(), other.unwrap_shared()) == 1 }
    }

    pub fn clone(&self) -> Option<$rust_name> {
        if self.unwrap_shared().is_null() {
            None
        } else {
            if let Some(mut v) = $rust_name::new(self.len()) {
                if v.copy_from(self).is_err() {
                    None
                } else {
                    Some(v)
                }
            } else {
                None
            }
        }
    }

    #[doc(alias = $name _subvector)]
    pub fn subvector(&mut self, offset: usize, n: usize) -> [<$rust_name View>]<'_> {
        [<$rust_name View>]::from_vector(self, offset, n)
    }
}

pub struct [<$rust_name View>]<'a> {
    v: sys::[<$name _view>],
    #[allow(dead_code)]
    phantom: PhantomData<&'a ()>,
}

impl<'a> [<$rust_name View>]<'a> {
    #[doc(hidden)]
    pub(crate) fn wrap<F: FnOnce(Option<Self>)>(v: sys::[<$name _view>], f: F) {
        let tmp = Self {
            v,
            phantom: PhantomData,
        };
        let is_none = {
            let v = &tmp.v.vector;
            let tmp = $rust_name::soft_wrap(v as *const _ as usize as *mut _);
            tmp.as_slice().is_none()
        };
        if is_none {
            f(None)
        } else {
            f(Some(tmp))
        }
    }

    /// These functions return a vector view of a subvector of another vector v. The start of the
    /// new vector is offset by offset elements from the start of the original vector. The new
    /// vector has n elements. Mathematically, the i-th element of the new vector v’ is given by,
    ///
    /// v'(i) = v->data[(offset + i)*v->stride]
    ///
    /// where the index i runs from 0 to n-1.
    ///
    /// The data pointer of the returned vector struct is set to null if the combined parameters
    /// (offset,n) overrun the end of the original vector.
    ///
    /// The new vector is only a view of the block underlying the original vector, v. The block
    /// containing the elements of v is not owned by the new vector. When the view goes out of scope
    /// the original vector v and its block will continue to exist. The original memory can only be
    /// deallocated by freeing the original vector. Of course, the original vector should not be
    /// deallocated while the view is still in use.
    ///
    /// The function gsl_vector_const_subvector is equivalent to gsl_vector_subvector but can be
    /// used for vectors which are declared const.
    #[doc(alias = $name _subvector)]
    pub fn from_vector(v: &'a mut $rust_name, offset: usize, n: usize) -> Self {
        unsafe {
            Self {
                v: sys::[<$name _subvector>](v.unwrap_unique(), offset, n),
                phantom: PhantomData,
            }
        }
    }

    /// These functions return a vector view of a subvector of another vector v with an additional
    /// stride argument. The subvector is formed in the same way as for gsl_vector_subvector but the
    /// new vector has n elements with a step-size of stride from one element to the next in the
    /// original vector. Mathematically, the i-th element of the new vector v’ is given by,
    ///
    /// v'(i) = v->data[(offset + i*stride)*v->stride]
    /// where the index i runs from 0 to n-1.
    ///
    /// Note that subvector views give direct access to the underlying elements of the original
    /// vector. For example, the following code will zero the even elements of the vector v of
    /// length n, while leaving the odd elements untouched,
    ///
    /// ```C
    /// gsl_vector_view v_even
    ///   = gsl_vector_subvector_with_stride (v, 0, 2, n/2);
    /// gsl_vector_set_zero (&v_even.vector);
    /// ```
    /// A vector view can be passed to any subroutine which takes a vector argument just as a
    /// directly allocated vector would be, using &view.vector.
    /// For example, the following code computes the norm of the odd elements of v using the BLAS
    /// routine DNRM2,
    ///
    /// ```C
    /// gsl_vector_view v_odd
    ///   = gsl_vector_subvector_with_stride (v, 1, 2, n/2);
    /// double r = gsl_blas_dnrm2 (&v_odd.vector);
    /// ```
    /// The function gsl_vector_const_subvector_with_stride is equivalent to
    /// gsl_vector_subvector_with_stride but can be used for vectors which are declared const.
    #[doc(alias = $name _subvector_with_stride)]
    pub fn from_vector_with_stride(
        v: &'a mut $rust_name,
        offset: usize,
        stride: usize,
        n: usize,
    ) -> Self {
        unsafe {
            Self {
                v: sys::[<$name _subvector_with_stride>](v.vec, offset, stride, n),
                phantom: PhantomData,
            }
        }
    }

    /// These functions return a vector view of an array. The start of the new vector is given by
    /// base and has n elements. Mathematically, the i-th element of the new vector v’ is given by,
    ///
    /// ```text
    /// v'(i) = base[i]
    /// ```
    ///
    /// where the index i runs from 0 to n-1.
    ///
    /// The array containing the elements of v is not owned by the new vector view. When the view
    /// goes out of scope the original array will continue to exist. The original memory can only be
    /// deallocated by freeing the original pointer base. Of course, the original array should not
    /// be deallocated while the view is still in use.
    ///
    /// The function gsl_vector_const_view_array is equivalent to gsl_vector_view_array but can be
    /// used for arrays which are declared const.
    #[doc(alias = $name _view_array)]
    pub fn from_array(base: &'a mut [f64]) -> Self {
        unsafe {
            Self {
                v: sys::[<$name _view_array>](base.as_mut_ptr() as _, base.len() as _),
                phantom: PhantomData,
            }
        }
    }

    /// These functions return a vector view of an array base with an additional stride argument.
    /// The subvector is formed in the same way as for gsl_vector_view_array but the new vector has
    /// n elements with a step-size of stride from one element to the next in the original
    /// array. Mathematically, the i-th element of the new vector v’ is given by,
    ///
    /// v'(i) = base[i*stride]
    ///
    /// where the index i runs from 0 to n-1.
    ///
    /// Note that the view gives direct access to the underlying elements of the original array. A
    /// vector view can be passed to any subroutine which takes a vector argument just as a directly
    /// allocated vector would be, using &view.vector.
    ///
    /// The function gsl_vector_const_view_array_with_stride is equivalent to
    /// gsl_vector_view_array_with_stride but can be used for arrays which are declared const.
    #[doc(alias = $name _view_array_with_stride)]
    pub fn from_array_with_stride(base: &'a mut [$rust_ty], stride: usize) -> Self {
        unsafe {
            Self {
                v: sys::[<$name _view_array_with_stride>](
                    base.as_mut_ptr(),
                    stride,
                    base.len() as _,
                ),
                phantom: PhantomData,
            }
        }
    }

    pub fn vector<F: FnOnce(Option<&$rust_name>)>(&self, f: F) {
        let v = &self.v.vector;
        let tmp = $rust_name::soft_wrap(v as *const _ as usize as *mut _);
        if tmp.as_slice().is_none() {
            f(None)
        } else {
            f(Some(&tmp))
        }
    }

    pub fn vector_mut<F: FnOnce(Option<&mut $rust_name>)>(&mut self, f: F) {
        let v = &mut self.v.vector;
        let mut tmp = $rust_name::soft_wrap(v as *mut _);
        if tmp.as_slice().is_none() {
            f(None)
        } else {
            f(Some(&mut tmp))
        }
    }
} // end of impl block

    unsafe impl Vector<$rust_ty> for $rust_name {
        #[inline]
        fn len(x: &Self) -> usize {
            $rust_name::len(x)
        }
        #[inline]
        fn stride(x: &Self) -> usize {
            let ptr = x.unwrap_shared();
            if ptr.is_null() {
                1
            } else {
                unsafe { (*ptr).stride }
            }
        }
        #[inline]
        fn as_slice(x: &Self) -> &[$rust_ty] {
            $rust_name::as_slice(x).unwrap_or(&[])
        }
    }
    unsafe impl VectorMut<$rust_ty> for $rust_name {
        #[inline]
        fn as_mut_slice(x: &mut Self) -> &mut [$rust_ty] {
            $rust_name::as_slice_mut(x).unwrap_or(&mut [])
        }
    }

} // end of paste! block
); // end of gsl_vec macro
}

gsl_vec!(VectorF32, gsl_vector_float, f32);
gsl_vec!(VectorF64, gsl_vector, f64);
gsl_vec!(VectorI32, gsl_vector_int, i32);
gsl_vec!(VectorU32, gsl_vector_uint, u32);

// Implement the `Vector` trait on standard vectors.

macro_rules! impl_AsRef {
    ($ty: ty) => {
        unsafe impl<T> Vector<$ty> for T
        where
            T: AsRef<[$ty]> + ?Sized,
        {
            #[inline]
            fn len(x: &Self) -> usize {
                x.as_ref().len()
            }
            #[inline]
            fn stride(_: &Self) -> usize {
                1
            }
            #[inline]
            fn as_slice(x: &Self) -> &[$ty] {
                x.as_ref()
            }
        }

        unsafe impl<T> VectorMut<$ty> for T
        where
            T: Vector<$ty> + AsMut<[$ty]> + ?Sized,
        {
            #[inline]
            fn as_mut_slice(x: &mut Self) -> &mut [$ty] {
                x.as_mut()
            }
        }
    };
}

impl_AsRef!(f32);
impl_AsRef!(f64);
#[cfg(feature = "complex")]
impl_AsRef!(Complex<f32>);
#[cfg(feature = "complex")]
impl_AsRef!(Complex<f64>);
