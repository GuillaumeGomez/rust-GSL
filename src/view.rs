//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use std::marker::PhantomData;
use std::ops::Deref;

/// A struct which binds a type to a lifetime and prevent mutable access.
pub struct View<'a, T> {
    inner: T,
    phantom: PhantomData<&'a ()>,
}

impl<T> View<'_, T> {
    pub(crate) fn new<P>(inner: *mut P) -> Self
    where
        T: FFI<P>,
    {
        Self {
            inner: FFI::soft_wrap(inner),
            phantom: PhantomData,
        }
    }
}

impl<T> Deref for View<'_, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}
