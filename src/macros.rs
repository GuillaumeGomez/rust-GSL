//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![macro_use]

#[doc(hidden)]
macro_rules! ffi_wrap {
    ($name:tt) => {
        unsafe { $crate::ffi::FFI::wrap(sys::$name as *mut _) }
    };
}

#[doc(hidden)]
macro_rules! wrap_callback {
    ($f:expr, $F:ident) => {{
        unsafe extern "C" fn trampoline<F: Fn(f64) -> f64>(
            x: f64,
            params: *mut ::std::os::raw::c_void,
        ) -> f64 {
            let f: &F = &*(params as *const F);
            f(x)
        }

        sys::gsl_function_struct {
            function: Some(trampoline::<$F>),
            params: &$f as *const _ as *mut _,
        }
    }};
}

#[doc(hidden)]
macro_rules! ffi_wrapper {
    ($name:ident, *mut $ty:ty, $drop:ident $(, $doc:expr)?) => {
        ffi_wrapper!($name, *mut $ty $(, $doc)?);

        impl Drop for $name {
            fn drop(&mut self) {
                unsafe { sys::$drop(self.inner) };
                self.inner = ::std::ptr::null_mut();
            }
        }
    };
    ($name:ident, *mut $ty:ty $(, $doc:expr)?) => {
        $(#[doc = $doc])?
        pub struct $name {
            inner: *mut $ty,
        }

        impl FFI<$ty> for $name {
            fn wrap(inner: *mut $ty) -> Self {
                Self { inner }
            }

            fn soft_wrap(r: *mut $ty) -> Self {
                Self::wrap(r)
            }

            #[inline]
            fn unwrap_shared(&self) -> *const $ty {
                self.inner as *const _
            }

            #[inline]
            fn unwrap_unique(&mut self) -> *mut $ty {
                self.inner
            }
        }
    };
    ($name:ident, *const $ty:ty $(, $doc:expr)?) => {
        $(#[doc = $doc])?
        #[derive(Clone, Copy)]
        pub struct $name {
            inner: *const $ty,
        }

        impl FFI<$ty> for $name {
            fn wrap(inner: *mut $ty) -> Self {
                Self { inner }
            }

            fn soft_wrap(inner: *mut $ty) -> Self {
                Self { inner }
            }

            #[inline]
            fn unwrap_shared(&self) -> *const $ty {
                self.inner
            }

            fn unwrap_unique(&mut self) -> *mut $ty {
                unimplemented!()
            }
        }
    };
}
