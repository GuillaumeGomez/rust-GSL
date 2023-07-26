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
    ($f:expr, $F:ident $(+ $lt:lifetime)?) => {{
        unsafe extern "C" fn trampoline<$($lt,)? F: Fn(f64) -> f64 $( + $lt)?>(
            x: f64,
            params: *mut std::os::raw::c_void,
        ) -> f64 {
            let f: &F = &*(params as *const F);
            let x = f(x);
            x
        }

        sys::gsl_function_struct {
            function: Some(trampoline::<$F>),
            params: &$f as *const _ as *mut _,
        }
    }};
}

#[doc(hidden)]
macro_rules! ffi_wrapper {
    ($name:ident $(<$($lt:lifetime),*>)?, *mut $ty:ty, $drop:ident $(;$extra_id:ident: $extra_ty:ty => $extra_expr:expr;)* $(, $doc:expr)?) => {
        ffi_wrapper!($name $(<$($lt),*>)?, *mut $ty $(;$extra_id: $extra_ty => $extra_expr;)* $(, $doc)?);

        impl$(<$($lt),*>)? Drop for $name$(<$($lt),*>)? {
            fn drop(&mut self) {
                unsafe { sys::$drop(self.inner) };
                self.inner = std::ptr::null_mut();
            }
        }
    };
    ($name:ident $(<$($lt:lifetime),*>)?, *mut $ty:ty $(;$extra_id:ident: $extra_ty:ty => $extra_expr:expr;)* $(, $doc:expr)?) => {
        $(#[doc = $doc])?
        pub struct $name$(<$($lt),*>)? {
            inner: *mut $ty,
            $($extra_id: $extra_ty,)*
        }

        impl$(<$($lt),*>)? FFI<$ty> for $name$(<$($lt),*>)? {
            fn wrap(inner: *mut $ty) -> Self {
                Self { inner $(, $extra_id: $extra_expr)* }
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
    ($name:ident $(<$($lt:lifetime),*>)?, *const $ty:ty $(;$extra_id:ident: $extra_ty:ty => $extra_expr:expr;)* $(, $doc:expr)?) => {
        $(#[doc = $doc])?
        #[derive(Clone, Copy)]
        pub struct $name$(<$($lt),*>)? {
            inner: *const $ty,
            $($extra_id: $extra_ty,)*
        }

        impl$(<$($lt),*>)? FFI<$ty> for $name$(<$($lt),*>)? {
            fn wrap(inner: *mut $ty) -> Self {
                Self { inner $(, $extra_id: $extra_expr)* }
            }

            fn soft_wrap(inner: *mut $ty) -> Self {
                Self { inner $(, $extra_id: $extra_expr)* }
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
    () => {}
}

#[doc(hidden)]
macro_rules! result_handler {
    ($ret:ident, $value:expr) => {{
        if $ret == crate::sys::GSL_SUCCESS {
            Ok($value)
        } else {
            Err(crate::Value::from($ret))
        }
    }};
}
