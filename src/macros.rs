//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![macro_use]

#[doc(hidden)]
macro_rules! ffi_wrap {
    ($name:tt, $cast:tt) => {
        unsafe { ffi::FFI::wrap(sys::$name as *mut sys::$cast) }
    };
}

#[doc(hidden)]
macro_rules! wrap_callback {
    ($f:expr, $F:ident) => {{
        unsafe extern "C" fn trampoline<F: Fn(f64) -> f64>(
            x: f64,
            params: *mut ::libc::c_void,
        ) -> f64 {
            let f: &F = &*(params as *const F);
            f(x)
        }

        let f: Box<$F> = Box::new($f);
        sys::gsl_function_struct {
            function: unsafe { ::std::mem::transmute(trampoline::<$F> as usize) },
            params: Box::into_raw(f) as *mut _,
        }
    }};
}
