//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![macro_use]

#[doc(hidden)]
macro_rules! rgsl_error(
    ($msg:expr, $err_value:expr) => (
        {
            let file = file!();

            unsafe {
                let c_msg = ::std::ffi::CString::new($msg.as_bytes()).unwrap();
                let c_file = ::std::ffi::CString::new(file.as_bytes()).unwrap();

                let v = $err_value.into();
                sys::gsl_error(c_msg.as_ptr(), c_file.as_ptr(), line!() as i32, v)
            }
        }
    );
);

#[doc(hidden)]
macro_rules! ffi_wrap {
    ($name:tt, $cast:tt) => {
        unsafe { ffi::FFI::wrap(sys::$name as *mut sys::$cast) }
    };
}

#[doc(hidden)]
macro_rules! result {
    ($value:expr, $wrap:expr) => {{
        let ret = ::Value::from($value);
        match ret {
            ::Value::Success => Ok($wrap),
            e => Err(e),
        }
    }};
    ($value:expr) => {{
        let ret = ::Value::from($value);
        match ret {
            ::Value::Success => Ok(()),
            e => Err(e),
        }
    }};
}

#[doc(hidden)]
macro_rules! wrap_callback {
    ($f:expr, $F:ident) => {{
        unsafe extern "C" fn trampoline<F: Fn(f64) -> f64>(x: f64, params: *mut c_void) -> f64 {
            let f: &F = &*(f as *const F);
            f(x)
        }

        let f: Box<$F> = Box::new($f);
        sys::gsl_function_struct {
            function: transmute(trampoline::<$F> as usize),
            params: Box::into_raw(f) as *mut _,
        }
    }};
}
