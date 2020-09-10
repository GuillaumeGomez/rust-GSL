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
#[macro_export]
macro_rules! ffi_wrap {
    ($name:tt, $cast:tt) => {
        unsafe { ffi::FFI::wrap(sys::$name as *mut sys::$cast) }
    };
}
