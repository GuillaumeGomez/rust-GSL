//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![macro_use]

macro_rules! rgsl_error(
    ($msg:expr, $err_value:expr) => (
        {
            let file = file!();

            unsafe {
                let c_msg = ::std::ffi::CString::from_slice($msg.as_bytes());
                let c_file = ::std::ffi::CString::from_slice(file.as_bytes());

                ffi::gsl_error(c_msg.as_ptr(), c_file.as_ptr(), line!() as i32, $err_value as i32)
            }
        }
    );
);