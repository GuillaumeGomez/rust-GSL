/// Utilities for interfacing with GSL/C
use std::ffi::CString;
use std::fs::OpenOptions;
use std::io;
use std::os::raw::c_char;
use std::path::Path;

use libc::{fclose, fopen, FILE};

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Write,
    Read,
}

/// A wrapper to handle I/O operations between GSL and rust
pub struct IOStream {
    inner: *mut FILE,
    mode: Mode,
}

impl IOStream {
    /// Open a file in write mode.
    pub fn fwrite_handle<P: AsRef<Path>>(file: &P) -> io::Result<IOStream> {
        {
            OpenOptions::new().write(true).open(file)?;
        }
        let path = CString::new(file.as_ref().to_str().unwrap()).unwrap();
        unsafe {
            Ok(IOStream {
                inner: fopen(path.as_ptr(), b"w\0".as_ptr() as *const c_char),
                mode: Mode::Write,
            })
        }
    }

    pub fn write_mode(&self) -> bool {
        self.mode == Mode::Write
    }

    #[doc(hidden)]
    pub fn as_raw(&mut self) -> *mut FILE {
        self.inner
    }
}

impl ::std::ops::Drop for IOStream {
    fn drop(&mut self) {
        unsafe {
            fclose(self.inner);
            self.inner = ::std::ptr::null_mut();
        }
    }
}
