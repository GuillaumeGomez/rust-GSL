/// Utilities for interfacing with GSL/C

use std::ffi::CString;
use std::path::Path;
use std::fs::File;

use libc::{FILE, fopen, fclose};

/// A wrapper to handle I/O operations between GSL and rust
pub struct IOStream {
    inner: *mut FILE,
    mode: Mode,
}

#[allow(dead_code)]
enum Mode {
    Write,
    Read,
}

impl IOStream {
    /// Open a file in write mode.
    pub fn fwrite_handle(file: &Path) -> Result<IOStream, ::std::io::ErrorKind> {
        {
            let f = File::open(file);
            if f.is_err() {
                return Err(f.unwrap_err().kind());
            }
        }
        let w = CString::new("w").unwrap();
        let path = CString::new(file.to_str().unwrap()).unwrap();
        unsafe {
            Ok(IOStream {
                inner: fopen(path.as_ptr(), w.as_ptr()),
                mode: Mode::Write,
            })
        }
    }

    pub fn write_mode(&self) -> bool {
        match self.mode {
            Mode::Write => true,
            Mode::Read => false,
        }
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