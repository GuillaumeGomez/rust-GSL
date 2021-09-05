//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# N-tuples

This chapter describes functions for creating and manipulating ntuples, sets of values associated with events. The ntuples are stored in
files. Their values can be extracted in any combination and booked in a histogram using a selection function.

The values to be stored are held in a user-defined data structure, and an ntuple is created associating this data structure with a file.
The values are then written to the file (normally inside a loop) using the ntuple functions described below.

A histogram can be created from ntuple data by providing a selection function and a value function. The selection function specifies
whether an event should be included in the subset to be analyzed or not. The value function computes the entry to be added to the histogram
for each event.

## Histogramming ntuple values

Once an ntuple has been created its contents can be histogrammed in various ways using the function gsl_ntuple_project. Two user-defined
functions must be provided, a function to select events and a function to compute scalar values. The selection function and the value
function both accept the ntuple row as a first argument and other parameters as a second argument.

The selection function determines which ntuple rows are selected for histogramming.

## References and Further Reading

Further information on the use of ntuples can be found in the documentation for the CERN packages PAW and HBOOK (available online).
!*/

use crate::Value;
use ffi::FFI;
use std::ffi::CString;
use std::mem::MaybeUninit;
use std::os::raw::{c_char, c_void};
use std::path::Path;

pub struct WriteNTuples {
    n: *mut sys::gsl_ntuple,
}

impl WriteNTuples {
    /// This function creates a new write-only ntuple file filename for ntuples of size size and
    /// returns a pointer to the newly created ntuple struct. Any existing file with the same name
    /// is truncated to zero length and overwritten. A pointer to memory for the current ntuple
    /// row ntuple_data must be supplied-this is used to copy ntuples in and out of the file.
    #[doc(alias = "gsl_ntuple_create")]
    pub fn create<P: AsRef<Path>>(filename: P) -> Option<WriteNTuples> {
        let filename = filename.as_ref();
        let filename = filename.to_str().expect("Failed to convert path to str");
        let c_str = CString::new(filename.as_bytes()).unwrap();
        let tmp = unsafe {
            sys::gsl_ntuple_create(c_str.as_ptr() as *mut c_char, ::std::ptr::null_mut(), 0)
        };

        if tmp.is_null() {
            None
        } else {
            Some(Self { n: tmp })
        }
    }

    /// This function writes the current ntuple ntuple->ntuple_data of size ntuple->size to the
    /// corresponding file.
    #[doc(alias = "gsl_ntuple_write")]
    pub fn write<T: Sized>(&mut self, data: &T) -> Value {
        Value::from(unsafe {
            (*self.n).ntuple_data = data as *const T as usize as *mut _;
            (*self.n).size = ::std::mem::size_of::<T>() as _;
            sys::gsl_ntuple_write(self.n)
        })
    }

    /// This function is a synonym for NTuples::write.
    #[doc(alias = "gsl_ntuple_bookdata")]
    pub fn bookdata<T: Sized>(&mut self, data: &T) -> Value {
        Value::from(unsafe {
            (*self.n).ntuple_data = data as *const T as usize as *mut _;
            (*self.n).size = ::std::mem::size_of::<T>() as _;
            sys::gsl_ntuple_bookdata(self.n)
        })
    }
}

impl Drop for WriteNTuples {
    #[doc(alias = "gsl_ntuple_close")]
    fn drop(&mut self) {
        unsafe { sys::gsl_ntuple_close(self.n) };
    }
}

pub struct ReadNTuples {
    n: *mut sys::gsl_ntuple,
}

impl ReadNTuples {
    /// This function opens an existing ntuple file filename for reading and returns a pointer to a
    /// corresponding ntuple struct. The ntuples in the file must have size size. A pointer to
    /// memory for the current ntuple row ntuple_data must be supplied—this is used to copy ntuples
    /// in and out of the file.
    #[doc(alias = "gsl_ntuple_open")]
    pub fn open<P: AsRef<Path>>(filename: P) -> Option<ReadNTuples> {
        let filename = filename.as_ref();
        let filename = filename.to_str().expect("Failed to convert path to str");
        let c_str = CString::new(filename.as_bytes()).unwrap();
        let tmp = unsafe {
            sys::gsl_ntuple_open(c_str.as_ptr() as *mut c_char, ::std::ptr::null_mut(), 0)
        };

        if tmp.is_null() {
            None
        } else {
            Some(Self { n: tmp })
        }
    }

    /// This function reads the current row of the ntuple file for ntuple and stores the values in
    /// ntuple->data.
    #[doc(alias = "gsl_ntuple_read")]
    pub fn read<T: Sized>(&mut self) -> (Value, T) {
        let mut data = MaybeUninit::<T>::uninit();

        let ret = unsafe {
            (*self.n).ntuple_data = data.as_mut_ptr() as *mut _;
            (*self.n).size = ::std::mem::size_of::<T>() as _;
            sys::gsl_ntuple_read(self.n)
        };
        (::Value::from(ret), unsafe { data.assume_init() })
    }
}

impl Drop for ReadNTuples {
    #[doc(alias = "gsl_ntuple_close")]
    fn drop(&mut self) {
        unsafe { sys::gsl_ntuple_close(self.n) };
    }
}

macro_rules! impl_project {
    ($name:ident) => {
        impl $name {
            /// This function updates the histogram `h` from the ntuple `ntuple` using the functions
            /// `value_func` and `select_func`. For each ntuple row where the selection function
            /// `select_func` is non-zero the corresponding value of that row is computed using the function
            /// `value_func` and added to the histogram. Those ntuple rows where `select_func` returns
            /// `false` are ignored. New entries are added to the histogram, so subsequent calls can be used
            /// to accumulate further data in the same histogram.
            #[doc(alias = "gsl_ntuple_project")]
            pub fn project<T: Sized, V: Fn(&T) -> f64, S: Fn(&T) -> bool>(
                &self,
                h: &mut ::Histogram,
                value_func: V,
                select_func: S,
            ) -> Value {
                unsafe extern "C" fn value_trampoline<T: Sized, F: Fn(&T) -> f64>(
                    x: *mut c_void,
                    params: *mut c_void,
                ) -> f64 {
                    let f: &F = &*(params as *const F);
                    let x: &T = &*(x as *const T);
                    f(x)
                }
                unsafe extern "C" fn select_trampoline<T: Sized, F: Fn(&T) -> bool>(
                    x: *mut c_void,
                    params: *mut c_void,
                ) -> i32 {
                    let f: &F = &*(params as *const F);
                    let x: &T = &*(x as *const T);
                    if f(x) {
                        1
                    } else {
                        0
                    }
                }

                let f: Box<V> = Box::new(value_func);
                let mut value_function = sys::gsl_ntuple_value_fn {
                    function: unsafe { ::std::mem::transmute(value_trampoline::<T, V> as usize) },
                    params: Box::into_raw(f) as *mut _,
                };
                let f: Box<S> = Box::new(select_func);
                let mut select_function = sys::gsl_ntuple_select_fn {
                    function: unsafe { ::std::mem::transmute(select_trampoline::<T, S> as usize) },
                    params: Box::into_raw(f) as *mut _,
                };
                Value::from(unsafe {
                    sys::gsl_ntuple_project(
                        h.unwrap_unique(),
                        self.n,
                        &mut value_function,
                        &mut select_function,
                    )
                })
            }
        }
    };
}

impl_project!(WriteNTuples);
impl_project!(ReadNTuples);
