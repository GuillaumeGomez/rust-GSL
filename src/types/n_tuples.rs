//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#N-tuples

This chapter describes functions for creating and manipulating ntuples, sets of values associated with events. The ntuples are stored in 
files. Their values can be extracted in any combination and booked in a histogram using a selection function.

The values to be stored are held in a user-defined data structure, and an ntuple is created associating this data structure with a file. 
The values are then written to the file (normally inside a loop) using the ntuple functions described below.

A histogram can be created from ntuple data by providing a selection function and a value function. The selection function specifies 
whether an event should be included in the subset to be analyzed or not. The value function computes the entry to be added to the histogram 
for each event.

##Histogramming ntuple values

Once an ntuple has been created its contents can be histogrammed in various ways using the function gsl_ntuple_project. Two user-defined 
functions must be provided, a function to select events and a function to compute scalar values. The selection function and the value 
function both accept the ntuple row as a first argument and other parameters as a second argument.

The selection function determines which ntuple rows are selected for histogramming.

##References and Further Reading

Further information on the use of ntuples can be found in the documentation for the CERN packages PAW and HBOOK (available online).
!*/

use ffi;
use enums;
use libc::funcs::c95::stdio::{feof, fread};
use std::marker::PhantomData;
use std::ffi::CString;

pub struct NTuples<T> {
    n: *mut ffi::gsl_ntuple,
    p: PhantomData<T>
}

impl<T> NTuples<T> {
    /// This function creates a new write-only ntuple file filename for ntuples of size size and returns a pointer to the newly created ntuple
    /// struct. Any existing file with the same name is truncated to zero length and overwritten. A pointer to memory for the current ntuple
    /// row ntuple_data must be supplied—this is used to copy ntuples in and out of the file.
    pub fn create(filename: &str, data: &mut T) -> Option<NTuples<T>> {
        let t_data = unsafe { ::std::mem::transmute(data) };
        let c_str = CString::from_slice(filename.as_bytes());
        let tmp = unsafe {
            ffi::gsl_ntuple_create(c_str.as_ptr() as *mut i8, t_data, ::std::mem::size_of::<T>() as u64)
        };

        if tmp.is_null() {
            None
        } else {
            Some(NTuples {
                n: tmp,
                p: PhantomData
            })
        }
    }

    /// This function opens an existing ntuple file filename for reading and returns a pointer to a corresponding ntuple struct. The ntuples
    /// in the file must have size size. A pointer to memory for the current ntuple row ntuple_data must be supplied—this is used to copy
    /// ntuples in and out of the file.
    pub fn open(filename: &str, data: &mut T) -> Option<NTuples<T>> {
        let t_data = unsafe { ::std::mem::transmute(data) };
        let c_str = CString::from_slice(filename.as_bytes());
        let tmp = unsafe {
            ffi::gsl_ntuple_open(c_str.as_ptr() as *mut i8, t_data, ::std::mem::size_of::<T>() as u64)
        };

        if tmp.is_null() {
            None
        } else {
            Some(NTuples {
                n: tmp,
                p: PhantomData
            })
        }
    }

    /// This function writes the current ntuple ntuple->ntuple_data of size ntuple->size to the corresponding file.
    pub fn write(&self) -> enums::value::Value {
        unsafe { ffi::gsl_ntuple_write(self.n) }
    }

    /// This function is a synonym for NTuples::write.
    pub fn bookdata(&self) -> enums::value::Value {
        unsafe { ffi::gsl_ntuple_bookdata(self.n) }
    }

    /// This function reads the current row of the ntuple file for ntuple and stores the values in ntuple->data.
    pub fn read(&self) -> enums::value::Value {
        unsafe { ffi::gsl_ntuple_read(self.n) }
    }

    pub fn project<U, V>(&self, h: &::Histogram, value_func: ::value_function<T, U>, value_arg: &mut U,
        select_func: ::select_function<T, V>, select_arg: &mut V) -> enums::value::Value {
        unsafe {
            loop {
                let nread = fread((*self.n).ntuple_data, (*self.n).size, 1, (*self.n).file);

                if nread == 0 && feof((*self.n).file) != 0 {
                    break;
                }
              
                if nread != 1 {
                    rgsl_error!("failed to read ntuple for projection", ::Value::Failed);
                }

                if select_func(::std::mem::transmute((*self.n).ntuple_data), select_arg) {
                    ffi::gsl_histogram_increment(ffi::FFI::unwrap(h), value_func(::std::mem::transmute((*self.n).ntuple_data), value_arg));
                }
            }

            ::Value::Success
        }
    }
}

#[unsafe_destructor]
impl<T> Drop for NTuples<T> {
    fn drop(&mut self) {
        unsafe { ffi::gsl_ntuple_close(self.n) };
        self.n = ::std::ptr::null_mut();
    }
}

impl<T> ffi::FFI<ffi::gsl_ntuple> for NTuples<T> {
    fn wrap(n: *mut ffi::gsl_ntuple) -> NTuples<T> {
        NTuples {
            n: n,
            p: PhantomData
        }
    }

    fn unwrap(n: &NTuples<T>) -> *mut ffi::gsl_ntuple {
        n.n
    }
}