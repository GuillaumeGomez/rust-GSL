//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Quasi-Random Sequences

This chapter describes functions for generating quasi-random sequences in arbitrary dimensions. A quasi-random sequence progressively 
covers a d-dimensional space with a set of points that are uniformly distributed. Quasi-random sequences are also known as low-discrepancy 
sequences. The quasi-random sequence generators use an interface that is similar to the interface for random number generators, except 
that seeding is not required—each generator produces a single sequence.

##References

The implementations of the quasi-random sequence routines are based on the algorithms described in the following paper,

P. Bratley and B.L. Fox and H. Niederreiter, “Algorithm 738: Programs to Generate Niederreiter’s Low-discrepancy Sequences”, ACM 
Transactions on Mathematical Software, Vol. 20, No. 4, December, 1994, p. 494–495.
!*/

use ffi;
use enums;
use c_vec::CSlice;

pub struct QRng {
    q: *mut ffi::gsl_qrng,
    data: CSlice<i8>
}

impl QRng {
    /// This function returns a pointer to a newly-created instance of a quasi-random sequence generator of type T and dimension d. If
    /// there is insufficient memory to create the generator then the function returns a null pointer and the error handler is invoked
    /// with an error code of ::NoMem.
    pub fn new(t: &QRngType, d: u32) -> Option<QRng> {
        let tmp = unsafe { ffi::gsl_qrng_alloc(t.t, d) };

        if tmp.is_null() {
            None
        } else {
            Some(QRng {
                q: tmp,
                data: unsafe { CSlice::new(tmp as *mut i8, 0) }
            })
        }
    }

    /// This function reinitializes the generator self to its starting point. Note that quasi-random sequences do not use a seed and always
    /// produce the same set of values.
    pub fn init(&self) {
        unsafe { ffi::gsl_qrng_init(self.q) }
    }

    /// This function stores the next point from the sequence generator self in the array x. The space available for x must match the
    /// dimension of the generator. The point x will lie in the range 0 < x_i < 1 for each x_i.
    pub fn get(&self, x: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_qrng_get(self.q, x.as_mut_ptr()) }
    }

    /// This function returns a pointer to the name of the generator.
    pub fn name(&self) -> Option<String> {
        let tmp = unsafe { ffi::gsl_qrng_name(self.q) };

        if tmp.is_null() {
            None
        } else {
            unsafe { Some(String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()) }
        }
    }

    /// These functions return a pointer to the state of generator r and its size.
    pub fn size(&self) -> usize {
        unsafe { ffi::gsl_qrng_size(self.q) }
    }

    /// These functions return a pointer to the state of generator r and its size.
    pub fn state<'r>(&'r mut self) -> &'r mut [i8] {
        let tmp = unsafe { ffi::gsl_qrng_state(self.q) };

        if !tmp.is_null() {
            self.data = unsafe { CSlice::new(tmp as *mut i8, self.size() as usize) };
        }
        self.data.as_mut()
    }

    /// This function copies the quasi-random sequence generator src into the pre-existing generator dest, making dest into an exact copy
    /// of src. The two generators must be of the same type.
    pub fn copy(&self, dest: &QRng) -> enums::Value {
        unsafe { ffi::gsl_qrng_memcpy(dest.q, self.q) }
    }
}

impl Clone for QRng {
    /// This function returns a pointer to a newly created generator which is an exact copy of the generator self.
    fn clone(&self) -> QRng {
        unsafe { ffi::FFI::wrap(ffi::gsl_qrng_clone(self.q)) }
    }
}

impl Drop for QRng {
    fn drop(&mut self) {
        unsafe { ffi::gsl_qrng_free(self.q) };
        self.q = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_qrng> for QRng {
    fn wrap(q: *mut ffi::gsl_qrng) -> QRng {
        QRng {
            q: q,
            data: unsafe { CSlice::new(q as *mut i8, 0) }
        }
    }

    fn unwrap(q: &QRng) -> *mut ffi::gsl_qrng {
        q.q
    }
}

#[derive(Clone, Copy)]
pub struct QRngType {
    t: *const ffi::gsl_qrng_type
}

impl QRngType {
    /// This generator uses the algorithm described in Bratley, Fox, Niederreiter, ACM Trans. Model. Comp. Sim. 2, 195 (1992). It is valid
    /// up to 12 dimensions.
    pub fn niederreiter_2() -> QRngType {
        unsafe {
            QRngType {
                t: ffi::gsl_qrng_niederreiter_2
            }
        }
    }

    /// This generator uses the Sobol sequence described in Antonov, Saleev, USSR Comput. Maths. Math. Phys. 19, 252 (1980). It is valid
    /// up to 40 dimensions.
    pub fn sobol() -> QRngType {
        unsafe {
            QRngType {
                t: ffi::gsl_qrng_sobol
            }
        }
    }

    /// These generators use the Halton and reverse Halton sequences described in J.H. Halton, Numerische Mathematik 2, 84-90 (1960) and
    /// B. Vandewoestyne and R. Cools Computational and Applied Mathematics 189, 1&2, 341-361 (2006). They are valid up to 1229 dimensions.
    pub fn halton() -> QRngType {
        unsafe {
            QRngType {
                t: ffi::gsl_qrng_halton
            }
        }
    }

    pub fn reversehalton() -> QRngType {
        unsafe {
            QRngType {
                t: ffi::gsl_qrng_reversehalton
            }
        }
    }
}

impl ffi::FFI<ffi::gsl_qrng_type> for QRngType {
    fn wrap(t: *mut ffi::gsl_qrng_type) -> QRngType {
        QRngType {
            t: t
        }
    }

    fn unwrap(t: &QRngType) -> *mut ffi::gsl_qrng_type {
        t.t as *mut ffi::gsl_qrng_type
    }
}