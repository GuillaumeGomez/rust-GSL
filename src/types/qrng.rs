//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Quasi-Random Sequences

This chapter describes functions for generating quasi-random sequences in arbitrary dimensions. A quasi-random sequence progressively
covers a d-dimensional space with a set of points that are uniformly distributed. Quasi-random sequences are also known as low-discrepancy
sequences. The quasi-random sequence generators use an interface that is similar to the interface for random number generators, except
that seeding is not required—each generator produces a single sequence.

## References

The implementations of the quasi-random sequence routines are based on the algorithms described in the following paper,

P. Bratley and B.L. Fox and H. Niederreiter, “Algorithm 738: Programs to Generate Niederreiter’s Low-discrepancy Sequences”, ACM
Transactions on Mathematical Software, Vol. 20, No. 4, December, 1994, p. 494–495.
!*/

use crate::Value;
use ffi::FFI;

ffi_wrapper!(QRng, *mut sys::gsl_qrng, gsl_qrng_free);

impl QRng {
    /// This function returns a pointer to a newly-created instance of a quasi-random sequence
    /// generator of type T and dimension d. If there is insufficient memory to create the generator
    /// then the function returns a null pointer and the error handler is invoked with an error code
    /// of [`Value::NoMemory`].
    #[doc(alias = "gsl_qrng_alloc")]
    pub fn new(t: QRngType, d: u32) -> Option<Self> {
        let tmp = unsafe { sys::gsl_qrng_alloc(t.unwrap_shared(), d) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function reinitializes the generator self to its starting point. Note that quasi-random
    /// sequences do not use a seed and always produce the same set of values.
    #[doc(alias = "gsl_qrng_init")]
    pub fn init(&mut self) {
        unsafe { sys::gsl_qrng_init(self.unwrap_unique()) }
    }

    /// This function stores the next point from the sequence generator self in the array x. The
    /// space available for x must match the dimension of the generator. The point x will lie in the
    /// range 0 < x_i < 1 for each x_i.
    #[doc(alias = "gsl_qrng_get")]
    pub fn get(&self, x: &mut [f64]) -> Value {
        Value::from(unsafe { sys::gsl_qrng_get(self.unwrap_shared(), x.as_mut_ptr()) })
    }

    /// This function returns a pointer to the name of the generator.
    #[doc(alias = "gsl_qrng_name")]
    pub fn name(&self) -> Option<String> {
        let tmp = unsafe { sys::gsl_qrng_name(self.unwrap_shared()) };

        if tmp.is_null() {
            None
        } else {
            unsafe {
                Some(
                    String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string(),
                )
            }
        }
    }

    /// These functions return a pointer to the state of generator r and its size.
    #[doc(alias = "gsl_qrng_size")]
    pub fn size(&self) -> usize {
        unsafe { sys::gsl_qrng_size(self.unwrap_shared()) }
    }

    /// This function returns a pointer to the state of generator `self`.
    #[doc(alias = "gsl_qrng_state")]
    pub fn state(&mut self) -> Option<&[i8]> {
        let tmp = unsafe { sys::gsl_qrng_state(self.unwrap_shared()) };

        if tmp.is_null() {
            None
        } else {
            Some(unsafe { ::std::slice::from_raw_parts(tmp as _, self.size()) })
        }
    }

    /// This function returns a pointer to the state of generator `self`.
    #[doc(alias = "gsl_qrng_state")]
    pub fn state_mut(&mut self) -> Option<&mut [i8]> {
        let tmp = unsafe { sys::gsl_qrng_state(self.unwrap_shared()) };

        if tmp.is_null() {
            None
        } else {
            Some(unsafe { ::std::slice::from_raw_parts_mut(tmp as _, self.size()) })
        }
    }

    /// This function copies the quasi-random sequence generator src into the pre-existing generator
    /// `dest`, making dest into an exact copy of `self`. The two generators must be of the same
    /// type.
    #[doc(alias = "gsl_qrng_memcpy")]
    pub fn copy(&self, dest: &mut QRng) -> Value {
        Value::from(unsafe { sys::gsl_qrng_memcpy(dest.unwrap_unique(), self.unwrap_shared()) })
    }
}

impl Clone for QRng {
    /// This function returns a pointer to a newly created generator which is an exact copy of the
    /// generator `self`.
    #[doc(alias = "gsl_qrng_clone")]
    fn clone(&self) -> Self {
        unsafe { Self::wrap(sys::gsl_qrng_clone(self.unwrap_shared())) }
    }
}

ffi_wrapper!(QRngType, *const sys::gsl_qrng_type);

impl QRngType {
    /// This generator uses the algorithm described in Bratley, Fox, Niederreiter, ACM Trans. Model.
    /// Comp. Sim. 2, 195 (1992). It is valid up to 12 dimensions.
    pub fn niederreiter_2() -> QRngType {
        ffi_wrap!(gsl_qrng_niederreiter_2)
    }

    /// This generator uses the Sobol sequence described in Antonov, Saleev, USSR Comput. Maths.
    /// Math. Phys. 19, 252 (1980). It is valid up to 40 dimensions.
    pub fn sobol() -> QRngType {
        ffi_wrap!(gsl_qrng_sobol)
    }

    /// These generators use the Halton and reverse Halton sequences described in J.H. Halton,
    /// Numerische Mathematik 2, 84-90 (1960) and B. Vandewoestyne and R. Cools Computational and
    /// Applied Mathematics 189, 1&2, 341-361 (2006). They are valid up to 1229 dimensions.
    pub fn halton() -> QRngType {
        ffi_wrap!(gsl_qrng_halton)
    }

    pub fn reversehalton() -> QRngType {
        ffi_wrap!(gsl_qrng_reversehalton)
    }
}
