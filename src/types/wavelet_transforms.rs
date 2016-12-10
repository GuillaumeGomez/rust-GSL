//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Wavelet Transforms

This chapter describes functions for performing Discrete Wavelet Transforms (DWTs). The library includes wavelets for real data in both 
one and two dimensions.

##Definitions

The continuous wavelet transform and its inverse are defined by the relations,

w(s,\tau) = \int f(t) * \psi^*_{s,\tau}(t) dt
and,

f(t) = \int \int_{-\infty}^\infty w(s, \tau) * \psi_{s,\tau}(t) d\tau ds
where the basis functions \psi_{s,\tau} are obtained by scaling and translation from a single function, referred to as the mother wavelet.

The discrete version of the wavelet transform acts on equally-spaced samples, with fixed scaling and translation steps (s, \tau). The 
frequency and time axes are sampled dyadically on scales of 2^j through a level parameter j. The resulting family of functions 
{\psi_{j,n}} constitutes an orthonormal basis for square-integrable signals.

The discrete wavelet transform is an O(N) algorithm, and is also referred to as the fast wavelet transform.

##References and Further Reading

The mathematical background to wavelet transforms is covered in the original lectures by Daubechies,

Ingrid Daubechies. Ten Lectures on Wavelets. CBMS-NSF Regional Conference Series in Applied Mathematics (1992), SIAM, ISBN 0898712742.
An easy to read introduction to the subject with an emphasis on the application of the wavelet transform in various branches of science is,

Paul S. Addison. The Illustrated Wavelet Transform Handbook. Institute of Physics Publishing (2002), ISBN 0750306920.
For extensive coverage of signal analysis by wavelets, wavelet packets and local cosine bases see,

S. G. Mallat. A wavelet tour of signal processing (Second edition). Academic Press (1999), ISBN 012466606X.
The concept of multiresolution analysis underlying the wavelet transform is described in,

S. G. Mallat. Multiresolution Approximations and Wavelet Orthonormal Bases of L^2(R). Transactions of the American Mathematical Society, 
315(1), 1989, 69–87.
S. G. Mallat. A Theory for Multiresolution Signal Decomposition—The Wavelet Representation. IEEE Transactions on Pattern Analysis and 
Machine Intelligence, 11, 1989, 674–693.
The coefficients for the individual wavelet families implemented by the library can be found in the following papers,

I. Daubechies. Orthonormal Bases of Compactly Supported Wavelets. Communications on Pure and Applied Mathematics, 41 (1988) 909–996.
A. Cohen, I. Daubechies, and J.-C. Feauveau. Biorthogonal Bases of Compactly Supported Wavelets. Communications on Pure and Applied 
Mathematics, 45 (1992) 485–560.
The PhysioNet archive of physiological datasets can be found online at http://www.physionet.org/ and is described in the following paper,

Goldberger et al. PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals. 
Circulation 101(23):e215-e220 2000.
!*/

use ffi;

/// The Wavelet structure contains the filter coefficients defining the wavelet and any associated offset parameters.
pub struct Wavelet {
    w: *mut ffi::gsl_wavelet
}

impl Wavelet {
    /// This function allocates and initializes a wavelet object of type T. The parameter k selects the specific member of the wavelet
    /// family. A null pointer is returned if insufficient memory is available or if a unsupported member is selected.
    pub fn new(t: &WaveletType, k: usize) -> Option<Wavelet> {
        let tmp = unsafe { ffi::gsl_wavelet_alloc(t.t, k) };

        if tmp.is_null() {
            None
        } else {
            Some(Wavelet {
                w: tmp
            })
        }
    }

    /// This function returns a pointer to the name of the wavelet family for w.
    pub fn name(&self) -> Option<String> {
        let tmp = unsafe { ffi::gsl_wavelet_name(self.w) };

        if tmp.is_null() {
            None
        } else {
            unsafe { Some(String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()) }
        }
    }
}

impl Drop for Wavelet {
    fn drop(&mut self) {
        unsafe { ffi::gsl_wavelet_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_wavelet> for Wavelet {
    fn wrap(w: *mut ffi::gsl_wavelet) -> Wavelet {
        Wavelet {
            w: w
        }
    }

    fn unwrap(w: &Wavelet) -> *mut ffi::gsl_wavelet {
        w.w
    }
}

/// The centered forms of the wavelets align the coefficients of the various sub-bands on edges. Thus the resulting visualization of the
/// coefficients of the wavelet transform in the phase plane is easier to understand.
#[derive(Clone, Copy)]
pub struct WaveletType {
    t: *const ffi::gsl_wavelet_type
}

impl WaveletType {
    /// This is the Daubechies wavelet family of maximum phase with k/2 vanishing moments. The implemented wavelets are k=4, 6, …, 20, with
    /// k even.
    pub fn daubechies() -> WaveletType {
        unsafe {
            WaveletType {
                t: ffi::gsl_wavelet_daubechies
            }
        }
    }

    /// This is the Daubechies wavelet family of maximum phase with k/2 vanishing moments. The implemented wavelets are k=4, 6, …, 20, with
    /// k even.
    pub fn daubechies_centered() -> WaveletType {
        unsafe {
            WaveletType {
                t: ffi::gsl_wavelet_daubechies_centered
            }
        }
    }

    /// This is the Haar wavelet. The only valid choice of k for the Haar wavelet is k=2.
    pub fn haar() -> WaveletType {
        unsafe {
            WaveletType {
                t: ffi::gsl_wavelet_haar
            }
        }
    }

    /// This is the Haar wavelet. The only valid choice of k for the Haar wavelet is k=2.
    pub fn haar_centered() -> WaveletType {
        unsafe {
            WaveletType {
                t: ffi::gsl_wavelet_haar_centered
            }
        }
    }

    /// This is the biorthogonal B-spline wavelet family of order (i,j). The implemented values of k = 100*i + j are 103, 105, 202, 204,
    /// 206, 208, 301, 303, 305 307, 309.
    pub fn bspline() -> WaveletType {
        unsafe {
            WaveletType {
                t: ffi::gsl_wavelet_bspline
            }
        }
    }

    /// This is the biorthogonal B-spline wavelet family of order (i,j). The implemented values of k = 100*i + j are 103, 105, 202, 204,
    /// 206, 208, 301, 303, 305 307, 309.
    pub fn bspline_centered() -> WaveletType {
        unsafe {
            WaveletType {
                t: ffi::gsl_wavelet_bspline_centered
            }
        }
    }
}

/// The WaveletWorkspace structure contains scratch space of the same size as the input data and is used to hold intermediate results
/// during the transform.
pub struct WaveletWorkspace {
    w: *mut ffi::gsl_wavelet_workspace
}

impl WaveletWorkspace {
    /// This function allocates a workspace for the discrete wavelet transform. To perform a one-dimensional transform on n elements, a
    /// workspace of size n must be provided. For two-dimensional transforms of n-by-n matrices it is sufficient to allocate a workspace
    /// of size n, since the transform operates on individual rows and columns. A null pointer is returned if insufficient memory is
    /// available.
    pub fn new(n: usize) -> Option<WaveletWorkspace> {
        let tmp = unsafe { ffi::gsl_wavelet_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(WaveletWorkspace {
                w: tmp
            })
        }
    }
}

impl Drop for WaveletWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_wavelet_workspace_free(self.w) };
        self.w = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_wavelet_workspace> for WaveletWorkspace {
    fn wrap(w: *mut ffi::gsl_wavelet_workspace) -> WaveletWorkspace {
        WaveletWorkspace {
            w: w
        }
    }

    fn unwrap(w: &WaveletWorkspace) -> *mut ffi::gsl_wavelet_workspace {
        w.w
    }
}