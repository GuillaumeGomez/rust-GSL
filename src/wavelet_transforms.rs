//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
##Transform Functions

This sections describes the actual functions performing the discrete wavelet transform. Note that the transforms use periodic boundary 
conditions. If the signal is not periodic in the sample length then spurious coefficients will appear at the beginning and end of each 
level of the transform.
!*/

/// These functions compute in-place forward and inverse discrete wavelet transforms of length n with stride stride on the array data. The
/// length of the transform n is restricted to powers of two. For the transform version of the function the argument dir can be either
/// forward (+1) or backward (-1). A workspace work of length n must be provided.
/// 
/// For the forward transform, the elements of the original array are replaced by the discrete wavelet transform f_i -> w_{j,k} in a
/// packed triangular storage layout, where j is the index of the level j = 0 ... J-1 and k is the index of the coefficient within each
/// level, k = 0 ... (2^j)-1. The total number of levels is J = \log_2(n). The output data has the following form,
/// 
/// (s_{-1,0}, d_{0,0}, d_{1,0}, d_{1,1}, d_{2,0}, ..., 
///   d_{j,k}, ..., d_{J-1,2^{J-1}-1}) 
/// where the first element is the smoothing coefficient s_{-1,0}, followed by the detail coefficients d_{j,k} for each level j. The
/// backward transform inverts these coefficients to obtain the original data.
/// 
/// These functions return a status of ::Value::Success upon successful completion. ::Inval is returned if n is not an integer power of
/// 2 or if insufficient workspace is provided.
pub mod one_dimension {
    use ffi;
    use enums;

    pub fn transform(w: &::Wavelet, data: &mut [f64], stride: u64, n: u64, dir: ::WaveletDirection,
        work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet_transform(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), stride, n, dir,
            ffi::FFI::unwrap(work)) }
    }

    pub fn transform_forward(w: &::Wavelet, data: &mut [f64], stride: u64, n: u64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet_transform_forward(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), stride, n,
            ffi::FFI::unwrap(work)) }
    }

    pub fn transform_inverse(w: &::Wavelet, data: &mut [f64], stride: u64, n: u64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet_transform_inverse(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), stride, n,
            ffi::FFI::unwrap(work)) }
    }
}

/// The library provides functions to perform two-dimensional discrete wavelet transforms on square matrices. The matrix dimensions must
/// be an integer power of two. There are two possible orderings of the rows and columns in the two-dimensional wavelet transform,
/// referred to as the “standard” and “non-standard” forms.
/// 
/// The “standard” transform performs a complete discrete wavelet transform on the rows of the matrix, followed by a separate complete
/// discrete wavelet transform on the columns of the resulting row-transformed matrix. This procedure uses the same ordering as a
/// two-dimensional Fourier transform.
/// 
/// The “non-standard” transform is performed in interleaved passes on the rows and columns of the matrix for each level of the transform.
/// The first level of the transform is applied to the matrix rows, and then to the matrix columns. This procedure is then repeated across
/// the rows and columns of the data for the subsequent levels of the transform, until the full discrete wavelet transform is complete.
/// The non-standard form of the discrete wavelet transform is typically used in image analysis.
pub mod two_dimension {
    use ffi;
    use enums;

    /// These functions compute two-dimensional in-place forward and inverse discrete wavelet transforms in standard form on the array
    /// data stored in row-major form with dimensions size1 and size2 and physical row length tda. The dimensions must be equal
    /// (square matrix) and are restricted to powers of two. For the transform version of the function the argument dir can be either
    /// forward (+1) or backward (-1). A workspace work of the appropriate size must be provided. On exit, the appropriate elements of
    /// the array data are replaced by their two-dimensional wavelet transform.
    /// 
    /// The functions return a status of ::Value::Success upon successful completion. ::Inval is returned if size1 and size2 are not
    /// equal and integer powers of 2, or if insufficient workspace is provided.
    pub fn transform(w: &::Wavelet, data: &mut [f64], tda: u64, size1: u64, size2: u64, dir: ::WaveletDirection,
        work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_transform(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), tda, size1, size2, dir,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute two-dimensional in-place forward and inverse discrete wavelet transforms in standard form on the array
    /// data stored in row-major form with dimensions size1 and size2 and physical row length tda. The dimensions must be equal
    /// (square matrix) and are restricted to powers of two. For the transform version of the function the argument dir can be either
    /// forward (+1) or backward (-1). A workspace work of the appropriate size must be provided. On exit, the appropriate elements of
    /// the array data are replaced by their two-dimensional wavelet transform.
    /// 
    /// The functions return a status of ::Value::Success upon successful completion. ::Inval is returned if size1 and size2 are not
    /// equal and integer powers of 2, or if insufficient workspace is provided.
    pub fn transform_forward(w: &::Wavelet, data: &mut [f64], tda: u64, size1: u64, size2: u64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_transform_forward(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), tda, size1, size2,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute two-dimensional in-place forward and inverse discrete wavelet transforms in standard form on the array
    /// data stored in row-major form with dimensions size1 and size2 and physical row length tda. The dimensions must be equal
    /// (square matrix) and are restricted to powers of two. For the transform version of the function the argument dir can be either
    /// forward (+1) or backward (-1). A workspace work of the appropriate size must be provided. On exit, the appropriate elements of
    /// the array data are replaced by their two-dimensional wavelet transform.
    /// 
    /// The functions return a status of ::Value::Success upon successful completion. ::Inval is returned if size1 and size2 are not
    /// equal and integer powers of 2, or if insufficient workspace is provided.
    pub fn transform_inverse(w: &::Wavelet, data: &mut [f64], tda: u64, size1: u64, size2: u64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_transform_inverse(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), tda, size1, size2,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the two-dimensional in-place wavelet transform on a matrix a.
    pub fn transform_matrix(w: &::Wavelet, m: &::MatrixF64, dir: ::WaveletDirection, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_transform_matrix(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, ffi::FFI::unwrap(m), dir,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the two-dimensional in-place wavelet transform on a matrix a.
    pub fn transform_matrix_forward(w: &::Wavelet, m: &::MatrixF64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_transform_matrix_forward(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, ffi::FFI::unwrap(m),
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the two-dimensional in-place wavelet transform on a matrix a.
    pub fn transform_matrix_inverse(w: &::Wavelet, m: &::MatrixF64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_transform_matrix_inverse(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, ffi::FFI::unwrap(m),
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the two-dimensional wavelet transform in non-standard form.
    pub fn nstransform(w: &::Wavelet, data: &mut [f64], tda: u64, size1: u64, size2: u64, dir: ::WaveletDirection,
        work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_nstransform(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), tda, size1, size2, dir,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the two-dimensional wavelet transform in non-standard form.
    pub fn nstransform_forward(w: &::Wavelet, data: &mut [f64], tda: u64, size1: u64, size2: u64,
        work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_nstransform_forward(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), tda, size1, size2,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the two-dimensional wavelet transform in non-standard form.
    pub fn nstransform_inverse(w: &::Wavelet, data: &mut [f64], tda: u64, size1: u64, size2: u64,
        work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_nstransform_inverse(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, data.as_mut_ptr(), tda, size1, size2,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the non-standard form of the two-dimensional in-place wavelet transform on a matrix a.
    pub fn nstransform_matrix(w: &::Wavelet, m: &::MatrixF64, dir: ::WaveletDirection, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_nstransform_matrix(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, ffi::FFI::unwrap(m), dir,
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the non-standard form of the two-dimensional in-place wavelet transform on a matrix a.
    pub fn nstransform_matrix_forward(w: &::Wavelet, m: &::MatrixF64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_nstransform_matrix_forward(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, ffi::FFI::unwrap(m),
            ffi::FFI::unwrap(work)) }
    }

    /// These functions compute the non-standard form of the two-dimensional in-place wavelet transform on a matrix a.
    pub fn nstransform_matrix_inverse(w: &::Wavelet, m: &::MatrixF64, work: &::WaveletWorkspace) -> enums::value::Value {
        unsafe { ffi::gsl_wavelet2d_nstransform_matrix_inverse(ffi::FFI::unwrap(w) as *const ffi::gsl_wavelet, ffi::FFI::unwrap(m),
            ffi::FFI::unwrap(work)) }
    }
}