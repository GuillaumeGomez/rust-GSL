//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Radix-2 FFT routines for complex data

The radix-2 algorithms described in this section are simple and compact, although not necessarily the most efficient. They use the Cooley-Tukey 
algorithm to compute in-place complex FFTs for lengths which are a power of 2—no additional storage is required. The corresponding self-sorting 
mixed-radix routines offer better performance at the expense of requiring additional working space.
!*/

use enums;
use ffi;

/// This function computes forward FFTs of length n with stride stride, on the packed complex array data using an in-place radix-2
/// decimation-in-time algorithm. The length of the transform is restricted to powers of two. For the transform version of the function
/// the sign argument can be either forward (-1) or backward (+1).
/// 
/// The functions return a value of Value::Success if no errors were detected, or Value::Dom if the length n is not a power of two.
pub fn forward(data: &mut[f64], stride: u64, n: u64) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_forward(data.as_mut_ptr(), stride, n) }
}

/// The length of the transform is restricted to powers of two. For the transform version of the function the sign argument can be either
/// forward (-1) or backward (+1).
/// 
/// The functions return a value of Value::Success if no errors were detected, or Value::Dom if the length n is not a power of two.
pub fn transform(data: &mut[f64], stride: u64, n: u64, sign: enums::FftDirection) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_transform(data.as_mut_ptr(), stride, n, sign) }
}

/// This function computes backward FFTs of length n with stride stride, on the packed complex array data using an in-place radix-2
/// decimation-in-time algorithm. The length of the transform is restricted to powers of two. For the transform version of the function
/// the sign argument can be either forward (-1) or backward (+1).
/// 
/// The functions return a value of Value::Success if no errors were detected, or Value::Dom if the length n is not a power of two.
pub fn backward(data: &mut[f64], stride: u64, n: u64) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_backward(data.as_mut_ptr(), stride, n) }
}

/// This function computes inverse FFTs of length n with stride stride, on the packed complex array data using an in-place radix-2
/// decimation-in-time algorithm. The length of the transform is restricted to powers of two. For the transform version of the function
/// the sign argument can be either forward (-1) or backward (+1).
/// 
/// The functions return a value of Value::Success if no errors were detected, or Value::Dom if the length n is not a power of two.
pub fn inverse(data: &mut[f64], stride: u64, n: u64) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_inverse(data.as_mut_ptr(), stride, n) }
}

/// This is decimation-in-frequency version of the radix-2 FFT function.
pub fn diff_forward(data: &mut[f64], stride: u64, n: u64) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_diff_forward(data.as_mut_ptr(), stride, n) }
}

/// This is decimation-in-frequency version of the radix-2 FFT function.
pub fn diff_transform(data: &mut[f64], stride: u64, n: u64, sign: enums::FftDirection) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_diff_transform(data.as_mut_ptr(), stride, n, sign) }
}

/// This is decimation-in-frequency version of the radix-2 FFT function.
pub fn diff_backward(data: &mut[f64], stride: u64, n: u64) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_diff_backward(data.as_mut_ptr(), stride, n) }
}

/// This is decimation-in-frequency version of the radix-2 FFT function.
pub fn diff_inverse(data: &mut[f64], stride: u64, n: u64) -> enums::Value {
    unsafe { ffi::gsl_fft_complex_radix2_diff_inverse(data.as_mut_ptr(), stride, n) }
}