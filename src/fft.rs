//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Fast Fourier Transforms (FFTs)

This chapter describes functions for performing Fast Fourier Transforms (FFTs). The library includes radix-2 routines (for lengths which are
a power of two) and mixed-radix routines (which work for any length). For efficiency there are separate versions of the routines for real data
and for complex data. The mixed-radix routines are a reimplementation of the FFTPACK library of Paul Swarztrauber. Fortran code for FFTPACK
is available on Netlib (FFTPACK also includes some routines for sine and cosine transforms but these are currently not available in GSL). For
details and derivations of the underlying algorithms consult the document GSL FFT Algorithms (see FFT References and Further Reading)

## Mathematical Definitions

Fast Fourier Transforms are efficient algorithms for calculating the discrete Fourier transform (DFT),

x_j = \sum_{k=0}^{n-1} z_k \exp(-2\pi i j k / n)

The DFT usually arises as an approximation to the continuous Fourier transform when functions are sampled at discrete intervals in space or time.
The naive evaluation of the discrete Fourier transform is a matrix-vector multiplication W\vec{z}. A general matrix-vector multiplication takes
O(n^2) operations for n data-points. Fast Fourier transform algorithms use a divide-and-conquer strategy to factorize the matrix W into smaller
sub-matrices, corresponding to the integer factors of the length n. If n can be factorized into a product of integers f_1 f_2 ... f_m then the
DFT can be computed in O(n \sum f_i) operations. For a radix-2 FFT this gives an operation count of O(n \log_2 n).

All the FFT functions offer three types of transform: forwards, inverse and backwards, based on the same mathematical definitions. The definition
of the forward Fourier transform, x = FFT(z), is,

x_j = \sum_{k=0}^{n-1} z_k \exp(-2\pi i j k / n)

and the definition of the inverse Fourier transform, x = IFFT(z), is,

z_j = {1 \over n} \sum_{k=0}^{n-1} x_k \exp(2\pi i j k / n).
The factor of 1/n makes this a true inverse. For example, a call to gsl_fft_complex_forward followed by a call to gsl_fft_complex_inverse should
return the original data (within numerical errors).

In general there are two possible choices for the sign of the exponential in the transform/ inverse-transform pair. GSL follows the same convention
as FFTPACK, using a negative exponential for the forward transform. The advantage of this convention is that the inverse transform recreates the
original function with simple Fourier synthesis. Numerical Recipes uses the opposite convention, a positive exponential in the forward transform.

The backwards FFT is simply our terminology for an unscaled version of the inverse FFT,

z^{backwards}_j = \sum_{k=0}^{n-1} x_k \exp(2\pi i j k / n).

When the overall scale of the result is unimportant it is often convenient to use the backwards FFT instead of the inverse to save unnecessary
divisions.

## Overview of complex data FFTs

The inputs and outputs for the complex FFT routines are packed arrays of floating point numbers. In a packed array the real and imaginary parts
of each complex number are placed in alternate neighboring elements. For example, the following definition of a packed array of length 6,

```C
double x[3*2];
gsl_complex_packed_array data = x;
// can be used to hold an array of three complex numbers, z[3], in the following way,

data[0] = Re(z[0])
data[1] = Im(z[0])
data[2] = Re(z[1])
data[3] = Im(z[1])
data[4] = Re(z[2])
data[5] = Im(z[2])
```

The array indices for the data have the same ordering as those in the definition of the DFT—i.e. there are no index transformations or
permutations of the data.

A stride parameter allows the user to perform transforms on the elements `z[stride*i]` instead of `z[i]`. A stride greater than 1 can be used
to take an in-place FFT of the column of a matrix. A stride of 1 accesses the array without any additional spacing between elements.

To perform an FFT on a vector argument, such as gsl_vector_complex * v, use the following definitions (or their equivalents) when calling
the functions described in this chapter:

gsl_complex_packed_array data = v->data;
size_t stride = v->stride;
size_t n = v->size;
For physical applications it is important to remember that the index appearing in the DFT does not correspond directly to a physical frequency.
If the time-step of the DFT is \Delta then the frequency-domain includes both positive and negative frequencies, ranging from -1/(2\Delta)
through 0 to +1/(2\Delta). The positive frequencies are stored from the beginning of the array up to the middle, and the negative frequencies
are stored backwards from the end of the array.

Here is a table which shows the layout of the array data, and the correspondence between the time-domain data z, and the frequency-domain
data x.

index    z               x = FFT(z)

0        z(t = 0)        x(f = 0)
1        z(t = 1)        x(f = 1/(n Delta))
2        z(t = 2)        x(f = 2/(n Delta))
.        ........        ..................
n/2      z(t = n/2)      x(f = +1/(2 Delta),
                               -1/(2 Delta))
.        ........        ..................
n-3      z(t = n-3)      x(f = -3/(n Delta))
n-2      z(t = n-2)      x(f = -2/(n Delta))
n-1      z(t = n-1)      x(f = -1/(n Delta))

When n is even the location n/2 contains the most positive and negative frequencies (+1/(2 \Delta), -1/(2 \Delta)) which are equivalent. If
n is odd then general structure of the table above still applies, but n/2 does not appear.

# Radix-2 FFT routines for complex data

The radix-2 algorithms described in this section are simple and compact, although not necessarily the most efficient. They use the Cooley-Tukey
algorithm to compute in-place complex FFTs for lengths which are a power of 2—no additional storage is required. The corresponding self-sorting
mixed-radix routines offer better performance at the expense of requiring additional working space.

## Mixed-radix FFT routines for complex data

This section describes mixed-radix FFT algorithms for complex data. The mixed-radix functions work for FFTs of any length. They are a
reimplementation of Paul Swarztrauber’s Fortran FFTPACK library. The theory is explained in the review article Self-sorting Mixed-radix FFTs
by Clive Temperton. The routines here use the same indexing scheme and basic algorithms as FFTPACK.

The mixed-radix algorithm is based on sub-transform modules—highly optimized small length FFTs which are combined to create larger FFTs. There
are efficient modules for factors of 2, 3, 4, 5, 6 and 7. The modules for the composite factors of 4 and 6 are faster than combining the modules
for 2*2 and 2*3.

For factors which are not implemented as modules there is a fall-back to a general length-n module which uses Singleton’s method for efficiently
computing a DFT. This module is O(n^2), and slower than a dedicated module would be but works for any length n. Of course, lengths which use the
general length-n module will still be factorized as much as possible. For example, a length of 143 will be factorized into 11*13. Large prime
factors are the worst case scenario, e.g. as found in n=2*3*99991, and should be avoided because their O(n^2) scaling will dominate the run-time
(consult the document GSL FFT Algorithms included in the GSL distribution if you encounter this problem).

The mixed-radix initialization function gsl_fft_complex_wavetable_alloc returns the list of factors chosen by the library for a given length n.
It can be used to check how well the length has been factorized, and estimate the run-time. To a first approximation the run-time scales as
n \sum f_i, where the f_i are the factors of n. For programs under user control you may wish to issue a warning that the transform will be slow
when the length is poorly factorized. If you frequently encounter data lengths which cannot be factorized using the existing small-prime modules
consult GSL FFT Algorithms for details on adding support for other factors.

## Overview of real data FFTs

The functions for real data are similar to those for complex data. However, there is an important difference between forward and inverse transforms.
The Fourier transform of a real sequence is not real. It is a complex sequence with a special symmetry:

z_k = z_{n-k}^*

A sequence with this symmetry is called conjugate-complex or half-complex. This different structure requires different storage layouts for the
forward transform (from real to half-complex) and inverse transform (from half-complex back to real). As a consequence the routines are divided
into two sets: functions in gsl_fft_real which operate on real sequences and functions in gsl_fft_halfcomplex which operate on half-complex sequences.

Functions in gsl_fft_real compute the frequency coefficients of a real sequence. The half-complex coefficients c of a real sequence x are given
by Fourier analysis,

c_k = \sum_{j=0}^{n-1} x_j \exp(-2 \pi i j k /n)

Functions in gsl_fft_halfcomplex compute inverse or backwards transforms. They reconstruct real sequences by Fourier synthesis from their half-complex
frequency coefficients, c,

x_j = {1 \over n} \sum_{k=0}^{n-1} c_k \exp(2 \pi i j k /n)

The symmetry of the half-complex sequence implies that only half of the complex numbers in the output need to be stored. The remaining half can be
reconstructed using the half-complex symmetry condition. This works for all lengths, even and odd—when the length is even the middle value where k=n/2
is also real. Thus only n real numbers are required to store the half-complex sequence, and the transform of a real sequence can be stored in the
same size array as the original data.

The precise storage arrangements depend on the algorithm, and are different for radix-2 and mixed-radix routines. The radix-2 function operates
in-place, which constrains the locations where each element can be stored. The restriction forces real and imaginary parts to be stored far apart.
The mixed-radix algorithm does not have this restriction, and it stores the real and imaginary parts of a given term in neighboring locations (which
is desirable for better locality of memory accesses).
!*/

/// These functions compute forward, backward and inverse FFTs of length n with stride stride, on the packed complex array data using an in-place radix-2
/// decimation-in-time algorithm. The length of the transform is restricted to powers of two. For the transform version of the function
/// the sign argument can be either forward (-1) or backward (+1).
///
/// The functions return a value of ::Value::Success if no errors were detected, or Value::Dom if the length n is not a power of two.
pub mod radix2 {
    use crate::Value;

    #[doc(alias = "gsl_fft_complex_radix2_forward")]
    pub fn forward(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe { sys::gsl_fft_complex_radix2_forward(data.as_mut_ptr(), stride, n) })
    }

    #[doc(alias = "gsl_fft_complex_radix2_transform")]
    pub fn transform(data: &mut [f64], stride: usize, n: usize, sign: ::FftDirection) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_complex_radix2_transform(data.as_mut_ptr(), stride, n, sign.into())
        })
    }

    #[doc(alias = "gsl_fft_complex_radix2_backward")]
    pub fn backward(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe { sys::gsl_fft_complex_radix2_backward(data.as_mut_ptr(), stride, n) })
    }

    #[doc(alias = "gsl_fft_complex_radix2_inverse")]
    pub fn inverse(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe { sys::gsl_fft_complex_radix2_inverse(data.as_mut_ptr(), stride, n) })
    }

    /// This is decimation-in-frequency version of the radix-2 FFT function.
    #[doc(alias = "gsl_fft_complex_radix2_dif_forward")]
    pub fn dif_forward(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_complex_radix2_dif_forward(data.as_mut_ptr(), stride, n)
        })
    }

    /// This is decimation-in-frequency version of the radix-2 FFT function.
    #[doc(alias = "gsl_fft_complex_radix2_dif_transform")]
    pub fn dif_transform(data: &mut [f64], stride: usize, n: usize, sign: ::FftDirection) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_complex_radix2_dif_transform(data.as_mut_ptr(), stride, n, sign.into())
        })
    }

    /// This is decimation-in-frequency version of the radix-2 FFT function.
    #[doc(alias = "gsl_fft_complex_radix2_dif_backward")]
    pub fn dif_backward(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_complex_radix2_dif_backward(data.as_mut_ptr(), stride, n)
        })
    }

    /// This is decimation-in-frequency version of the radix-2 FFT function.
    #[doc(alias = "gsl_fft_complex_radix2_dif_inverse")]
    pub fn dif_inverse(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_complex_radix2_dif_inverse(data.as_mut_ptr(), stride, n)
        })
    }
}

/// This section describes radix-2 FFT algorithms for real data. They use the Cooley-Tukey algorithm to compute in-place FFTs for lengths which
/// are a power of 2.
pub mod real_radix2 {
    use crate::Value;

    /// This function computes an in-place radix-2 FFT of length n and stride stride on the real array data. The output is a half-complex sequence,
    /// which is stored in-place. The arrangement of the half-complex terms uses the following scheme: for k < n/2 the real part of the k-th term
    /// is stored in location k, and the corresponding imaginary part is stored in location n-k. Terms with k > n/2 can be reconstructed using the
    /// symmetry z_k = z^*_{n-k}. The terms for k=0 and k=n/2 are both purely real, and count as a special case. Their real parts are stored in
    /// locations 0 and n/2 respectively, while their imaginary parts which are zero are not stored.
    ///
    /// The following table shows the correspondence between the output data and the equivalent results obtained by considering the input data as
    /// a complex sequence with zero imaginary part (assuming stride=1),
    ///
    /// ```text
    /// complex[0].real    =    data[0]
    /// complex[0].imag    =    0
    /// complex[1].real    =    data[1]
    /// complex[1].imag    =    data[n-1]
    /// ...............         ................
    /// complex[k].real    =    data[k]
    /// complex[k].imag    =    data[n-k]
    /// ...............         ................
    /// complex[n/2].real  =    data[n/2]
    /// complex[n/2].imag  =    0
    /// ...............         ................
    /// complex[k'].real   =    data[k]        k' = n - k
    /// complex[k'].imag   =   -data[n-k]
    /// ...............         ................
    /// complex[n-1].real  =    data[1]
    /// complex[n-1].imag  =   -data[n-1]
    /// ```
    ///
    /// Note that the output data can be converted into the full complex sequence using the function gsl_fft_halfcomplex_radix2_unpack described
    /// below.
    #[doc(alias = "gsl_fft_real_radix2_transform")]
    pub fn transform(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe { sys::gsl_fft_real_radix2_transform(data.as_mut_ptr(), stride, n) })
    }

    /// This function computes the inverse or backwards in-place radix-2 FFT of length n and stride stride on the half-complex sequence data
    /// stored according the output scheme used by gsl_fft_real_radix2. The result is a real array stored in natural order.
    #[doc(alias = "gsl_fft_halfcomplex_radix2_inverse")]
    pub fn inverse(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_halfcomplex_radix2_inverse(data.as_mut_ptr(), stride, n)
        })
    }

    /// This function computes the inverse or backwards in-place radix-2 FFT of length n and stride stride on the half-complex sequence data
    /// stored according the output scheme used by gsl_fft_real_radix2. The result is a real array stored in natural order.
    #[doc(alias = "gsl_fft_halfcomplex_radix2_backward")]
    pub fn backward(data: &mut [f64], stride: usize, n: usize) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_halfcomplex_radix2_backward(data.as_mut_ptr(), stride, n)
        })
    }

    /// This function converts halfcomplex_coefficient, an array of half-complex coefficients as returned by gsl_fft_real_radix2_transform,
    /// into an ordinary complex array, complex_coefficient. It fills in the complex array using the symmetry z_k = z_{n-k}^* to reconstruct
    /// the redundant elements. The algorithm for the conversion is,
    ///
    /// ```C
    /// complex_coefficient[0].real
    ///   = halfcomplex_coefficient[0];
    /// complex_coefficient[0].imag
    ///   = 0.0;
    ///
    /// for (i = 1; i < n - i; i++)
    ///   {
    ///     double hc_real
    ///       = halfcomplex_coefficient[i*stride];
    ///     double hc_imag
    ///       = halfcomplex_coefficient[(n-i)*stride];
    ///     complex_coefficient[i*stride].real = hc_real;
    ///     complex_coefficient[i*stride].imag = hc_imag;
    ///     complex_coefficient[(n - i)*stride].real = hc_real;
    ///     complex_coefficient[(n - i)*stride].imag = -hc_imag;
    ///   }
    ///
    /// if (i == n - i)
    ///   {
    ///     complex_coefficient[i*stride].real
    ///       = halfcomplex_coefficient[(n - 1)*stride];
    ///     complex_coefficient[i*stride].imag
    ///       = 0.0;
    ///   }
    /// ```
    #[doc(alias = "gsl_fft_halfcomplex_radix2_unpack")]
    pub fn unpack(
        halfcomplex_coefficient: &mut [f64],
        complex_coefficient: &mut [f64],
        stride: usize,
        n: usize,
    ) -> Value {
        Value::from(unsafe {
            sys::gsl_fft_halfcomplex_radix2_unpack(
                halfcomplex_coefficient.as_mut_ptr(),
                complex_coefficient.as_mut_ptr(),
                stride,
                n,
            )
        })
    }
}
