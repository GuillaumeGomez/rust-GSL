//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::ffi::FFI;
use crate::{
    vector::VectorMut,
    Value
};
use paste::paste;

macro_rules! gsl_fft_wavetable {
    ($rust_name:ident, $name:ident, $complex_rust_name:ident, $complex_name:ident, $ty:ident $(, $extra:ident)?) => (
paste! {

ffi_wrapper!(
    $rust_name,
    *mut sys::[<$name _wavetable $($extra)?>],
    [<$name _wavetable $($extra)? _free>]
);

impl $rust_name {
    /// This function prepares a trigonometric lookup table for a complex FFT of length n. The
    /// function returns a pointer to the newly allocated gsl_fft_complex_wavetable if no errors
    /// were detected, and a null pointer in the case of error. The length n is factorized into a
    /// product of subtransforms, and the factors and their trigonometric coefficients are stored in
    /// the wavetable. The trigonometric coefficients are computed using direct calls to sin and
    /// cos, for accuracy. Recursion relations could be used to compute the lookup table faster, but
    /// if an application performs many FFTs of the same length then this computation is a one-off
    /// overhead which does not affect the final throughput.
    ///
    /// The wavetable structure can be used repeatedly for any transform of the same length. The
    /// table is not modified by calls to any of the other FFT functions. The same wavetable can be
    /// used for both forward and backward (or inverse) transforms of a given length.
    #[doc(alias = $name _wavetable $($extra)? _alloc)]
    pub fn new(n: usize) -> Option<Self> {
        let tmp = unsafe { sys::[<$name _wavetable $($extra)? _alloc>](n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    pub fn nf(&self) -> usize {
        unsafe { (*self.unwrap_shared()).nf }
    }

    pub fn factor(&self) -> &[usize] {
        unsafe { &(*self.unwrap_shared()).factor }
    }

    pub fn factor_mut(&mut self) -> &mut [usize] {
        unsafe { &mut (*self.unwrap_unique()).factor }
    }
}


ffi_wrapper!(
    $complex_rust_name,
    *mut sys::$complex_name,
    [<$complex_name _free>]
);

impl $complex_rust_name {
    /// This function allocates a workspace for a complex transform of length n.
    #[doc(alias = $complex_name _alloc)]
    pub fn new(n: usize) -> Option<Self> {
        let tmp = unsafe { sys::[<$complex_name _alloc>](n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    #[doc(alias = $name $($extra)? _forward)]
    pub fn forward<V: VectorMut<$ty> + ?Sized>(
        &mut self,
        data: &mut V,
        wavetable: &$rust_name,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::[<$name $($extra)? _forward>](
                V::as_mut_slice(data).as_mut_ptr(),
                V::stride(data),
                V::len(data),
                wavetable.unwrap_shared(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = $name $($extra)? _transform)]
    pub fn transform<V: VectorMut<$ty> + ?Sized>(
        &mut self,
        data: &mut V,
        wavetable: &$rust_name,
        sign: crate::FftDirection,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::[<$name $($extra)? _transform>](
                V::as_mut_slice(data).as_mut_ptr(),
                V::stride(data),
                V::len(data),
                wavetable.unwrap_shared(),
                self.unwrap_unique(),
                sign.into(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = $name $($extra)? _backward)]
    pub fn backward<V: VectorMut<$ty> + ?Sized>(
        &mut self,
        data: &mut V,
        wavetable: &$rust_name,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::[<$name $($extra)? _backward>](
                V::as_mut_slice(data).as_mut_ptr(),
                V::stride(data),
                V::len(data),
                wavetable.unwrap_shared(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }

    #[doc(alias = $name $($extra)? _inverse)]
    pub fn inverse<V: VectorMut<$ty> + ?Sized>(
        &mut self,
        data: &mut V,
        wavetable: &$rust_name,
    ) -> Result<(), Value> {
        let ret = unsafe {
            sys::[<$name $($extra)? _inverse>](
                V::as_mut_slice(data).as_mut_ptr(),
                V::stride(data),
                V::len(data),
                wavetable.unwrap_shared(),
                self.unwrap_unique(),
            )
        };
        result_handler!(ret, ())
    }
}

} // end of paste! block
); // end of macro block
}

gsl_fft_wavetable!(
    FftComplexF64WaveTable,
    gsl_fft_complex,
    FftComplexF64Workspace,
    gsl_fft_complex_workspace,
    f64
);
gsl_fft_wavetable!(
    FftComplexF32WaveTable,
    gsl_fft_complex,
    FftComplexF32Workspace,
    gsl_fft_complex_workspace_float,
    f32,
    _float
);
