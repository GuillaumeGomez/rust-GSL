//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod Blas {
    pub mod Level1 {
        use Gsl;

        pub fn sdsdot(alpha: f32, x: &Gsl::VectorFloat, y: &Gsl::VectorFloat, result: &mut f32) -> i32 {
            unsafe { ::ffi::gsl_blas_sdsdot(alpha, y.get_ffi() as *const ::ffi::gsl_vector_float,
                x.get_ffi() as *const ::ffi::gsl_vector_float, result) }
        }
    }
}