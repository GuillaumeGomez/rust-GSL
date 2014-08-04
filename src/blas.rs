//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod Blas {
    pub mod Level1 {
        use Gsl;

        pub fn sdsdot(alpha: f32, x: &Gsl::VectorFloat, y: &Gsl::VectorFloat, result: &mut f32) -> i32 {
            unsafe { ::ffi::gsl_blas_sdsdot(alpha, x.get_ffi() as *const ::ffi::gsl_vector_float,
                y.get_ffi() as *const ::ffi::gsl_vector_float, result) }
        }

        pub fn sdot(x: &Gsl::VectorFloat, y: &Gsl::VectorFloat, result: &mut f32) -> i32 {
            unsafe { ::ffi::gsl_blas_sdot(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi() as *const ::ffi::gsl_vector_float,
                result) }
        }

        pub fn dsdot(x: &Gsl::VectorFloat, y: &Gsl::VectorFloat, result: &mut f64) -> i32 {
            unsafe { ::ffi::gsl_blas_dsdot(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi() as *const ::ffi::gsl_vector_float,
                result) }
        }

        pub fn ddot(x: &Gsl::Vector, y: &Gsl::Vector, result: &mut f64) -> i32 {
            unsafe { ::ffi::gsl_blas_ddot(x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi() as *const ::ffi::gsl_vector, result) }
        }

        pub fn cdotu(x: &Gsl::VectorComplexFloat, y: &Gsl::VectorComplexFloat, dotu: &mut Gsl::ComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_cdotu(x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
                y.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(dotu)) }
        }

        pub fn zdotu(x: &Gsl::VectorComplex, y: &Gsl::VectorComplex, dotu: &mut Gsl::Complex) -> i32 {
            unsafe { ::ffi::gsl_blas_zdotu(x.get_ffi() as *const ::ffi::gsl_vector_complex,
                y.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(dotu)) }
        }
    }
}