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

        pub fn cdotc(x: &Gsl::VectorComplexFloat, y: &Gsl::VectorComplexFloat, dotc: &mut Gsl::ComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_cdotc(x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
                y.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(dotc)) }
        }

        pub fn zdotc(x: &Gsl::VectorComplex, y: &Gsl::VectorComplex, dotc: &mut Gsl::Complex) -> i32 {
            unsafe { ::ffi::gsl_blas_zdotc(x.get_ffi() as *const ::ffi::gsl_vector_complex,
                y.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(dotc)) }
        }
        
        pub fn snrm2(x: &Gsl::VectorFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_snrm2(x.get_ffi() as *const ::ffi::gsl_vector_float) }
        }

        pub fn dnrm2(x: &Gsl::Vector) -> f64 {
            unsafe { ::ffi::gsl_blas_dnrm2(x.get_ffi() as *const ::ffi::gsl_vector) }
        }

        pub fn scnrm2(x: &Gsl::VectorComplexFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_scnrm2(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
        }

        pub fn dznrm2(x: &Gsl::VectorComplex) -> f64 {
            unsafe { ::ffi::gsl_blas_dznrm2(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
        }

        pub fn sasum(x: &Gsl::VectorFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_sasum(x.get_ffi() as *const ::ffi::gsl_vector_float) }
        }

        pub fn dasum(x: &Gsl::Vector) -> f64 {
            unsafe { ::ffi::gsl_blas_dasum(x.get_ffi() as *const ::ffi::gsl_vector) }
        }

        pub fn scasum(x: &Gsl::VectorComplexFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_scasum(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
        }

        pub fn dzasum(x: &Gsl::VectorComplex) -> f64 {
            unsafe { ::ffi::gsl_blas_dzasum(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
        }

        pub fn isamax(x: &Gsl::VectorFloat) -> u32 {
            unsafe { ::ffi::gsl_blas_isamax(x.get_ffi() as *const ::ffi::gsl_vector_float) }
        }

        pub fn idamax(x: &Gsl::Vector) -> u32 {
            unsafe { ::ffi::gsl_blas_idamax(x.get_ffi() as *const ::ffi::gsl_vector) }
        }

        pub fn icamax(x: &Gsl::VectorComplexFloat) -> u32 {
            unsafe { ::ffi::gsl_blas_icamax(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
        }

        pub fn izamax(x: &Gsl::VectorComplex) -> u32 {
            unsafe { ::ffi::gsl_blas_izamax(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
        }

        pub fn sswap(x: &mut Gsl::VectorFloat, y: &mut Gsl::VectorFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_sswap(x.get_ffi(), y.get_ffi()) }
        }

        pub fn dswap(x: &mut Gsl::Vector, y: &mut Gsl::Vector) -> i32 {
            unsafe { ::ffi::gsl_blas_dswap(x.get_ffi(), y.get_ffi()) }
        }

        pub fn cswap(x: &mut Gsl::VectorComplexFloat, y: &mut Gsl::VectorComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_cswap(x.get_ffi(), y.get_ffi()) }
        }

        pub fn zswap(x: &mut Gsl::VectorComplex, y: &mut Gsl::VectorComplex) -> i32 {
            unsafe { ::ffi::gsl_blas_zswap(x.get_ffi(), y.get_ffi()) }
        }

        pub fn scopy(x: &mut Gsl::VectorFloat, y: &mut Gsl::VectorFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_scopy(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi()) }
        }

        pub fn dcopy(x: &mut Gsl::Vector, y: &mut Gsl::Vector) -> i32 {
            unsafe { ::ffi::gsl_blas_dcopy(x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi()) }
        }

        pub fn ccopy(x: &mut Gsl::VectorComplexFloat, y: &mut Gsl::VectorComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_ccopy(x.get_ffi() as *const ::ffi::gsl_vector_complex_float, y.get_ffi()) }
        }

        pub fn zcopy(x: &mut Gsl::VectorComplex, y: &mut Gsl::VectorComplex) -> i32 {
            unsafe { ::ffi::gsl_blas_zcopy(x.get_ffi() as *const ::ffi::gsl_vector_complex, y.get_ffi()) }
        }

        pub fn saxpy(alpha: f32, x: &Gsl::VectorFloat, y: &mut Gsl::VectorFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_saxpy(alpha, x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi()) }
        }

        pub fn daxpy(alpha: f64, x: &Gsl::Vector, y: &mut Gsl::Vector) -> i32 {
            unsafe { ::ffi::gsl_blas_daxpy(alpha, x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi()) }
        }

        pub fn caxpy(alpha: &Gsl::ComplexFloat, x: &Gsl::VectorComplexFloat, y: &mut Gsl::VectorComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_caxpy(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex_float, y.get_ffi()) }
        }

        pub fn zaxpy(alpha: &Gsl::Complex, x: &Gsl::VectorComplex, y: &mut Gsl::VectorComplex) -> i32 {
            unsafe { ::ffi::gsl_blas_zaxpy(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex, y.get_ffi()) }
        }

        pub fn sscal(alpha: f32, x: &mut Gsl::VectorFloat) {
            unsafe { ::ffi::gsl_blas_sscal(alpha, x.get_ffi()) }
        }

        pub fn dscal(alpha: f64, x: &mut Gsl::Vector) {
            unsafe { ::ffi::gsl_blas_dscal(alpha, x.get_ffi()) }
        }

        pub fn cscal(alpha: &Gsl::ComplexFloat, x: &mut Gsl::VectorComplexFloat) {
            unsafe { ::ffi::gsl_blas_cscal(::std::mem::transmute(*alpha), x.get_ffi()) }
        }

        pub fn zscal(alpha: &Gsl::Complex, x: &mut Gsl::VectorComplex) {
            unsafe { ::ffi::gsl_blas_zscal(::std::mem::transmute(*alpha), x.get_ffi()) }
        }

        pub fn csscal(alpha: f32, x: &mut Gsl::VectorComplexFloat) {
            unsafe { ::ffi::gsl_blas_csscal(alpha, x.get_ffi()) }
        }

        pub fn zdscal(alpha: f64, x: &mut Gsl::VectorComplex) {
            unsafe { ::ffi::gsl_blas_zdscal(alpha, x.get_ffi()) }
        }

        pub fn srotg(a: &mut [f32], b: &mut [f32], c: &mut [f32], d: &mut [f32]) -> i32 {
            unsafe { ::ffi::gsl_blas_srotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), d.as_mut_ptr()) }
        }

        pub fn drotg(a: &mut [f64], b: &mut [f64], c: &mut [f64], d: &mut [f64]) -> i32 {
            unsafe { ::ffi::gsl_blas_drotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), d.as_mut_ptr()) }
        }

        pub fn srot(a: &mut Gsl::VectorFloat, b: &mut Gsl::VectorFloat, c: f32, d: f32) -> i32 {
            unsafe { ::ffi::gsl_blas_srot(a.get_ffi(), b.get_ffi(), c, d) }
        }

        pub fn drot(a: &mut Gsl::Vector, b: &mut Gsl::Vector, c: f64, d: f64) -> i32 {
            unsafe { ::ffi::gsl_blas_drot(a.get_ffi(), b.get_ffi(), c, d) }
        }

        pub fn srotmg(d1: &mut [f32], d2: &mut [f32], b1: &mut [f32], b2: f32, P: &mut [f32]) -> i32 {
            unsafe { ::ffi::gsl_blas_srotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2, P.as_mut_ptr()) }
        }

        pub fn drotmg(d1: &mut [f64], d2: &mut [f64], b1: &mut [f64], b2: f64, P: &mut [f64]) -> i32 {
            unsafe { ::ffi::gsl_blas_drotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2, P.as_mut_ptr()) }
        }

        pub fn srotm(x: &mut Gsl::VectorFloat, y: &mut Gsl::VectorFloat, P: &mut [f32]) -> i32 {
            unsafe { ::ffi::gsl_blas_srotm(x. get_ffi(), y.get_ffi(), P.as_mut_ptr()) }
        }

        pub fn drotm(x: &mut Gsl::Vector, y: &mut Gsl::Vector, P: &mut [f64]) -> i32 {
            unsafe { ::ffi::gsl_blas_drotm(x. get_ffi(), y.get_ffi(), P.as_mut_ptr()) }
        }
    }
}