//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod Blas {
    pub mod Level1 {
        use Gsl;

        /// This function computes the sum \alpha + x^T y for the vectors x and y, returning the result in result.
        pub fn sdsdot(alpha: f32, x: &Gsl::VectorFloat, y: &Gsl::VectorFloat, result: &mut f32) -> i32 {
            unsafe { ::ffi::gsl_blas_sdsdot(alpha, x.get_ffi() as *const ::ffi::gsl_vector_float,
                y.get_ffi() as *const ::ffi::gsl_vector_float, result) }
        }

        /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
        pub fn sdot(x: &Gsl::VectorFloat, y: &Gsl::VectorFloat, result: &mut f32) -> i32 {
            unsafe { ::ffi::gsl_blas_sdot(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi() as *const ::ffi::gsl_vector_float,
                result) }
        }

        /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
        pub fn dsdot(x: &Gsl::VectorFloat, y: &Gsl::VectorFloat, result: &mut f64) -> i32 {
            unsafe { ::ffi::gsl_blas_dsdot(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi() as *const ::ffi::gsl_vector_float,
                result) }
        }

        /// This function computes the scalar product x^T y for the vectors x and y, returning the result in result.
        pub fn ddot(x: &Gsl::Vector, y: &Gsl::Vector, result: &mut f64) -> i32 {
            unsafe { ::ffi::gsl_blas_ddot(x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi() as *const ::ffi::gsl_vector, result) }
        }

        /// This function computes the complex scalar product x^T y for the vectors x and y, returning the result in dotu.
        pub fn cdotu(x: &Gsl::VectorComplexFloat, y: &Gsl::VectorComplexFloat, dotu: &mut Gsl::ComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_cdotu(x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
                y.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(dotu)) }
        }

        /// This function computes the complex scalar product x^T y for the vectors x and y, returning the result in dotu.
        pub fn zdotu(x: &Gsl::VectorComplex, y: &Gsl::VectorComplex, dotu: &mut Gsl::Complex) -> i32 {
            unsafe { ::ffi::gsl_blas_zdotu(x.get_ffi() as *const ::ffi::gsl_vector_complex,
                y.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(dotu)) }
        }

        /// This function computes the complex conjugate scalar product x^H y for the vectors x and y, returning the result in dotc.
        pub fn cdotc(x: &Gsl::VectorComplexFloat, y: &Gsl::VectorComplexFloat, dotc: &mut Gsl::ComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_cdotc(x.get_ffi() as *const ::ffi::gsl_vector_complex_float,
                y.get_ffi() as *const ::ffi::gsl_vector_complex_float, ::std::mem::transmute(dotc)) }
        }

        /// This function computes the complex conjugate scalar product x^H y for the vectors x and y, returning the result in dotc.
        pub fn zdotc(x: &Gsl::VectorComplex, y: &Gsl::VectorComplex, dotc: &mut Gsl::Complex) -> i32 {
            unsafe { ::ffi::gsl_blas_zdotc(x.get_ffi() as *const ::ffi::gsl_vector_complex,
                y.get_ffi() as *const ::ffi::gsl_vector_complex, ::std::mem::transmute(dotc)) }
        }

        /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
        pub fn snrm2(x: &Gsl::VectorFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_snrm2(x.get_ffi() as *const ::ffi::gsl_vector_float) }
        }

        /// This function computes the Euclidean norm ||x||_2 = \sqrt {\sum x_i^2} of the vector x.
        pub fn dnrm2(x: &Gsl::Vector) -> f64 {
            unsafe { ::ffi::gsl_blas_dnrm2(x.get_ffi() as *const ::ffi::gsl_vector) }
        }

        /// This function computes the Euclidean norm of the complex vector x,
        /// 
        /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
        pub fn scnrm2(x: &Gsl::VectorComplexFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_scnrm2(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
        }

        /// This function computes the Euclidean norm of the complex vector x,
        /// 
        /// ||x||_2 = \sqrt {\sum (\Re(x_i)^2 + \Im(x_i)^2)}.
        pub fn dznrm2(x: &Gsl::VectorComplex) -> f64 {
            unsafe { ::ffi::gsl_blas_dznrm2(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
        }

        /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
        pub fn sasum(x: &Gsl::VectorFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_sasum(x.get_ffi() as *const ::ffi::gsl_vector_float) }
        }

        /// This function computes the absolute sum \sum |x_i| of the elements of the vector x.
        pub fn dasum(x: &Gsl::Vector) -> f64 {
            unsafe { ::ffi::gsl_blas_dasum(x.get_ffi() as *const ::ffi::gsl_vector) }
        }

        /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
        pub fn scasum(x: &Gsl::VectorComplexFloat) -> f32 {
            unsafe { ::ffi::gsl_blas_scasum(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
        }

        /// This function computes the sum of the magnitudes of the real and imaginary parts of the complex vector x, \sum |\Re(x_i)| + |\Im(x_i)|.
        pub fn dzasum(x: &Gsl::VectorComplex) -> f64 {
            unsafe { ::ffi::gsl_blas_dzasum(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
        }

        /// This function returns the index of the largest element of the vector x.
        /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
        /// If the largest value occurs several times then the index of the first occurrence is returned.
        pub fn isamax(x: &Gsl::VectorFloat) -> u32 {
            unsafe { ::ffi::gsl_blas_isamax(x.get_ffi() as *const ::ffi::gsl_vector_float) }
        }

        /// This function returns the index of the largest element of the vector x.
        /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
        /// If the largest value occurs several times then the index of the first occurrence is returned.
        pub fn idamax(x: &Gsl::Vector) -> u32 {
            unsafe { ::ffi::gsl_blas_idamax(x.get_ffi() as *const ::ffi::gsl_vector) }
        }

        /// This function returns the index of the largest element of the vector x.
        /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
        /// If the largest value occurs several times then the index of the first occurrence is returned.
        pub fn icamax(x: &Gsl::VectorComplexFloat) -> u32 {
            unsafe { ::ffi::gsl_blas_icamax(x.get_ffi() as *const ::ffi::gsl_vector_complex_float) }
        }

        /// This function returns the index of the largest element of the vector x.
        /// The largest element is determined by its absolute magnitude for real vectors and by the sum of the magnitudes of the real and imaginary parts |\Re(x_i)| + |\Im(x_i)| for complex vectors.
        /// If the largest value occurs several times then the index of the first occurrence is returned.
        pub fn izamax(x: &Gsl::VectorComplex) -> u32 {
            unsafe { ::ffi::gsl_blas_izamax(x.get_ffi() as *const ::ffi::gsl_vector_complex) }
        }

        /// This function exchanges the elements of the vectors x and y.
        pub fn sswap(x: &mut Gsl::VectorFloat, y: &mut Gsl::VectorFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_sswap(x.get_ffi(), y.get_ffi()) }
        }

        /// This function exchanges the elements of the vectors x and y.
        pub fn dswap(x: &mut Gsl::Vector, y: &mut Gsl::Vector) -> i32 {
            unsafe { ::ffi::gsl_blas_dswap(x.get_ffi(), y.get_ffi()) }
        }

        /// This function exchanges the elements of the vectors x and y.
        pub fn cswap(x: &mut Gsl::VectorComplexFloat, y: &mut Gsl::VectorComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_cswap(x.get_ffi(), y.get_ffi()) }
        }

        /// This function exchanges the elements of the vectors x and y.
        pub fn zswap(x: &mut Gsl::VectorComplex, y: &mut Gsl::VectorComplex) -> i32 {
            unsafe { ::ffi::gsl_blas_zswap(x.get_ffi(), y.get_ffi()) }
        }

        /// This function copy the elements of the vector x into the vector y.
        pub fn scopy(x: &mut Gsl::VectorFloat, y: &mut Gsl::VectorFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_scopy(x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi()) }
        }

        /// This function copy the elements of the vector x into the vector y.
        pub fn dcopy(x: &mut Gsl::Vector, y: &mut Gsl::Vector) -> i32 {
            unsafe { ::ffi::gsl_blas_dcopy(x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi()) }
        }

        /// This function copy the elements of the vector x into the vector y.
        pub fn ccopy(x: &mut Gsl::VectorComplexFloat, y: &mut Gsl::VectorComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_ccopy(x.get_ffi() as *const ::ffi::gsl_vector_complex_float, y.get_ffi()) }
        }

        /// This function copy the elements of the vector x into the vector y.
        pub fn zcopy(x: &mut Gsl::VectorComplex, y: &mut Gsl::VectorComplex) -> i32 {
            unsafe { ::ffi::gsl_blas_zcopy(x.get_ffi() as *const ::ffi::gsl_vector_complex, y.get_ffi()) }
        }

        /// This function computes the sum y = \alpha x + y for the vectors x and y.
        pub fn saxpy(alpha: f32, x: &Gsl::VectorFloat, y: &mut Gsl::VectorFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_saxpy(alpha, x.get_ffi() as *const ::ffi::gsl_vector_float, y.get_ffi()) }
        }

        /// This function computes the sum y = \alpha x + y for the vectors x and y.
        pub fn daxpy(alpha: f64, x: &Gsl::Vector, y: &mut Gsl::Vector) -> i32 {
            unsafe { ::ffi::gsl_blas_daxpy(alpha, x.get_ffi() as *const ::ffi::gsl_vector, y.get_ffi()) }
        }

        /// This function computes the sum y = \alpha x + y for the vectors x and y.
        pub fn caxpy(alpha: &Gsl::ComplexFloat, x: &Gsl::VectorComplexFloat, y: &mut Gsl::VectorComplexFloat) -> i32 {
            unsafe { ::ffi::gsl_blas_caxpy(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex_float, y.get_ffi()) }
        }

        /// This function computes the sum y = \alpha x + y for the vectors x and y.
        pub fn zaxpy(alpha: &Gsl::Complex, x: &Gsl::VectorComplex, y: &mut Gsl::VectorComplex) -> i32 {
            unsafe { ::ffi::gsl_blas_zaxpy(::std::mem::transmute(*alpha), x.get_ffi() as *const ::ffi::gsl_vector_complex, y.get_ffi()) }
        }

        /// This function rescales the vector x by the multiplicative factor alpha.
        pub fn sscal(alpha: f32, x: &mut Gsl::VectorFloat) {
            unsafe { ::ffi::gsl_blas_sscal(alpha, x.get_ffi()) }
        }

        /// This function rescales the vector x by the multiplicative factor alpha.
        pub fn dscal(alpha: f64, x: &mut Gsl::Vector) {
            unsafe { ::ffi::gsl_blas_dscal(alpha, x.get_ffi()) }
        }

        /// This function rescales the vector x by the multiplicative factor alpha.
        pub fn cscal(alpha: &Gsl::ComplexFloat, x: &mut Gsl::VectorComplexFloat) {
            unsafe { ::ffi::gsl_blas_cscal(::std::mem::transmute(*alpha), x.get_ffi()) }
        }

        /// This function rescales the vector x by the multiplicative factor alpha.
        pub fn zscal(alpha: &Gsl::Complex, x: &mut Gsl::VectorComplex) {
            unsafe { ::ffi::gsl_blas_zscal(::std::mem::transmute(*alpha), x.get_ffi()) }
        }

        /// This function rescales the vector x by the multiplicative factor alpha.
        pub fn csscal(alpha: f32, x: &mut Gsl::VectorComplexFloat) {
            unsafe { ::ffi::gsl_blas_csscal(alpha, x.get_ffi()) }
        }

        /// This function rescales the vector x by the multiplicative factor alpha.
        pub fn zdscal(alpha: f64, x: &mut Gsl::VectorComplex) {
            unsafe { ::ffi::gsl_blas_zdscal(alpha, x.get_ffi()) }
        }

        /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
        /// 
        /// [  c  s ] [ a ] = [ r ]
        /// 
        /// [ -s  c ] [ b ]   [ 0 ]
        /// 
        /// The variables a and b are overwritten by the routine.
        pub fn srotg(a: &mut [f32], b: &mut [f32], c: &mut [f32], d: &mut [f32]) -> i32 {
            unsafe { ::ffi::gsl_blas_srotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), d.as_mut_ptr()) }
        }

        /// This function computes a Givens rotation (c,s) which zeroes the vector (a,b),
        /// 
        /// [  c  s ] [ a ] = [ r ]
        /// 
        /// [ -s  c ] [ b ]   [ 0 ]
        /// 
        /// The variables a and b are overwritten by the routine.
        pub fn drotg(a: &mut [f64], b: &mut [f64], c: &mut [f64], d: &mut [f64]) -> i32 {
            unsafe { ::ffi::gsl_blas_drotg(a.as_mut_ptr(), b.as_mut_ptr(), c.as_mut_ptr(), d.as_mut_ptr()) }
        }

        /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
        pub fn srot(a: &mut Gsl::VectorFloat, b: &mut Gsl::VectorFloat, c: f32, d: f32) -> i32 {
            unsafe { ::ffi::gsl_blas_srot(a.get_ffi(), b.get_ffi(), c, d) }
        }

        /// This function applies a Givens rotation (x', y') = (c x + s y, -s x + c y) to the vectors x, y.
        pub fn drot(a: &mut Gsl::Vector, b: &mut Gsl::Vector, c: f64, d: f64) -> i32 {
            unsafe { ::ffi::gsl_blas_drot(a.get_ffi(), b.get_ffi(), c, d) }
        }

        /// This function computes a modified Givens transformation.
        /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
        pub fn srotmg(d1: &mut [f32], d2: &mut [f32], b1: &mut [f32], b2: f32, P: &mut [f32]) -> i32 {
            unsafe { ::ffi::gsl_blas_srotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2, P.as_mut_ptr()) }
        }

        /// This function computes a modified Givens transformation.
        /// The modified Givens transformation is defined in the original Level-1 BLAS specification, given in the references.
        pub fn drotmg(d1: &mut [f64], d2: &mut [f64], b1: &mut [f64], b2: f64, P: &mut [f64]) -> i32 {
            unsafe { ::ffi::gsl_blas_drotmg(d1.as_mut_ptr(), d2.as_mut_ptr(), b1.as_mut_ptr(), b2, P.as_mut_ptr()) }
        }

        /// This function applies a modified Givens transformation.
        pub fn srotm(x: &mut Gsl::VectorFloat, y: &mut Gsl::VectorFloat, P: &mut [f32]) -> i32 {
            unsafe { ::ffi::gsl_blas_srotm(x. get_ffi(), y.get_ffi(), P.as_mut_ptr()) }
        }

        /// This function applies a modified Givens transformation.
        pub fn drotm(x: &mut Gsl::Vector, y: &mut Gsl::Vector, P: &mut [f64]) -> i32 {
            unsafe { ::ffi::gsl_blas_drotm(x. get_ffi(), y.get_ffi(), P.as_mut_ptr()) }
        }
    }
}