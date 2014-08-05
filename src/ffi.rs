//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use libc::{c_double, c_int, c_uint, c_float, c_void};
use types;
use enums;

pub type CBLAS_INDEX = c_uint;
pub type CBLAS_INDEX_t = CBLAS_INDEX;

extern "C" {
    // Airy functions
    pub fn gsl_sf_airy_Ai(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Ai_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    // Derivatives of Airy Functions
    pub fn gsl_sf_airy_Ai_deriv(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi_deriv(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Ai_deriv_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi_deriv_scaled(x: c_double, mode: types::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_scaled_e(x: c_double, mode: types::gsl_mode_t, result: *mut gsl_sf_result) -> c_int;
    //  Zeros of Airy Functions
    pub fn gsl_sf_airy_zero_Ai(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Ai_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_zero_Bi(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Bi_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    // Zeros of Derivatives of Airy Functions
    pub fn gsl_sf_airy_zero_Ai_deriv(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Ai_deriv_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_zero_Bi_deriv(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Bi_deriv_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;

    // Bessel functions
    // Regular Modified Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_I0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_I1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_In_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_I0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_I1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_In_scaled_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In_scaled_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Regular Modified Spherical Bessel Functions
    // The regular modified spherical Bessel functions i_l(x) are related to the modified Bessel functions of fractional order, i_l(x) = \sqrt{\pi/(2x)} I_{l+1/2}(x)
    pub fn gsl_sf_bessel_i0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_i1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_i2_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i2_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_il_scaled(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_il_scaled_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_il_scaled_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_Inu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Inu_scaled(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_scaled_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Regular Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_J0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_J1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Jn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Jn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Regular Spherical Bessel Functions
    pub fn gsl_sf_bessel_j0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_j1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_j2(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_jl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_jl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_jl_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_jl_steed_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Regular Bessel Function—Fractional Order
    pub fn gsl_sf_bessel_Jnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jnu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_sequence_Jnu_e(nu: c_double, mode: types::gsl_mode_t, size: i64, v: *mut c_double) -> c_int;
    // Irregular Modified Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_K0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_K1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    pub fn gsl_sf_bessel_K0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_K1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_scaled_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_scaled_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Irregular Modified Spherical Bessel Functions
    pub fn gsl_sf_bessel_k0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_k1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_k2_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k2_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_kl_scaled(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_kl_scaled_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_kl_scaled_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Irregular Modified Bessel Functions—Fractional Order
    pub fn gsl_sf_bessel_Knu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Knu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_InKnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_InKnu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Knu_scaled(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Knu_scaled_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Irregular Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_Y0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Y0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Y1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Y1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Yn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Yn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Yn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Irregular Spherical Bessel Functions
    pub fn gsl_sf_bessel_y0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_y1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_y2(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_yl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_yl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_yl_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> c_int;
    // Irregular Bessel Functions—Fractional Order
    pub fn gsl_sf_bessel_Ynu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Ynu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Zeros of Regular Bessel Functions
    pub fn gsl_sf_bessel_zero_J0(s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_J0_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_zero_J1(s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_J1_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_zero_Jnu(nu: f64, s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_Jnu_e(nu: f64, s: c_uint, result: *mut gsl_sf_result) -> c_int;

    // Conical Functions
    pub fn gsl_sf_conicalP_half(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_half_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_conicalP_mhalf(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_mhalf_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_conicalP_0(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_0_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_conicalP_1(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_1_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_conicalP_sph_reg(l: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_sph_reg_e(l: c_int, lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_conicalP_cyl_reg(m: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_cyl_reg_e(m: c_int, lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Level 1 CBLAS functions
    pub fn cblas_sdsdot(N: c_int, alpha: c_float, x: *const c_float, incx: c_int, y: *const c_float, incy: c_int) -> c_float;
    pub fn cblas_dsdot(N: c_int, x: *const c_float, incx: c_int, y: *const c_float, incy: c_int) -> c_double;
    pub fn cblas_sdot(N: c_int, x: *const c_float, incx: c_int, y: *const c_float, incy: c_int) -> c_float;
    pub fn cblas_ddot(N: c_int, x: *const c_float, incx: c_int, y: *const c_float, incy: c_int) -> c_double;
    pub fn cblas_cdotu_sub(N: c_int, x: *const c_void, incx: c_int, y: *const c_void, incy: c_int, dotu: *mut c_void);
    pub fn cblas_cdotc_sub(N: c_int, x: *const c_void, incx: c_int, y: *const c_void, incy: c_int, dotc: *mut c_void);
    pub fn cblas_zdotu_sub(N: c_int, x: *const c_void, incx: c_int, y: *const c_void, incy: c_int, dotu: *mut c_void);
    pub fn cblas_zdotc_sub(N: c_int, x: *const c_void, incx: c_int, y: *const c_void, incy: c_int, dotc: *mut c_void);
    pub fn cblas_snrm2(N: c_int, x: *const c_float, incx: c_int) -> c_float;
    pub fn cblas_sasum(N: c_int, x: *const c_float, incx: c_int) -> c_float;
    pub fn cblas_dnrm2(N: c_int, x: *const c_double, incx: c_int) -> c_double;
    pub fn cblas_dasum(N: c_int, x: *const c_double, incx: c_int) -> c_double;
    pub fn cblas_scnrm2(N: c_int, x: *const c_void, incx: c_int) -> c_float;
    pub fn cblas_scasum(N: c_int, x: *const c_void, incx: c_int) -> c_float;
    pub fn cblas_dznrm2(N: c_int, x: *const c_void, incx: c_int) -> c_double;
    pub fn cblas_dzasum(N: c_int, x: *const c_void, incx: c_int) -> c_double;
    pub fn cblas_isamax(N: c_int, x: *const c_float, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_idamax(N: c_int, x: *const c_double, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_icamax(N: c_int, x: *const c_void, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_izamax(N: c_int, x: *const c_void, incx: c_int) -> CBLAS_INDEX;
    pub fn cblas_sswap(N: c_int, x: *mut c_float, incx: c_int, y: *mut c_float, incy: c_int);
    pub fn cblas_scopy(N: c_int, x: *const c_float, incx: c_int, y: *mut c_float, incy: c_int);
    pub fn cblas_saxpy(N: c_int, alpha: c_float, x: *const c_float, incx: c_int, y: *mut c_float, incy: c_int);
    pub fn cblas_dswap(N: c_int, x: *mut c_double, incx: c_int, y: *mut c_double, incy: c_int);
    pub fn cblas_dcopy(N: c_int, x: *const c_double, incx: c_int, y: *mut c_double, incy: c_int);
    pub fn cblas_daxpy(N: c_int, alpha: c_double, x: *const c_double, incx: c_int, y: *mut c_double, incy: c_int);
    pub fn cblas_cswap(N: c_int, x: *mut c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_ccopy(N: c_int, x: *const c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_caxpy(N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_zswap(N: c_int, x: *mut c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_zcopy(N: c_int, x: *const c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_zaxpy(N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *mut c_void, incy: c_int);
    pub fn cblas_srotg(a: *mut c_float, b: *mut c_float, c: *mut c_float, s: *mut c_float);
    pub fn cblas_srotmg(d1: *mut c_float, d2: *mut c_float, b1: *mut c_float, b2: *const c_float, P: *mut c_float);
    pub fn cblas_srot(N: c_int, x: *mut c_float, incx: c_int, y: *mut c_float, incy: c_int, c: c_float, s: c_float);
    pub fn cblas_srotm(N: c_int, x: *mut c_float, incx: c_int, y: *mut c_float, incy: c_int, P: *const c_float);
    pub fn cblas_drotg(a: *mut c_double, b: *mut c_double, c: *mut c_double, s: *mut c_double);
    pub fn cblas_drotmg(d1: *mut c_double, d2: *mut c_double, b1: *mut c_double, b2: *const c_double, P: *mut c_double);
    pub fn cblas_drot(N: c_int, x: *mut c_double, incx: c_int, y: *mut c_double, incy: c_int, c: c_double, s: c_double);
    pub fn cblas_drotm(N: c_int, x: *mut c_double, incx: c_int, y: *mut c_double, incy: c_int, P: *const c_double);
    pub fn cblas_sscal(N: c_int, alpha: c_float, x: *mut c_float, incx: c_int);
    pub fn cblas_dscal(N: c_int, alpha: c_double, x: *mut c_double, incx: c_int);
    pub fn cblas_cscal(N: c_int, alpha: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_zscal(N: c_int, alpha: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_csscal(N: c_int, alpha: c_float, x: *mut c_void, incx: c_int);
    pub fn cblas_zdscal(N: c_int, alpha: c_double, x: *mut c_void, incx: c_int);
    // Level 2 CBLAS functions
    pub fn cblas_sgemv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, alpha: c_float,
        A: *const c_float, lda: c_int, x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_sgbmv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_strmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stpmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_float, x: *mut c_float, incx: c_int);
    pub fn cblas_strsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stbsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stpsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_float, x: *mut c_float, incx: c_int);
    pub fn cblas_dgemv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, alpha: c_double,
        A: *const c_double, lda: c_int, x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dgbmv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dtrmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtpmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_double, x: *mut c_double, incx: c_int);
    pub fn cblas_dtrsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtbsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtpsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_double, x: *mut c_double, incx: c_int);
    pub fn cblas_cgemv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, alpha: *const c_void,
        A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_cgbmv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_ctrmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctpmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_ctrsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctbsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctpsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_zgemv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, alpha: *const c_void,
        A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zgbmv(order: enums::CblasOrder, transA: enums::CblasTranspose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_ztrmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztpmv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_ztrsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztbsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztpsv(order: enums::CblasOrder, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_ssymv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_float, A: *const c_float, lda: c_int,
        x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_ssbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, K: c_int, alpha: c_float, A: *const c_float,
        lda: c_int, x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_sspmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_float, Ap: *const c_float,
        x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_sger(order: enums::CblasOrder, M: c_int, N: c_int, alpha: c_float, x: *const c_float, incx: c_int, y: *const c_float,
        incy: c_int, A: *mut c_float, lda: c_int);
    pub fn cblas_ssyr(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        A: *mut c_float, lda: c_int);
    pub fn cblas_sspr(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        Ap: *mut c_float);
    pub fn cblas_ssyr2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        y: *const c_float, incy: c_int, A: *mut c_float, lda: c_int);
    pub fn cblas_sspr2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        y: *const c_float, incy: c_int, A: *mut c_float);
    pub fn cblas_dsymv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_double, A: *const c_double, lda: c_int,
        x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dsbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, K: c_int, alpha: c_double, A: *const c_double,
        lda: c_int, x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dspmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_double, Ap: *const c_double,
        x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dger(order: enums::CblasOrder, M: c_int, N: c_int, alpha: c_double, x: *const c_double, incx: c_int, y: *const c_double,
        incy: c_int, A: *mut c_double, lda: c_int);
    pub fn cblas_dsyr(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        A: *mut c_double, lda: c_int);
    pub fn cblas_dspr(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        Ap: *mut c_double);
    pub fn cblas_dsyr2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        y: *const c_double, incy: c_int, A: *mut c_double, lda: c_int);
    pub fn cblas_dspr2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        y: *const c_double, incy: c_int, A: *mut c_double);
    pub fn cblas_csymv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_csbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, K: c_int, alpha: *const c_void, A: *const c_void,
        lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_cspmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, Ap: *const c_void,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_cgeru(order: enums::CblasOrder, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_cgerc(order: enums::CblasOrder, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_cher(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        A: *mut c_void, lda: c_int);
    pub fn cblas_chpr(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        Ap: *mut c_void);
    pub fn cblas_cher2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_chpr2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, Ap: *mut c_void);
    pub fn cblas_zsymv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zsbmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, K: c_int, alpha: *const c_void, A: *const c_void,
        lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zspmv(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, Ap: *const c_void,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zgeru(order: enums::CblasOrder, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_zgerc(order: enums::CblasOrder, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_zher(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        A: *mut c_void, lda: c_int);
    pub fn cblas_zhpr(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        Ap: *mut c_void);
    pub fn cblas_zher2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_zhpr2(order: enums::CblasOrder, uplo: enums::CblasUplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, Ap: *mut c_void);
    // Level 3 CBLAS functions
    pub fn cblas_sgemm(order: enums::CblasOrder, transA: enums::CblasTranspose, transB: enums::CblasTranspose, M: c_int, N: c_int,
        K: c_int, alpha: c_float, A: *const c_float, lda: c_int, B: *const c_float, ldb: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_ssymm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, M: c_int, N: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, B: *const c_float, ldb: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_ssyrk(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_ssyr2k(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, B: *const c_float, ldb: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_strmm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: c_float, A: *const c_float, lda: c_int, B: *mut c_float, ldb: c_int);
    pub fn cblas_strsm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: c_float, A: *const c_float, lda: c_int, B: *mut c_float, ldb: c_int);
    pub fn cblas_dgemm(order: enums::CblasOrder, transA: enums::CblasTranspose, transB: enums::CblasTranspose, M: c_int, N: c_int,
        K: c_int, alpha: c_double, A: *const c_double, lda: c_int, B: *const c_double, ldb: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dsymm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, M: c_int, N: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, B: *const c_double, ldb: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dsyrk(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dsyr2k(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, B: *const c_double, ldb: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dtrmm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: c_double, A: *const c_double, lda: c_int, B: *mut c_double, ldb: c_int);
    pub fn cblas_dtrsm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: c_double, A: *const c_double, lda: c_int, B: *mut c_double, ldb: c_int);
    pub fn cblas_cgemm(order: enums::CblasOrder, transA: enums::CblasTranspose, transB: enums::CblasTranspose, M: c_int, N: c_int,
        K: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void,
        ldc: c_int);
    pub fn cblas_csymm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_csyrk(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_csyr2k(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_ctrmm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_ctrsm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_zgemm(order: enums::CblasOrder, transA: enums::CblasTranspose, transB: enums::CblasTranspose, M: c_int, N: c_int,
        K: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void,
        ldc: c_int);
    pub fn cblas_zsymm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zsyrk(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zsyr2k(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_ztrmm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_ztrsm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, transA: enums::CblasTranspose,
        diag: enums::CblasDiag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_chemm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_cherk(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_cher2k(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zhemm(order: enums::CblasOrder, side: enums::CblasSide, uplo: enums::CblasUplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zherk(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: c_double, A: *const c_void, lda: c_int, beta: c_double, C: *mut c_void, ldc: c_int);
    pub fn cblas_zher2k(order: enums::CblasOrder, uplo: enums::CblasUplo, trans: enums::CblasTranspose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    //to bind later
    //pub fn cblas_xerbla(p: c_int, rout: *const c_char, form: *const c_char, ...);

    // Elementary functions
    pub fn gsl_log1p(x: c_double) -> c_double;
    pub fn gsl_expm1(x: c_double) -> c_double;
    pub fn gsl_hypot(x: c_double, y: c_double) -> c_double;
    pub fn gsl_hypot3(x: c_double, y: c_double, z: c_double) -> c_double;
    pub fn gsl_acosh(x: c_double) -> c_double;
    pub fn gsl_asinh(x: c_double) -> c_double;
    pub fn gsl_atanh(x: c_double) -> c_double;
    pub fn gsl_ldexp(x: c_double, e: c_int) -> c_double;
    pub fn gsl_frexp(x: c_double, e: *mut c_int) -> c_double;

    // Vector functions
    pub fn gsl_vector_alloc(size: u32) -> *mut gsl_vector;
    pub fn gsl_vector_calloc(size: u32) -> *mut gsl_vector;
    pub fn gsl_vector_free(vector: *mut gsl_vector);
    pub fn gsl_vector_get(vector: *mut gsl_vector, i: c_uint) -> c_double;
    pub fn gsl_vector_set(vector: *mut gsl_vector, i: c_uint, x: c_double);
    pub fn gsl_vector_set_all(vector: *mut gsl_vector, x: c_double);
    pub fn gsl_vector_set_zero(vector: *mut gsl_vector);
    pub fn gsl_vector_set_basis(vector: *mut gsl_vector, i: c_uint);
    pub fn gsl_vector_memcpy(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_swap(v: *mut gsl_vector, w: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_swap_elements(vector: *mut gsl_vector, i: c_uint, j: c_uint) -> c_int;
    pub fn gsl_vector_reverse(vector: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_add(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_sub(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_mul(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_div(dest: *mut gsl_vector, src: *const gsl_vector) -> c_int;
    pub fn gsl_vector_scale(dest: *mut gsl_vector, x: c_double) -> c_int;
    pub fn gsl_vector_add_constant(dest: *mut gsl_vector, x: c_double) -> c_int;
    pub fn gsl_vector_max(vector: *mut gsl_vector) -> c_double;
    pub fn gsl_vector_min(vector: *mut gsl_vector) -> c_double;
    pub fn gsl_vector_minmax(vector: *mut gsl_vector, min_out: *mut c_double, max_out: *mut c_double);
    pub fn gsl_vector_max_index(vector: *mut gsl_vector) -> c_uint;
    pub fn gsl_vector_min_index(vector: *mut gsl_vector) -> c_uint;
    pub fn gsl_vector_minmax_index(vector: *mut gsl_vector, imin: *mut c_uint, imax: *mut c_uint);
    pub fn gsl_vector_isnull(vector: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_ispos(vector: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_isneg(vector: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_isnonneg(vector: *mut gsl_vector) -> c_int;
    pub fn gsl_vector_equal(u: *const gsl_vector, v: *const gsl_vector) -> c_int;

    // VectorComplex functions
    pub fn gsl_vector_complex_alloc(size: u32) -> *mut gsl_vector_complex;
    pub fn gsl_vector_complex_calloc(size: u32) -> *mut gsl_vector_complex;
    pub fn gsl_vector_complex_free(vector: *mut gsl_vector_complex);
    pub fn gsl_vector_complex_get(vector: *mut gsl_vector_complex, i: c_uint) -> c_double;
    pub fn gsl_vector_complex_set(vector: *mut gsl_vector_complex, i: c_uint, x: c_double);
    pub fn gsl_vector_complex_set_all(vector: *mut gsl_vector_complex, x: c_double);
    pub fn gsl_vector_complex_set_zero(vector: *mut gsl_vector_complex);
    pub fn gsl_vector_complex_set_basis(vector: *mut gsl_vector_complex, i: c_uint);
    pub fn gsl_vector_complex_memcpy(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_swap(v: *mut gsl_vector_complex, w: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_swap_elements(vector: *mut gsl_vector_complex, i: c_uint, j: c_uint) -> c_int;
    pub fn gsl_vector_complex_reverse(vector: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_add(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_sub(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_mul(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_div(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_scale(dest: *mut gsl_vector_complex, x: c_double) -> c_int;
    pub fn gsl_vector_complex_add_constant(dest: *mut gsl_vector_complex, x: c_double) -> c_int;
    pub fn gsl_vector_complex_max(vector: *mut gsl_vector_complex) -> c_double;
    pub fn gsl_vector_complex_min(vector: *mut gsl_vector_complex) -> c_double;
    pub fn gsl_vector_complex_minmax(vector: *mut gsl_vector_complex, min_out: *mut c_double, max_out: *mut c_double);
    pub fn gsl_vector_complex_max_index(vector: *mut gsl_vector_complex) -> c_uint;
    pub fn gsl_vector_complex_min_index(vector: *mut gsl_vector_complex) -> c_uint;
    pub fn gsl_vector_complex_minmax_index(vector: *mut gsl_vector_complex, imin: *mut c_uint, imax: *mut c_uint);
    pub fn gsl_vector_complex_isnull(vector: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_ispos(vector: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_isneg(vector: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_isnonneg(vector: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_equal(u: *const gsl_vector_complex, v: *const gsl_vector_complex) -> c_int;

    // VectorFloat functions
    pub fn gsl_vector_float_alloc(size: u32) -> *mut gsl_vector_float;
    pub fn gsl_vector_float_calloc(size: u32) -> *mut gsl_vector_float;
    pub fn gsl_vector_float_free(vector: *mut gsl_vector_float);
    pub fn gsl_vector_float_get(vector: *mut gsl_vector_float, i: c_uint) -> c_float;
    pub fn gsl_vector_float_set(vector: *mut gsl_vector_float, i: c_uint, x: c_float);
    pub fn gsl_vector_float_set_all(vector: *mut gsl_vector_float, x: c_float);
    pub fn gsl_vector_float_set_zero(vector: *mut gsl_vector_float);
    pub fn gsl_vector_float_set_basis(vector: *mut gsl_vector_float, i: c_uint);
    pub fn gsl_vector_float_memcpy(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_swap(v: *mut gsl_vector_float, w: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_swap_elements(vector: *mut gsl_vector_float, i: c_uint, j: c_uint) -> c_int;
    pub fn gsl_vector_float_reverse(vector: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_add(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_sub(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_mul(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_div(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_scale(dest: *mut gsl_vector_float, x: c_float) -> c_int;
    pub fn gsl_vector_float_add_constant(dest: *mut gsl_vector_float, x: c_float) -> c_int;
    pub fn gsl_vector_float_max(vector: *mut gsl_vector_float) -> c_float;
    pub fn gsl_vector_float_min(vector: *mut gsl_vector_float) -> c_float;
    pub fn gsl_vector_float_minmax(vector: *mut gsl_vector_float, min_out: *mut c_float, max_out: *mut c_float);
    pub fn gsl_vector_float_max_index(vector: *mut gsl_vector_float) -> c_uint;
    pub fn gsl_vector_float_min_index(vector: *mut gsl_vector_float) -> c_uint;
    pub fn gsl_vector_float_minmax_index(vector: *mut gsl_vector_float, imin: *mut c_uint, imax: *mut c_uint);
    pub fn gsl_vector_float_isnull(vector: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_ispos(vector: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_isneg(vector: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_isnonneg(vector: *mut gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_equal(u: *const gsl_vector_float, v: *const gsl_vector_float) -> c_int;

    // VectorComplexFloat functions
    pub fn gsl_vector_complex_float_alloc(size: u32) -> *mut gsl_vector_complex_float;
    pub fn gsl_vector_complex_float_calloc(size: u32) -> *mut gsl_vector_complex_float;
    pub fn gsl_vector_complex_float_free(vector: *mut gsl_vector_complex_float);
    pub fn gsl_vector_complex_float_get(vector: *mut gsl_vector_complex_float, i: c_uint) -> c_float;
    pub fn gsl_vector_complex_float_set(vector: *mut gsl_vector_complex_float, i: c_uint, x: c_float);
    pub fn gsl_vector_complex_float_set_all(vector: *mut gsl_vector_complex_float, x: c_float);
    pub fn gsl_vector_complex_float_set_zero(vector: *mut gsl_vector_complex_float);
    pub fn gsl_vector_complex_float_set_basis(vector: *mut gsl_vector_complex_float, i: c_uint);
    pub fn gsl_vector_complex_float_memcpy(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_swap(v: *mut gsl_vector_complex_float, w: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_swap_elements(vector: *mut gsl_vector_complex_float, i: c_uint, j: c_uint) -> c_int;
    pub fn gsl_vector_complex_float_reverse(vector: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_add(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_sub(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_mul(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_div(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_scale(dest: *mut gsl_vector_complex_float, x: c_float) -> c_int;
    pub fn gsl_vector_complex_float_add_constant(dest: *mut gsl_vector_complex_float, x: c_float) -> c_int;
    pub fn gsl_vector_complex_float_max(vector: *mut gsl_vector_complex_float) -> c_float;
    pub fn gsl_vector_complex_float_min(vector: *mut gsl_vector_complex_float) -> c_float;
    pub fn gsl_vector_complex_float_minmax(vector: *mut gsl_vector_complex_float, min_out: *mut c_float, max_out: *mut c_float);
    pub fn gsl_vector_complex_float_max_index(vector: *mut gsl_vector_complex_float) -> c_uint;
    pub fn gsl_vector_complex_float_min_index(vector: *mut gsl_vector_complex_float) -> c_uint;
    pub fn gsl_vector_complex_float_minmax_index(vector: *mut gsl_vector_complex_float, imin: *mut c_uint, imax: *mut c_uint);
    pub fn gsl_vector_complex_float_isnull(vector: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_ispos(vector: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_isneg(vector: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_isnonneg(vector: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_equal(u: *const gsl_vector_complex_float, v: *const gsl_vector_complex_float) -> c_int;

    // Level 1 BLAS functions
    pub fn gsl_blas_sdsdot(alpha: c_float, x: *const gsl_vector_float, y: *const gsl_vector_float, result: *mut c_float) -> c_int;
    pub fn gsl_blas_sdot(x: *const gsl_vector_float, y: *const gsl_vector_float, result: *mut c_float) -> c_int;
    pub fn gsl_blas_dsdot(x: *const gsl_vector_float, y: *const gsl_vector_float, result: *mut c_double) -> c_int;
    pub fn gsl_blas_ddot(x: *const gsl_vector, y: *const gsl_vector, result: *mut c_double) -> c_int;
    pub fn gsl_blas_cdotu(x: *const gsl_vector_complex_float, y: *const gsl_vector_complex_float, dotu: *mut gsl_complex_float) -> c_int;
    pub fn gsl_blas_zdotu(x: *const gsl_vector_complex, y: *const gsl_vector_complex, dotu: *mut gsl_complex) -> c_int;
    pub fn gsl_blas_cdotc(x: *const gsl_vector_complex_float, y: *const gsl_vector_complex_float, dotc: *mut gsl_complex_float) -> c_int;
    pub fn gsl_blas_zdotc(x: *const gsl_vector_complex, y: *const gsl_vector_complex, dotc: *mut gsl_complex) -> c_int;
    pub fn gsl_blas_snrm2(x: *const gsl_vector_float) -> c_float;
    pub fn gsl_blas_dnrm2(x: *const gsl_vector) -> c_double;
    pub fn gsl_blas_scnrm2(x: *const gsl_vector_complex_float) -> c_float;
    pub fn gsl_blas_dznrm2(x: *const gsl_vector_complex) -> c_double;
    pub fn gsl_blas_sasum(x: *const gsl_vector_float) -> c_float;
    pub fn gsl_blas_dasum(x: *const gsl_vector) -> c_double;
    pub fn gsl_blas_scasum(x: *const gsl_vector_complex_float) -> c_float;
    pub fn gsl_blas_dzasum(x: *const gsl_vector_complex) -> c_double;

    pub fn gsl_blas_isamax(x: *const gsl_vector_float) -> CBLAS_INDEX_t;
    pub fn gsl_blas_idamax(x: *const gsl_vector) -> CBLAS_INDEX_t;
    pub fn gsl_blas_icamax(x: *const gsl_vector_complex_float) -> CBLAS_INDEX_t;
    pub fn gsl_blas_izamax(x: *const gsl_vector_complex) -> CBLAS_INDEX_t;
    pub fn gsl_blas_sswap(x: *mut gsl_vector_float, y: *mut gsl_vector_float) -> c_int;
    pub fn gsl_blas_dswap(x: *mut gsl_vector, y: *mut gsl_vector) -> c_int;
    pub fn gsl_blas_cswap(x: *mut gsl_vector_complex_float, y: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_blas_zswap(x: *mut gsl_vector_complex, y: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_blas_scopy(x: *const gsl_vector_float, y: *mut gsl_vector_float) -> c_int;
    pub fn gsl_blas_dcopy(x: *const gsl_vector, y: *mut gsl_vector) -> c_int;
    pub fn gsl_blas_ccopy(x: *const gsl_vector_complex_float, y: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_blas_zcopy(x: *const gsl_vector_complex, y: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_blas_saxpy(alpha: c_float, x: *const gsl_vector_float, y: *mut gsl_vector_float) -> c_int;
    pub fn gsl_blas_daxpy(alpha: f64, x: *const gsl_vector, y: *mut gsl_vector) -> c_int;
    pub fn gsl_blas_caxpy(alpha: gsl_complex_float, x: *const gsl_vector_complex_float, y: *mut gsl_vector_complex_float) -> c_int;
    pub fn gsl_blas_zaxpy(alpha: gsl_complex, x: *const gsl_vector_complex, y: *mut gsl_vector_complex) -> c_int;
    pub fn gsl_blas_sscal(alpha: c_float, x: *mut gsl_vector_float);
    pub fn gsl_blas_dscal(alpha: c_double, x: *mut gsl_vector);
    pub fn gsl_blas_cscal(alpha: gsl_complex_float, x: *mut gsl_vector_complex_float);
    pub fn gsl_blas_zscal(alpha: gsl_complex, x: *mut gsl_vector_complex);
    pub fn gsl_blas_csscal(alpha: c_float, x: *mut gsl_vector_complex_float);
    pub fn gsl_blas_zdscal(alpha: c_double, x: *mut gsl_vector_complex);
    pub fn gsl_blas_srotg(a: *mut c_float, b: *mut c_float, c: *mut c_float, d: *mut c_float) -> c_int;
    pub fn gsl_blas_drotg(a: *mut c_double, b: *mut c_double, c: *mut c_double, d: *mut c_double) -> c_int;
    pub fn gsl_blas_srot(a: *mut gsl_vector_float, b: *mut gsl_vector_float, c: c_float, d: c_float) -> c_int;
    pub fn gsl_blas_drot(a: *mut gsl_vector, b: *mut gsl_vector, c: c_double, d: c_double) -> c_int;
    pub fn gsl_blas_srotmg(d1: *mut c_float, d2: *mut c_float, b1: *mut c_float, b2: c_float, P: *mut c_float) -> c_int;
    pub fn gsl_blas_drotmg(d1: *mut c_double, d2: *mut c_double, b1: *mut c_double, b2: c_double, P: *mut c_double) -> c_int;
    pub fn gsl_blas_srotm(x: *mut gsl_vector_float, y: *mut gsl_vector_float, P: *mut c_float) -> c_int;
    pub fn gsl_blas_drotm(x: *mut gsl_vector, y: *mut gsl_vector, P: *mut c_double) -> c_int;
}

pub struct gsl_sf_result {
    pub val: c_double,
    pub err: c_double
}

pub struct gsl_vector_float {
    pub size: u32,
    pub stride: u32,
    pub data: *mut c_float,
    pub block: *mut gsl_block_float,
    pub owner: c_int
}

pub struct gsl_block_float {
    pub size: u32,
    pub data: *mut c_float
}

pub struct gsl_vector {
    pub size: u32,
    pub stride: u32,
    pub data: *mut c_double,
    pub block: *mut gsl_block,
    pub owner: c_int
}

pub struct gsl_block {
    pub size: u32,
    pub data: *mut c_double
}

pub struct gsl_vector_complex_float {
    pub size: u32,
    pub stride: u32,
    pub data: *mut c_float,
    pub block: *mut gsl_block_complex_float,
    pub owner: c_int
}

pub struct gsl_block_complex_float {
    pub size: u32,
    pub data: *mut c_float
}

pub struct gsl_vector_complex {
    pub size: u32,
    pub stride: u32,
    pub data: *mut c_double,
    pub block: *mut gsl_block_complex,
    pub owner: c_int
}

pub struct gsl_block_complex {
    pub size: u32,
    pub data: *mut c_double
}

pub struct gsl_complex {
    pub data: [c_double, ..2]
}

pub struct gsl_complex_float {
    pub data: [c_float, ..2]
}