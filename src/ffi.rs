//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![allow(improper_ctypes)]

use libc::{c_double, c_int, c_uint, c_float, c_void, size_t, c_ulong, c_char, FILE};
use enums;
use cblas;

pub type CBLAS_INDEX = c_uint;
pub type CBLAS_INDEX_t = CBLAS_INDEX;
pub type CBLAS_TRANSPOSE_t = cblas::Transpose;
pub type CBLAS_UPLO_t = cblas::Uplo;
pub type CBLAS_DIAG_t = cblas::Diag;
pub type CBLAS_SIDE_t = cblas::Side;
pub type gsl_complex_packed_ptr = *mut c_double;
pub type gsl_complex_packed_array = *mut c_double;
pub type coord = c_int;

pub trait FFI<T> {
    fn wrap(r: *mut T) -> Self;
    fn unwrap(&Self) -> *mut T;
}

extern "C" {
    pub static gsl_rng_mt19937 : *const gsl_rng_type;
    pub static gsl_rng_ranlxs0 : *const gsl_rng_type;
    pub static gsl_rng_ranlxs1 : *const gsl_rng_type;
    pub static gsl_rng_ranlxs2 : *const gsl_rng_type;
    pub static gsl_rng_ranlxd1 : *const gsl_rng_type;
    pub static gsl_rng_ranlxd2 : *const gsl_rng_type;
    pub static gsl_rng_ranlux : *const gsl_rng_type;
    pub static gsl_rng_ranlux389 : *const gsl_rng_type;
    pub static gsl_rng_cmrg : *const gsl_rng_type;
    pub static gsl_rng_mrg : *const gsl_rng_type;
    pub static gsl_rng_taus : *const gsl_rng_type;
    pub static gsl_rng_taus2 : *const gsl_rng_type;
    pub static gsl_rng_gfsr4 : *const gsl_rng_type;

    pub static gsl_rng_rand : *const gsl_rng_type;
    pub static gsl_rng_random_bsd : *const gsl_rng_type;
    pub static gsl_rng_random_libc5 : *const gsl_rng_type;
    pub static gsl_rng_random_glibc2 : *const gsl_rng_type;
    pub static gsl_rng_rand48 : *const gsl_rng_type;

    pub static gsl_rng_default : *const gsl_rng_type;
    pub static gsl_rng_ranf : *const gsl_rng_type;
    pub static gsl_rng_ranmar : *const gsl_rng_type;
    pub static gsl_rng_r250 : *const gsl_rng_type;
    pub static gsl_rng_tt800 : *const gsl_rng_type;
    pub static gsl_rng_vax : *const gsl_rng_type;
    pub static gsl_rng_transputer : *const gsl_rng_type;
    pub static gsl_rng_randu : *const gsl_rng_type;
    pub static gsl_rng_minstd : *const gsl_rng_type;
    pub static gsl_rng_uni : *const gsl_rng_type;
    pub static gsl_rng_uni32 : *const gsl_rng_type;
    pub static gsl_rng_slatec : *const gsl_rng_type;
    pub static gsl_rng_zuf : *const gsl_rng_type;
    pub static gsl_rng_knuthran2 : *const gsl_rng_type;
    pub static gsl_rng_knuthran2002 : *const gsl_rng_type;
    pub static gsl_rng_knuthran : *const gsl_rng_type;
    pub static gsl_rng_borosh13 : *const gsl_rng_type;
    pub static gsl_rng_fishman18 : *const gsl_rng_type;
    pub static gsl_rng_fishman20 : *const gsl_rng_type;
    pub static gsl_rng_lecuyer21 : *const gsl_rng_type;
    pub static gsl_rng_waterman14 : *const gsl_rng_type;
    pub static gsl_rng_fishman2x : *const gsl_rng_type;
    pub static gsl_rng_coveyou : *const gsl_rng_type;

    pub static gsl_rng_default_seed : c_ulong;

    pub static gsl_interp_linear : *const gsl_interp_type;
    pub static gsl_interp_polynomial : *const gsl_interp_type;
    pub static gsl_interp_cspline : *const gsl_interp_type;
    pub static gsl_interp_cspline_periodic : *const gsl_interp_type;
    pub static gsl_interp_akima : *const gsl_interp_type;
    pub static gsl_interp_akima_periodic : *const gsl_interp_type;

    pub static gsl_odeiv2_step_rk2 : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk4 : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rkf45 : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rkck : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk8pd : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk1imp : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk2imp : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk4imp : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_bsimp : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_msadams : *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_msbdf : *const gsl_odeiv2_step_type;

    pub static gsl_odeiv2_control_scaled : *const gsl_odeiv2_control_type;
    pub static gsl_odeiv2_control_standard : *const gsl_odeiv2_control_type;

    pub static gsl_qrng_niederreiter_2 : *const gsl_qrng_type;
    pub static gsl_qrng_sobol : *const gsl_qrng_type;
    pub static gsl_qrng_halton : *const gsl_qrng_type;
    pub static gsl_qrng_reversehalton : *const gsl_qrng_type;

    pub static gsl_wavelet_daubechies : *const gsl_wavelet_type;
    pub static gsl_wavelet_daubechies_centered : *const gsl_wavelet_type;
    pub static gsl_wavelet_haar : *const gsl_wavelet_type;
    pub static gsl_wavelet_haar_centered : *const gsl_wavelet_type;
    pub static gsl_wavelet_bspline : *const gsl_wavelet_type;
    pub static gsl_wavelet_bspline_centered : *const gsl_wavelet_type;

    // Airy functions
    pub fn gsl_sf_airy_Ai(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Ai_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_scaled_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_scaled_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    // Derivatives of Airy Functions
    pub fn gsl_sf_airy_Ai_deriv(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi_deriv(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Ai_deriv_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_scaled_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi_deriv_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_scaled_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    //  Zeros of Airy Functions
    pub fn gsl_sf_airy_zero_Ai(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Ai_e(s: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_zero_Bi(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Bi_e(s: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    // Zeros of Derivatives of Airy Functions
    pub fn gsl_sf_airy_zero_Ai_deriv(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Ai_deriv_e(s: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_zero_Bi_deriv(s: c_uint) -> c_double;
    pub fn gsl_sf_airy_zero_Bi_deriv_e(s: c_uint, result: *mut gsl_sf_result) -> enums::Value;

    // Bessel functions
    // Regular Modified Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_I0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_I1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_In(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_In_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_In_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_bessel_I0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_I1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_In_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_In_scaled_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_In_scaled_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    // Regular Modified Spherical Bessel Functions
    // The regular modified spherical Bessel functions i_l(x) are related to the modified Bessel functions of fractional order, i_l(x) = \sqrt{\pi/(2x)} I_{l+1/2}(x)
    pub fn gsl_sf_bessel_i0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_i1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_i2_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i2_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_il_scaled(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_il_scaled_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_il_scaled_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_bessel_Inu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Inu_scaled(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_scaled_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Regular Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_J0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_J1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Jn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Jn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    // Regular Spherical Bessel Functions
    pub fn gsl_sf_bessel_j0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_j1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_j2(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_jl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_jl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_jl_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_bessel_jl_steed_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    // Regular Bessel Function—Fractional Order
    pub fn gsl_sf_bessel_Jnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jnu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_sequence_Jnu_e(nu: c_double, mode: ::Mode, size: i64, v: *mut c_double) -> enums::Value;
    // Irregular Modified Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_K0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_K1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Kn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Kn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_bessel_K0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_K1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Kn_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_scaled_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Kn_scaled_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    // Irregular Modified Spherical Bessel Functions
    pub fn gsl_sf_bessel_k0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_k1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_k2_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k2_scaled_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_kl_scaled(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_kl_scaled_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_kl_scaled_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    // Irregular Modified Bessel Functions—Fractional Order
    pub fn gsl_sf_bessel_Knu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Knu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_lnKnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_lnKnu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Knu_scaled(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Knu_scaled_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Irregular Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_Y0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Y0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Y1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Y1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Yn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Yn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_Yn_array(nmin: c_int, nmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    // Irregular Spherical Bessel Functions
    pub fn gsl_sf_bessel_y0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_y1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_y2(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_yl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_yl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_yl_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    // Irregular Bessel Functions—Fractional Order
    pub fn gsl_sf_bessel_Ynu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Ynu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Zeros of Regular Bessel Functions
    pub fn gsl_sf_bessel_zero_J0(s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_J0_e(s: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_zero_J1(s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_J1_e(s: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_bessel_zero_Jnu(nu: f64, s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_Jnu_e(nu: f64, s: c_uint, result: *mut gsl_sf_result) -> enums::Value;

    // Trigonometric Functions
    pub fn gsl_sf_sin(x: c_double) -> c_double;
    pub fn gsl_sf_sin_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_cos(x: c_double) -> c_double;
    pub fn gsl_sf_cos_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hypot(x: c_double) -> c_double;
    pub fn gsl_sf_hypot_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_sinc(x: c_double) -> c_double;
    pub fn gsl_sf_sinc_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_complex_sin_e(zr: c_double, zi: c_double, szr: *mut gsl_sf_result, szi: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_complex_cos_e(zr: c_double, zi: c_double, czr: *mut gsl_sf_result, czi: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_complex_logsin_e(zr: c_double, zi: c_double, lszr: *mut gsl_sf_result, lszi: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lnsinh(x: c_double) -> c_double;
    pub fn gsl_sf_lnsinh_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lncosh(x: c_double) -> c_double;
    pub fn gsl_sf_lncosh_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_polar_to_rect(r: c_double, theta: c_double, x: *mut gsl_sf_result, y: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_rect_to_polar(x: c_double, y: c_double, r: *mut gsl_sf_result, theta: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_angle_restrict_symm(theta: c_double) -> c_double;
    pub fn gsl_sf_angle_restrict_symm_e(theta: *mut c_double) -> enums::Value;
    pub fn gsl_sf_angle_restrict_pos(theta: c_double) -> c_double;
    pub fn gsl_sf_angle_restrict_pos_e(theta: *mut c_double) -> enums::Value;
    pub fn gsl_sf_sin_err_e(x: c_double, dx: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_cos_err_e(x: c_double, dx: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Exponential Integrals functions
    pub fn gsl_sf_expint_E1(x: c_double) -> c_double;
    pub fn gsl_sf_expint_E1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_expint_E2(x: c_double) -> c_double;
    pub fn gsl_sf_expint_E2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_expint_En(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_expint_En_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_expint_Ei(x: c_double) -> c_double;
    pub fn gsl_sf_expint_Ei_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_Shi(x: c_double) -> c_double;
    pub fn gsl_sf_Shi_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_Chi(x: c_double) -> c_double;
    pub fn gsl_sf_Chi_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_expint_3(x: c_double) -> c_double;
    pub fn gsl_sf_expint_3_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_Si(x: c_double) -> c_double;
    pub fn gsl_sf_Si_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_Ci(x: c_double) -> c_double;
    pub fn gsl_sf_Ci_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_atanint(x: c_double) -> c_double;
    pub fn gsl_sf_atanint_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Clausen functions
    pub fn gsl_sf_clausen(x: c_double) -> c_double;
    pub fn gsl_sf_clausen_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Coulomb functions
    // Normalized Hydrogenic Bound States
    pub fn gsl_sf_hydrogenicR_1(Z: c_double, r: c_double) -> c_double;
    pub fn gsl_sf_hydrogenicR_1_e(Z: c_double, r: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hydrogenicR(n: c_int, l: c_int, Z: c_double, r: c_double) -> c_double;
    pub fn gsl_sf_hydrogenicR_e(n: c_int, l: c_int, Z: c_double, r: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Coulomb Wave Functions
    // The Coulomb wave functions F_L(\eta,x), G_L(\eta,x) are described in Abramowitz & Stegun, Chapter 14. Because there can be a large dynamic range of values for these functions, overflows are handled gracefully. If an overflow occurs, GSL_EOVRFLW is signalled and exponent(s) are returned through the modifiable parameters exp_F, exp_G. The full solution can be reconstructed from the following relations,
    //
    // F_L(eta,x)  =  fc[k_L] * exp(exp_F)
    // G_L(eta,x)  =  gc[k_L] * exp(exp_G)
    //
    // F_L'(eta,x) = fcp[k_L] * exp(exp_F)
    // G_L'(eta,x) = gcp[k_L] * exp(exp_G)
    pub fn gsl_sf_coulomb_wave_FG_e(eta: c_double, x: c_double, L_F: c_double, k: c_int, F: *mut gsl_sf_result, Fp: *mut gsl_sf_result,
        G: *mut gsl_sf_result, Gp: *mut gsl_sf_result, exp_F: *mut c_double, exp_G: *mut c_double) -> enums::Value;
    pub fn gsl_sf_coulomb_wave_F_array(L_min: c_double, kmax: c_int, eta: c_double, x: c_double, fc_array: *mut c_double,
        F_exponent: *mut c_double) -> enums::Value;
    pub fn gsl_sf_coulomb_wave_FG_array(L_min: c_double, kmax: c_int, eta: c_double, x: c_double, fc_array: *mut c_double,
        gc_array: *mut c_double, F_exponent: *mut c_double, G_exponent: *mut c_double) -> enums::Value;
    pub fn gsl_sf_coulomb_wave_FGp_array(L_min: c_double, kmax: c_int, eta: c_double, x: c_double, fc_array: *mut c_double,
        fcp_array: *mut c_double, gc_array: *mut c_double, gcp_array: *mut c_double, F_exponent: *mut c_double,
        G_exponent: *mut c_double) -> enums::Value;
    pub fn gsl_sf_coulomb_wave_sphF_array(L_min: c_double, kmax: c_int, eta: c_double, x: c_double, fc_array: *mut c_double,
        f_exponent: *mut c_double) -> enums::Value;
    // Coulomb Wave Function Normalization Constant
    pub fn gsl_sf_coulomb_CL_e(L: c_double, eta: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_coulomb_CL_array(Lmin: c_double, kmax: c_int, eta: c_double, cl: *mut c_double) -> enums::Value;

    // Coupling Coefficients functions
    pub fn gsl_sf_coupling_3j(two_ja: c_int, two_jb: c_int, two_jc: c_int, two_ma: c_int, two_mc: c_int, two_mc: c_int) -> c_double;
    pub fn gsl_sf_coupling_3j_e(two_ja: c_int, two_jb: c_int, two_jc: c_int, two_ma: c_int, two_mc: c_int, two_mc: c_int,
        result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_coupling_6j(two_ja: c_int, two_jb: c_int, two_jc: c_int, two_jd: c_int, two_je: c_int, two_jf: c_int) -> c_double;
    pub fn gsl_sf_coupling_6j_e(two_ja: c_int, two_jb: c_int, two_jc: c_int, two_jd: c_int, two_je: c_int, two_jf: c_int,
        result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_coupling_9j(two_ja: c_int, two_jb: c_int, two_jc: c_int, two_jd: c_int, two_je: c_int, two_jf: c_int,
        two_jg: c_int, two_jh: c_int, two_ji: c_int) -> c_double;
    pub fn gsl_sf_coupling_9j_e(two_ja: c_int, two_jb: c_int, two_jc: c_int, two_jd: c_int, two_je: c_int, two_jf: c_int,
        two_jg: c_int, two_jh: c_int, two_ji: c_int, result: *mut gsl_sf_result) -> enums::Value;

    // Dawson functions
    pub fn gsl_sf_dawson(x: c_double) -> c_double;
    pub fn gsl_sf_dawson_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Debye functions
    pub fn gsl_sf_debye_1(x: c_double) -> c_double;
    pub fn gsl_sf_debye_1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_debye_2(x: c_double) -> c_double;
    pub fn gsl_sf_debye_2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_debye_3(x: c_double) -> c_double;
    pub fn gsl_sf_debye_3_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_debye_4(x: c_double) -> c_double;
    pub fn gsl_sf_debye_4_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_debye_5(x: c_double) -> c_double;
    pub fn gsl_sf_debye_5_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_debye_6(x: c_double) -> c_double;
    pub fn gsl_sf_debye_6_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Dilogarithm functions
    // real argument
    pub fn gsl_sf_dilog(x: c_double) -> c_double;
    pub fn gsl_sf_dilog_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // complex argument
    pub fn gsl_sf_complex_dilog_e(r: c_double, theta: c_double, result: *mut gsl_sf_result, result_im: *mut gsl_sf_result) -> enums::Value;

    // Elementary Operations functions
    pub fn gsl_sf_multiply_e(x: c_double, y: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_multiply_err_e(x: c_double, dx: c_double, y: c_double, dy: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Elliptic functions (Jacobi)
    pub fn gsl_sf_elljac_e(u: c_double, m: c_double, sn: *mut c_double, cn: *mut c_double, dn: *mut c_double) -> enums::Value;

    // Error functions
    pub fn gsl_sf_erf(x: c_double) -> c_double;
    pub fn gsl_sf_erf_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Complementary Error functions
    pub fn gsl_sf_erfc(x: c_double) -> c_double;
    pub fn gsl_sf_erfc_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Log Complementary Error functions
    pub fn gsl_sf_log_erfc(x: c_double) -> c_double;
    pub fn gsl_sf_log_erfc_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Probability functions
    // The probability functions for the Normal or Gaussian distribution are described in Abramowitz & Stegun, Section 26.2.
    pub fn gsl_sf_erf_Z(x: c_double) -> c_double;
    pub fn gsl_sf_erf_Z_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_erf_Q(x: c_double) -> c_double;
    pub fn gsl_sf_erf_Q_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hazard(x: c_double) -> c_double;
    pub fn gsl_sf_hazard_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Exponential functions
    pub fn gsl_sf_exp(x: c_double) -> c_double;
    pub fn gsl_sf_exp_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_exp_e10_e(x: c_double, result: *mut gsl_sf_result_e10) -> enums::Value;
    pub fn gsl_sf_exp_mult(x: c_double, y: c_double) -> c_double;
    pub fn gsl_sf_exp_mult_e(x: c_double, y: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_exp_mult_e10_e(x: c_double, y: c_double, result: *mut gsl_sf_result_e10) -> enums::Value;
    // Relative Exponential functions
    pub fn gsl_sf_expm1(x: c_double) -> c_double;
    pub fn gsl_sf_expm1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_exprel(x: c_double) -> c_double;
    pub fn gsl_sf_exprel_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_exprel_2(x: c_double) -> c_double;
    pub fn gsl_sf_exprel_2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_exprel_n(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_exprel_n_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Exponentiation With Error Estimate
    pub fn gsl_sf_exp_err_e(x: c_double, dx: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_exp_err_e10_e(x: c_double, dx: c_double, result: *mut gsl_sf_result_e10) -> enums::Value;
    pub fn gsl_sf_exp_mult_err_e(x: c_double, dx: c_double, y: c_double, dy: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_exp_mult_err_e10_e(x: c_double, dx: c_double, y: c_double, dy: c_double, result: *mut gsl_sf_result_e10) -> enums::Value;

    // Gamma Beta functions
    // Gamma functions
    pub fn gsl_sf_gamma(x: c_double) -> c_double;
    pub fn gsl_sf_gamma_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lngamma(x: c_double) -> c_double;
    pub fn gsl_sf_lngamma_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lngamma_sgn_e(x: c_double, result_lg: *mut gsl_sf_result, sgn: *mut c_double) -> enums::Value;
    pub fn gsl_sf_gammastar(x: c_double) -> c_double;
    pub fn gsl_sf_gammastar_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_gammainv(x: c_double) -> c_double;
    pub fn gsl_sf_gammainv_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lngamma_complex_e(zr: c_double, zi: c_double, lnr: *mut gsl_sf_result, arg: *mut gsl_sf_result) -> enums::Value;
    // Factorials
    pub fn gsl_sf_fact(n: c_uint) -> c_double;
    pub fn gsl_sf_fact_e(n: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_doublefact(n: c_uint) -> c_double;
    pub fn gsl_sf_doublefact_e(n: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lnfact(n: c_uint) -> c_double;
    pub fn gsl_sf_lnfact_e(n: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lndoublefact(n: c_uint) -> c_double;
    pub fn gsl_sf_lndoublefact_e(n: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_choose(n: c_uint, m: c_uint) -> c_double;
    pub fn gsl_sf_choose_e(n: c_uint, m: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lnchoose(n: c_uint, m: c_uint) -> c_double;
    pub fn gsl_sf_lnchoose_e(n: c_uint, m: c_uint, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_taylorcoeff(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_taylorcoeff_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Pochhammer Symbol
    pub fn gsl_sf_poch(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_poch_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lnpoch(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_lnpoch_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lnpoch_sgn_e(a: c_double, x: c_double, result: *mut gsl_sf_result, sgn: *mut c_double) -> enums::Value;
    pub fn gsl_sf_pochrel(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_pochrel_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Beta functions
    pub fn gsl_sf_beta(a: c_double, b: c_double) -> c_double;
    pub fn gsl_sf_beta_e(a: c_double, b: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lnbeta(a: c_double, b: c_double) -> c_double;
    pub fn gsl_sf_lnbeta_e(a: c_double, b: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Incomplete Gamma functions
    pub fn gsl_sf_gamma_inc(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gamma_inc_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_gamma_inc_Q(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gamma_inc_Q_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_gamma_inc_P(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gamma_inc_P_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Incomplete Beta functions
    pub fn gsl_sf_beta_inc(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_beta_inc_e(a: c_double, b: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Gegenbauer functions
    pub fn gsl_sf_gegenpoly_1(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_2(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_3(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_1_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_gegenpoly_2_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_gegenpoly_3_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_gegenpoly_n(n: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_n_e(n: c_int, lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_gegenpoly_array(nmax: c_int, lambda: c_double, x: c_double, result_array: *mut c_double) -> enums::Value;

    // Hypergeometric functions
    pub fn gsl_sf_hyperg_0F1(c: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_0F1_e(c: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_1F1_int(m: c_int, n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_1F1_int_e(m: c_int, n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_1F1(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_1F1_e(a: c_double, b: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_U_int(m: c_int, n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_U_int_e(m: c_int, n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_U_int_e10_e(m: c_int, n: c_int, x: c_double, result: *mut gsl_sf_result_e10) -> enums::Value;
    pub fn gsl_sf_hyperg_U(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_U_e(a: c_double, b: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_U_e10_e(a: c_double, b: c_double, x: c_double, result: *mut gsl_sf_result_e10) -> enums::Value;
    pub fn gsl_sf_hyperg_2F1(a: c_double, b: c_double, c: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_2F1_e(a: c_double, b: c_double, c: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_2F1_conj(aR: c_double, aI: c_double, c: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_2F1_conj_e(aR: c_double, aI: c_double, c: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_2F1_renorm(a: c_double, b: c_double, c: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_2F1_renorm_e(a: c_double, b: c_double, c: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_2F1_conj_renorm(aR: c_double, aI: c_double, c: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_2F1_conj_renorm_e(aR: c_double, aI: c_double, c: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_hyperg_2F0(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_2F0_e(a: c_double, b: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    /// Laguerre functions
    pub fn gsl_sf_laguerre_1(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_2(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_3(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_1_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_laguerre_2_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_laguerre_3_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_laguerre_n(n: c_int, a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_n_e(n: c_int, a: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Lambert W functions
    pub fn gsl_sf_lambert_W0(x: c_double) -> c_double;
    pub fn gsl_sf_lambert_W0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_lambert_Wm1(x: c_double) -> c_double;
    pub fn gsl_sf_lambert_Wm1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Legendre functions
    // Legendre Polynomials
    pub fn gsl_sf_legendre_P1(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_P2(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_P3(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_P1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_P2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_P3_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_Pl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Pl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_Pl_array(lmax: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_legendre_Pl_deriv_array(lmax: c_int, x: c_double, result_array: *mut c_double,
        result_deriv_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_legendre_Q0(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Q0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_Q1(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Q1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_Ql(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Ql_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Associated Legendre Polynomials and Spherical Harmonics
    pub fn gsl_sf_legendre_Plm(l: c_int, m: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Plm_e(l: c_int, m: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_Plm_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_Plm_deriv_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double,
        result_deriv_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_legendre_sphPlm(l: c_int, m: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_sphPlm_e(l: c_int, m: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_sphPlm_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_sphPlm_deriv_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double,
        result_deriv_array: *mut c_double) -> enums::Value;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_array_size(lmax: c_int, m: c_int) -> enums::Value;
    // Conical functions
    pub fn gsl_sf_conicalP_half(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_half_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_conicalP_mhalf(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_mhalf_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_conicalP_0(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_0_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_conicalP_1(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_1_e(lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_conicalP_sph_reg(l: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_sph_reg_e(l: c_int, lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_conicalP_cyl_reg(m: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_cyl_reg_e(m: c_int, lambda: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Radial Functions for Hyperbolic Space
    pub fn gsl_sf_legendre_H3d_0(lambda: c_double, eta: c_double) -> c_double;
    pub fn gsl_sf_legendre_H3d_0_e(lambda: c_double, eta: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_H3d_1(lambda: c_double, eta: c_double) -> c_double;
    pub fn gsl_sf_legendre_H3d_1_e(lambda: c_double, eta: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_H3d(l: c_int, lambda: c_double, eta: c_double) -> c_double;
    pub fn gsl_sf_legendre_H3d_e(l: c_int, lambda: c_double, eta: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_H3d_array(lmax: c_int, lambda: c_double, eta: c_double, result_array: *mut c_double) -> enums::Value;

    // Logarithm and Related Functions
    pub fn gsl_sf_log(x: c_double) -> c_double;
    pub fn gsl_sf_log_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_log_abs(x: c_double) -> c_double;
    pub fn gsl_sf_log_abs_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_complex_log_e(zr: c_double, zi: c_double, lnr: *mut gsl_sf_result, theta: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_log_1plusx(x: c_double) -> c_double;
    pub fn gsl_sf_log_1plusx_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_log_1plusx_mx(x: c_double) -> c_double;
    pub fn gsl_sf_log_1plusx_mx_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Power functions
    pub fn gsl_sf_pow_int(x: c_double, n: c_int) -> c_double;
    pub fn gsl_sf_pow_int_e(x: c_double, n: c_int, result: *mut gsl_sf_result) -> enums::Value;

    // Psi (Digamma) functions
    // Digamma functions
    pub fn gsl_sf_psi_int(n: c_int) -> c_double;
    pub fn gsl_sf_psi_int_e(n: c_int, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_psi(x: c_double) -> c_double;
    pub fn gsl_sf_psi_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_psi_1piy(y: c_double) -> c_double;
    pub fn gsl_sf_psi_1piy_e(y: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Trigamma functions
    pub fn gsl_sf_psi_1_int(n: c_int) -> c_double;
    pub fn gsl_sf_psi_1_int_e(n: c_int, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_psi_1(x: c_double) -> c_double;
    pub fn gsl_sf_psi_1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Polygamma functions
    pub fn gsl_sf_psi_n(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_psi_n_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Synchrotron functions
    pub fn gsl_sf_synchrotron_1(x: c_double) -> c_double;
    pub fn gsl_sf_synchrotron_1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_synchrotron_2(x: c_double) -> c_double;
    pub fn gsl_sf_synchrotron_2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Transport functions
    pub fn gsl_sf_transport_2(x: c_double) -> c_double;
    pub fn gsl_sf_transport_2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_transport_3(x: c_double) -> c_double;
    pub fn gsl_sf_transport_3_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_transport_4(x: c_double) -> c_double;
    pub fn gsl_sf_transport_4_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_transport_5(x: c_double) -> c_double;
    pub fn gsl_sf_transport_5_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Zeta functions
    // Riemann Zeta functions
    pub fn gsl_sf_zeta_int(n: c_int) -> c_double;
    pub fn gsl_sf_zeta_int_e(n: c_int, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_zeta(s: c_double) -> c_double;
    pub fn gsl_sf_zeta_e(s: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Riemann Zeta functions Minus One
    pub fn gsl_sf_zetam1_int(n: c_int) -> c_double;
    pub fn gsl_sf_zetam1_int_e(n: c_int, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_zetam1(s: c_double) -> c_double;
    pub fn gsl_sf_zetam1_e(s: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Hurwitz Zeta functions
    pub fn gsl_sf_hzeta(s: c_double, q: c_double) -> c_double;
    pub fn gsl_sf_hzeta_e(s: c_double, q: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Eta functions
    pub fn gsl_sf_eta_int(n: c_int) -> c_double;
    pub fn gsl_sf_eta_int_e(n: c_int, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_eta(s: c_double) -> c_double;
    pub fn gsl_sf_eta_e(s: c_double, result: *mut gsl_sf_result) -> enums::Value;

    // Elliptic Integrals
    // Legendre Form of Complete Elliptic Integrals
    pub fn gsl_sf_ellint_Kcomp(k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_Kcomp_e(k: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_Ecomp(k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_Ecomp_e(k: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_Pcomp(k: c_double, n: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_Pcomp_e(k: c_double, n: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    // Legendre Form of Incomplete Elliptic Integrals
    pub fn gsl_sf_ellint_F(phi: c_double, k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_F_e(phi: c_double, k: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_E(phi: c_double, k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_E_e(phi: c_double, k: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_P(phi: c_double, k: c_double, n: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_P_e(phi: c_double, k: c_double, n: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_D(phi: c_double, k: c_double, n: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_D_e(phi: c_double, k: c_double, n: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    // Carlson Forms
    pub fn gsl_sf_ellint_RC(x: c_double, y: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_RC_e(x: c_double, y: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_RD(x: c_double, y: c_double, z: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_RD_e(x: c_double, y: c_double, z: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_RF(x: c_double, y: c_double, z: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_RF_e(x: c_double, y: c_double, z: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_ellint_RJ(x: c_double, y: c_double, z: c_double, p: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_RJ_e(x: c_double, y: c_double, z: c_double, p: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> enums::Value;

    // Fermi-Dirac functions
    // Complete Fermi-Dirac Integrals
    pub fn gsl_sf_fermi_dirac_m1(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_m1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_fermi_dirac_0(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_0_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_fermi_dirac_1(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_1_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_fermi_dirac_2(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_2_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_fermi_dirac_int(j: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_int_e(j: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_fermi_dirac_mhalf(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_mhalf_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_fermi_dirac_half(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_half_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_fermi_dirac_3half(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_3half_e(x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    // Incomplete Fermi-Dirac Integrals
    pub fn gsl_sf_fermi_dirac_inc_0(x: c_double, b: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_inc_0_e(x: c_double, b: c_double, result: *mut gsl_sf_result) -> enums::Value;

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
    pub fn cblas_sgemv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, alpha: c_float,
        A: *const c_float, lda: c_int, x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_sgbmv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_strmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stbmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stpmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_float, x: *mut c_float, incx: c_int);
    pub fn cblas_strsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stbsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_float, lda: c_int, x: *mut c_float, incx: c_int);
    pub fn cblas_stpsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_float, x: *mut c_float, incx: c_int);
    pub fn cblas_dgemv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, alpha: c_double,
        A: *const c_double, lda: c_int, x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dgbmv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dtrmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtbmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtpmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_double, x: *mut c_double, incx: c_int);
    pub fn cblas_dtrsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtbsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_double, lda: c_int, x: *mut c_double, incx: c_int);
    pub fn cblas_dtpsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_double, x: *mut c_double, incx: c_int);
    pub fn cblas_cgemv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, alpha: *const c_void,
        A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_cgbmv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_ctrmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctbmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctpmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_ctrsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctbsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ctpsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_zgemv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, alpha: *const c_void,
        A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zgbmv(order: cblas::Order, transA: cblas::Transpose, M: c_int, N: c_int, KL: c_int, KU: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_ztrmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztbmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztpmv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_ztrsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztbsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, K: c_int, A: *const c_void, lda: c_int, x: *mut c_void, incx: c_int);
    pub fn cblas_ztpsv(order: cblas::Order, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, N: c_int, Ap: *const c_void, x: *mut c_void, incx: c_int);
    pub fn cblas_ssymv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_float, A: *const c_float, lda: c_int,
        x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_ssbmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, K: c_int, alpha: c_float, A: *const c_float,
        lda: c_int, x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_sspmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_float, Ap: *const c_float,
        x: *const c_float, incx: c_int, beta: c_float, y: *mut c_float, incy: c_int);
    pub fn cblas_sger(order: cblas::Order, M: c_int, N: c_int, alpha: c_float, x: *const c_float, incx: c_int, y: *const c_float,
        incy: c_int, A: *mut c_float, lda: c_int);
    pub fn cblas_ssyr(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        A: *mut c_float, lda: c_int);
    pub fn cblas_sspr(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        Ap: *mut c_float);
    pub fn cblas_ssyr2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        y: *const c_float, incy: c_int, A: *mut c_float, lda: c_int);
    pub fn cblas_sspr2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_float, x: *const c_float, incx: c_int,
        y: *const c_float, incy: c_int, A: *mut c_float);
    pub fn cblas_dsymv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_double, A: *const c_double, lda: c_int,
        x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dsbmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, K: c_int, alpha: c_double, A: *const c_double,
        lda: c_int, x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dspmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_double, Ap: *const c_double,
        x: *const c_double, incx: c_int, beta: c_double, y: *mut c_double, incy: c_int);
    pub fn cblas_dger(order: cblas::Order, M: c_int, N: c_int, alpha: c_double, x: *const c_double, incx: c_int, y: *const c_double,
        incy: c_int, A: *mut c_double, lda: c_int);
    pub fn cblas_dsyr(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        A: *mut c_double, lda: c_int);
    pub fn cblas_dspr(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        Ap: *mut c_double);
    pub fn cblas_dsyr2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        y: *const c_double, incy: c_int, A: *mut c_double, lda: c_int);
    pub fn cblas_dspr2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: c_double, x: *const c_double, incx: c_int,
        y: *const c_double, incy: c_int, A: *mut c_double);
    pub fn cblas_csymv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_csbmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, K: c_int, alpha: *const c_void, A: *const c_void,
        lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_cspmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, Ap: *const c_void,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_cgeru(order: cblas::Order, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_cgerc(order: cblas::Order, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_cher(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        A: *mut c_void, lda: c_int);
    pub fn cblas_chpr(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        Ap: *mut c_void);
    pub fn cblas_cher2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_chpr2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, Ap: *mut c_void);
    pub fn cblas_zsymv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zsbmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, K: c_int, alpha: *const c_void, A: *const c_void,
        lda: c_int, x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zspmv(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, Ap: *const c_void,
        x: *const c_void, incx: c_int, beta: *const c_void, y: *mut c_void, incy: c_int);
    pub fn cblas_zgeru(order: cblas::Order, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_zgerc(order: cblas::Order, M: c_int, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int, y: *const c_void,
        incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_zher(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        A: *mut c_void, lda: c_int);
    pub fn cblas_zhpr(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        Ap: *mut c_void);
    pub fn cblas_zher2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, A: *mut c_void, lda: c_int);
    pub fn cblas_zhpr2(order: cblas::Order, uplo: cblas::Uplo, N: c_int, alpha: *const c_void, x: *const c_void, incx: c_int,
        y: *const c_void, incy: c_int, Ap: *mut c_void);
    // Level 3 CBLAS functions
    pub fn cblas_sgemm(order: cblas::Order, transA: cblas::Transpose, transB: cblas::Transpose, M: c_int, N: c_int,
        K: c_int, alpha: c_float, A: *const c_float, lda: c_int, B: *const c_float, ldb: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_ssymm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, M: c_int, N: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, B: *const c_float, ldb: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_ssyrk(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_ssyr2k(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: c_float, A: *const c_float, lda: c_int, B: *const c_float, ldb: c_int, beta: c_float, C: *mut c_float, ldc: c_int);
    pub fn cblas_strmm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: c_float, A: *const c_float, lda: c_int, B: *mut c_float, ldb: c_int);
    pub fn cblas_strsm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: c_float, A: *const c_float, lda: c_int, B: *mut c_float, ldb: c_int);
    pub fn cblas_dgemm(order: cblas::Order, transA: cblas::Transpose, transB: cblas::Transpose, M: c_int, N: c_int,
        K: c_int, alpha: c_double, A: *const c_double, lda: c_int, B: *const c_double, ldb: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dsymm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, M: c_int, N: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, B: *const c_double, ldb: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dsyrk(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dsyr2k(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: c_double, A: *const c_double, lda: c_int, B: *const c_double, ldb: c_int, beta: c_double, C: *mut c_double, ldc: c_int);
    pub fn cblas_dtrmm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: c_double, A: *const c_double, lda: c_int, B: *mut c_double, ldb: c_int);
    pub fn cblas_dtrsm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: c_double, A: *const c_double, lda: c_int, B: *mut c_double, ldb: c_int);
    pub fn cblas_cgemm(order: cblas::Order, transA: cblas::Transpose, transB: cblas::Transpose, M: c_int, N: c_int,
        K: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void,
        ldc: c_int);
    pub fn cblas_csymm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_csyrk(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_csyr2k(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_ctrmm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_ctrsm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_zgemm(order: cblas::Order, transA: cblas::Transpose, transB: cblas::Transpose, M: c_int, N: c_int,
        K: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void,
        ldc: c_int);
    pub fn cblas_zsymm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zsyrk(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zsyr2k(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_ztrmm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_ztrsm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, transA: cblas::Transpose,
        diag: cblas::Diag, M: c_int, N: c_int, alpha: *const c_void, A: *const c_void, lda: c_int, B: *mut c_void, ldb: c_int);
    pub fn cblas_chemm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_cherk(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_cher2k(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zhemm(order: cblas::Order, side: cblas::Side, uplo: cblas::Uplo, M: c_int, N: c_int,
        alpha: *const c_void, A: *const c_void, lda: c_int, B: *const c_void, ldb: c_int, beta: *const c_void, C: *mut c_void, ldc: c_int);
    pub fn cblas_zherk(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
        alpha: c_double, A: *const c_void, lda: c_int, beta: c_double, C: *mut c_void, ldc: c_int);
    pub fn cblas_zher2k(order: cblas::Order, uplo: cblas::Uplo, trans: cblas::Transpose, N: c_int, K: c_int,
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
    pub fn gsl_vector_alloc(size: size_t) -> *mut gsl_vector;
    pub fn gsl_vector_calloc(size: size_t) -> *mut gsl_vector;
    pub fn gsl_vector_free(vector: *mut gsl_vector);
    pub fn gsl_vector_get(vector: *const gsl_vector, i: size_t) -> c_double;
    pub fn gsl_vector_set(vector: *mut gsl_vector, i: size_t, x: c_double);
    pub fn gsl_vector_set_all(vector: *mut gsl_vector, x: c_double);
    pub fn gsl_vector_set_zero(vector: *mut gsl_vector);
    pub fn gsl_vector_set_basis(vector: *mut gsl_vector, i: size_t);
    pub fn gsl_vector_memcpy(dest: *mut gsl_vector, src: *const gsl_vector) -> enums::Value;
    pub fn gsl_vector_swap(v: *mut gsl_vector, w: *mut gsl_vector) -> enums::Value;
    pub fn gsl_vector_swap_elements(vector: *mut gsl_vector, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_vector_reverse(vector: *mut gsl_vector) -> enums::Value;
    pub fn gsl_vector_add(dest: *mut gsl_vector, src: *const gsl_vector) -> enums::Value;
    pub fn gsl_vector_sub(dest: *mut gsl_vector, src: *const gsl_vector) -> enums::Value;
    pub fn gsl_vector_mul(dest: *mut gsl_vector, src: *const gsl_vector) -> enums::Value;
    pub fn gsl_vector_div(dest: *mut gsl_vector, src: *const gsl_vector) -> enums::Value;
    pub fn gsl_vector_scale(dest: *mut gsl_vector, x: c_double) -> enums::Value;
    pub fn gsl_vector_add_constant(dest: *mut gsl_vector, x: c_double) -> enums::Value;
    pub fn gsl_vector_max(vector: *const gsl_vector) -> c_double;
    pub fn gsl_vector_min(vector: *const gsl_vector) -> c_double;
    pub fn gsl_vector_minmax(vector: *const gsl_vector, min_out: *mut c_double, max_out: *mut c_double);
    pub fn gsl_vector_max_index(vector: *const gsl_vector) -> size_t;
    pub fn gsl_vector_min_index(vector: *const gsl_vector) -> size_t;
    pub fn gsl_vector_minmax_index(vector: *const gsl_vector, imin: *mut size_t, imax: *mut size_t);
    pub fn gsl_vector_isnull(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_ispos(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_isneg(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_isnonneg(vector: *const gsl_vector) -> c_int;
    pub fn gsl_vector_equal(u: *const gsl_vector, v: *const gsl_vector) -> c_int;

    // Vector views
    pub fn gsl_vector_subvector(v: *mut gsl_vector, offset: size_t, n: size_t) -> gsl_vector_view;
    pub fn gsl_vector_subvector_with_stride(v: *mut gsl_vector, offset: size_t, stride: size_t, n: size_t) -> gsl_vector_view;
    //pub fn gsl_vector_complex_real(v: *mut gsl_vector_complex) -> gsl_vector_view;
    //pub fn gsl_vector_complex_imag(v: *mut gsl_vector_complex) -> gsl_vector_view;
    pub fn gsl_vector_view_array(base: *mut c_double, n: size_t) -> gsl_vector_view;
    pub fn gsl_vector_view_array_with_stride(base: *mut c_double, stride: size_t, n: size_t) -> gsl_vector_view;

    // VectorComplex functions
    pub fn gsl_vector_complex_alloc(size: size_t) -> *mut gsl_vector_complex;
    pub fn gsl_vector_complex_calloc(size: size_t) -> *mut gsl_vector_complex;
    pub fn gsl_vector_complex_free(vector: *mut gsl_vector_complex);
    pub fn gsl_vector_complex_get(vector: *const gsl_vector_complex, i: size_t) -> gsl_complex;
    pub fn gsl_vector_complex_set(vector: *mut gsl_vector_complex, i: size_t, x: gsl_complex);
    pub fn gsl_vector_complex_set_all(vector: *mut gsl_vector_complex, x: gsl_complex);
    pub fn gsl_vector_complex_set_zero(vector: *mut gsl_vector_complex);
    pub fn gsl_vector_complex_set_basis(vector: *mut gsl_vector_complex, i: size_t);
    pub fn gsl_vector_complex_memcpy(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> enums::Value;
    pub fn gsl_vector_complex_swap(v: *mut gsl_vector_complex, w: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_vector_complex_swap_elements(vector: *mut gsl_vector_complex, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_vector_complex_reverse(vector: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_vector_complex_add(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> enums::Value;
    pub fn gsl_vector_complex_sub(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> enums::Value;
    pub fn gsl_vector_complex_mul(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> enums::Value;
    pub fn gsl_vector_complex_div(dest: *mut gsl_vector_complex, src: *const gsl_vector_complex) -> enums::Value;
    pub fn gsl_vector_complex_scale(dest: *mut gsl_vector_complex, x: gsl_complex) -> enums::Value;
    pub fn gsl_vector_complex_add_constant(dest: *mut gsl_vector_complex, x: gsl_complex) -> enums::Value;
    pub fn gsl_vector_complex_isnull(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_ispos(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_isneg(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_isnonneg(vector: *const gsl_vector_complex) -> c_int;
    pub fn gsl_vector_complex_equal(u: *const gsl_vector_complex, v: *const gsl_vector_complex) -> c_int;

    // VectorFloat functions
    pub fn gsl_vector_float_alloc(size: size_t) -> *mut gsl_vector_float;
    pub fn gsl_vector_float_calloc(size: size_t) -> *mut gsl_vector_float;
    pub fn gsl_vector_float_free(vector: *mut gsl_vector_float);
    pub fn gsl_vector_float_get(vector: *const gsl_vector_float, i: size_t) -> c_float;
    pub fn gsl_vector_float_set(vector: *mut gsl_vector_float, i: size_t, x: c_float);
    pub fn gsl_vector_float_set_all(vector: *mut gsl_vector_float, x: c_float);
    pub fn gsl_vector_float_set_zero(vector: *mut gsl_vector_float);
    pub fn gsl_vector_float_set_basis(vector: *mut gsl_vector_float, i: size_t);
    pub fn gsl_vector_float_memcpy(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> enums::Value;
    pub fn gsl_vector_float_swap(v: *mut gsl_vector_float, w: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_vector_float_swap_elements(vector: *mut gsl_vector_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_vector_float_reverse(vector: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_vector_float_add(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> enums::Value;
    pub fn gsl_vector_float_sub(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> enums::Value;
    pub fn gsl_vector_float_mul(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> enums::Value;
    pub fn gsl_vector_float_div(dest: *mut gsl_vector_float, src: *const gsl_vector_float) -> enums::Value;
    pub fn gsl_vector_float_scale(dest: *mut gsl_vector_float, x: c_float) -> enums::Value;
    pub fn gsl_vector_float_add_constant(dest: *mut gsl_vector_float, x: c_float) -> enums::Value;
    pub fn gsl_vector_float_max(vector: *const gsl_vector_float) -> c_float;
    pub fn gsl_vector_float_min(vector: *const gsl_vector_float) -> c_float;
    pub fn gsl_vector_float_minmax(vector: *const gsl_vector_float, min_out: *mut c_float, max_out: *mut c_float);
    pub fn gsl_vector_float_max_index(vector: *const gsl_vector_float) -> size_t;
    pub fn gsl_vector_float_min_index(vector: *const gsl_vector_float) -> size_t;
    pub fn gsl_vector_float_minmax_index(vector: *const gsl_vector_float, imin: *mut size_t, imax: *mut size_t);
    pub fn gsl_vector_float_isnull(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_ispos(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_isneg(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_isnonneg(vector: *const gsl_vector_float) -> c_int;
    pub fn gsl_vector_float_equal(u: *const gsl_vector_float, v: *const gsl_vector_float) -> c_int;

    // VectorComplexFloat functions
    pub fn gsl_vector_complex_float_alloc(size: size_t) -> *mut gsl_vector_complex_float;
    pub fn gsl_vector_complex_float_calloc(size: size_t) -> *mut gsl_vector_complex_float;
    pub fn gsl_vector_complex_float_free(vector: *mut gsl_vector_complex_float);
    pub fn gsl_vector_complex_float_get(vector: *const gsl_vector_complex_float, i: size_t) -> gsl_complex_float;
    pub fn gsl_vector_complex_float_set(vector: *mut gsl_vector_complex_float, i: size_t, x: gsl_complex_float);
    pub fn gsl_vector_complex_float_set_all(vector: *mut gsl_vector_complex_float, x: gsl_complex_float);
    pub fn gsl_vector_complex_float_set_zero(vector: *mut gsl_vector_complex_float);
    pub fn gsl_vector_complex_float_set_basis(vector: *mut gsl_vector_complex_float, i: size_t);
    pub fn gsl_vector_complex_float_memcpy(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_swap(v: *mut gsl_vector_complex_float, w: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_swap_elements(vector: *mut gsl_vector_complex_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_vector_complex_float_reverse(vector: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_add(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_sub(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_mul(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_div(dest: *mut gsl_vector_complex_float, src: *const gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_scale(dest: *mut gsl_vector_complex_float, x: gsl_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_add_constant(dest: *mut gsl_vector_complex_float, x: gsl_complex_float) -> enums::Value;
    pub fn gsl_vector_complex_float_isnull(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_ispos(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_isneg(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_isnonneg(vector: *const gsl_vector_complex_float) -> c_int;
    pub fn gsl_vector_complex_float_equal(u: *const gsl_vector_complex_float, v: *const gsl_vector_complex_float) -> c_int;

    // Matrix functions
    pub fn gsl_matrix_alloc(size1: size_t, size2: size_t) -> *mut gsl_matrix;
    pub fn gsl_matrix_calloc(size1: size_t, size2: size_t) -> *mut gsl_matrix;
    pub fn gsl_matrix_free(m: *mut gsl_matrix);
    pub fn gsl_matrix_get(m: *const gsl_matrix, i: size_t, j: size_t) -> c_double;
    pub fn gsl_matrix_set(m: *mut gsl_matrix, i: size_t, j: size_t, x: c_double);
    pub fn gsl_matrix_set_all(m: *mut gsl_matrix, x: c_double);
    pub fn gsl_matrix_set_zero(m: *mut gsl_matrix);
    pub fn gsl_matrix_set_identity(m: *mut gsl_matrix);
    pub fn gsl_matrix_memcpy(dest: *mut gsl_matrix, src: *const gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_swap(m: *mut gsl_matrix, w: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_get_row(vector: *mut gsl_vector, m: *const gsl_matrix, i: size_t) -> enums::Value;
    pub fn gsl_matrix_get_col(vector: *mut gsl_vector, m: *const gsl_matrix, j: size_t) -> enums::Value;
    pub fn gsl_matrix_set_row(m: *mut gsl_matrix, i: size_t, v: *const gsl_vector) -> enums::Value;
    pub fn gsl_matrix_set_col(m: *mut gsl_matrix, j: size_t, v: *const gsl_vector) -> enums::Value;
    pub fn gsl_matrix_swap_rows(m: *mut gsl_matrix, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_swap_columns(m: *mut gsl_matrix, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_swap_rowcol(m: *mut gsl_matrix, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_transpose_memcpy(dest: *mut gsl_matrix, src: *const gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_transpose(m: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_add(dest: *mut gsl_matrix, src: *const gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_sub(dest: *mut gsl_matrix, src: *const gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_mul_elements(dest: *mut gsl_matrix, src: *const gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_div_elements(dest: *mut gsl_matrix, src: *const gsl_matrix) -> enums::Value;
    pub fn gsl_matrix_scale(dest: *mut gsl_matrix, x: c_double) -> enums::Value;
    pub fn gsl_matrix_add_constant(dest: *mut gsl_matrix, x: c_double) -> enums::Value;
    pub fn gsl_matrix_max(m: *const gsl_matrix) -> c_double;
    pub fn gsl_matrix_min(m: *const gsl_matrix) -> c_double;
    pub fn gsl_matrix_minmax(m: *const gsl_matrix, min_out: *mut c_double, max_out: *mut c_double);
    pub fn gsl_matrix_max_index(m: *const gsl_matrix, imax: *mut size_t, jmax: *mut size_t);
    pub fn gsl_matrix_min_index(m: *const gsl_matrix, imin: *mut size_t, jmin: *mut size_t);
    pub fn gsl_matrix_minmax_index(m: *const gsl_matrix, imin: *mut size_t, jmin: *mut size_t, imax: *mut size_t, jmax: *mut size_t);
    pub fn gsl_matrix_isnull(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_ispos(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_isneg(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_isnonneg(m: *const gsl_matrix) -> c_int;
    pub fn gsl_matrix_equal(u: *const gsl_matrix, v: *const gsl_matrix) -> c_int;

    // Matrix views
    pub fn gsl_matrix_submatrix(m: *mut gsl_matrix, k1: size_t, k2: size_t, n1: size_t, n2: size_t) -> gsl_matrix_view;
    pub fn gsl_matrix_view_array(base: *mut c_double, n1: size_t, n2: size_t) -> gsl_matrix_view;
    pub fn gsl_matrix_view_array_with_tda(base: *mut c_double, n1: size_t, n2: size_t, tda: size_t) -> gsl_matrix_view;
    pub fn gsl_matrix_view_vector(v: *mut gsl_vector, n1: size_t, n2: size_t) -> gsl_matrix_view;
    pub fn gsl_matrix_view_vector_with_tda(v: *mut gsl_vector, n1: size_t, n2: size_t, tda: size_t) -> gsl_matrix_view;

    // MatrixFloat functions
    pub fn gsl_matrix_float_alloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_float;
    pub fn gsl_matrix_float_calloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_float;
    pub fn gsl_matrix_float_free(m: *mut gsl_matrix_float);
    pub fn gsl_matrix_float_get(m: *const gsl_matrix_float, i: size_t, j: size_t) -> c_float;
    pub fn gsl_matrix_float_set(m: *mut gsl_matrix_float, i: size_t, j: size_t, x: c_float);
    pub fn gsl_matrix_float_set_all(m: *mut gsl_matrix_float, x: c_float);
    pub fn gsl_matrix_float_set_zero(m: *mut gsl_matrix_float);
    pub fn gsl_matrix_float_set_identity(m: *mut gsl_matrix_float);
    pub fn gsl_matrix_float_memcpy(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_swap(m: *mut gsl_matrix_float, w: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_get_row(vector: *mut gsl_vector_float, m: *const gsl_matrix_float, i: size_t) -> enums::Value;
    pub fn gsl_matrix_float_get_col(vector: *mut gsl_vector_float, m: *const gsl_matrix_float, j: size_t) -> enums::Value;
    pub fn gsl_matrix_float_set_row(m: *mut gsl_matrix_float, i: size_t, v: *const gsl_vector_float) -> enums::Value;
    pub fn gsl_matrix_float_set_col(m: *mut gsl_matrix_float, j: size_t, v: *const gsl_vector_float) -> enums::Value;
    pub fn gsl_matrix_float_swap_rows(m: *mut gsl_matrix_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_float_swap_columns(m: *mut gsl_matrix_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_float_swap_rowcol(m: *mut gsl_matrix_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_float_transpose_memcpy(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_transpose(m: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_add(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_sub(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_mul_elements(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_div_elements(dest: *mut gsl_matrix_float, src: *const gsl_matrix_float) -> enums::Value;
    pub fn gsl_matrix_float_scale(dest: *mut gsl_matrix_float, x: c_float) -> enums::Value;
    pub fn gsl_matrix_float_add_constant(dest: *mut gsl_matrix_float, x: c_float) -> enums::Value;
    pub fn gsl_matrix_float_max(m: *const gsl_matrix_float) -> c_float;
    pub fn gsl_matrix_float_min(m: *const gsl_matrix_float) -> c_float;
    pub fn gsl_matrix_float_minmax(m: *const gsl_matrix_float, min_out: *mut c_float, max_out: *mut c_float);
    pub fn gsl_matrix_float_max_index(m: *const gsl_matrix_float, imax: *mut size_t, jmax: *mut size_t);
    pub fn gsl_matrix_float_min_index(m: *const gsl_matrix_float, imin: *mut size_t, jmin: *mut size_t);
    pub fn gsl_matrix_float_minmax_index(m: *const gsl_matrix_float, imin: *mut size_t, jmin: *mut size_t, imax: *mut size_t, jmax: *mut size_t);
    pub fn gsl_matrix_float_isnull(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_ispos(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_isneg(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_isnonneg(m: *const gsl_matrix_float) -> c_int;
    pub fn gsl_matrix_float_equal(u: *const gsl_matrix_float, v: *const gsl_matrix_float) -> c_int;

    // MatrixComplex functions
    pub fn gsl_matrix_complex_alloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_complex;
    pub fn gsl_matrix_complex_calloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_complex;
    pub fn gsl_matrix_complex_free(m: *mut gsl_matrix_complex);
    pub fn gsl_matrix_complex_get(m: *const gsl_matrix_complex, i: size_t, j: size_t) -> gsl_complex;
    pub fn gsl_matrix_complex_set(m: *mut gsl_matrix_complex, i: size_t, j: size_t, x: gsl_complex);
    pub fn gsl_matrix_complex_set_all(m: *mut gsl_matrix_complex, x: gsl_complex);
    pub fn gsl_matrix_complex_set_zero(m: *mut gsl_matrix_complex);
    pub fn gsl_matrix_complex_set_identity(m: *mut gsl_matrix_complex);
    pub fn gsl_matrix_complex_memcpy(dest: *mut gsl_matrix_complex, src: *const gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_swap(m: *mut gsl_matrix_complex, w: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_get_row(vector: *mut gsl_vector_complex, m: *const gsl_matrix_complex, i: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_get_col(vector: *mut gsl_vector_complex, m: *const gsl_matrix_complex, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_set_row(m: *mut gsl_matrix_complex, i: size_t, v: *const gsl_vector_complex) -> enums::Value;
    pub fn gsl_matrix_complex_set_col(m: *mut gsl_matrix_complex, j: size_t, v: *const gsl_vector_complex) -> enums::Value;
    pub fn gsl_matrix_complex_swap_rows(m: *mut gsl_matrix_complex, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_swap_columns(m: *mut gsl_matrix_complex, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_swap_rowcol(m: *mut gsl_matrix_complex, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_transpose_memcpy(dest: *mut gsl_matrix_complex, src: *const gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_transpose(m: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_add(dest: *mut gsl_matrix_complex, src: *const gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_sub(dest: *mut gsl_matrix_complex, src: *const gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_mul_elements(dest: *mut gsl_matrix_complex, src: *const gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_div_elements(dest: *mut gsl_matrix_complex, src: *const gsl_matrix_complex) -> enums::Value;
    pub fn gsl_matrix_complex_scale(dest: *mut gsl_matrix_complex, x: gsl_complex) -> enums::Value;
    pub fn gsl_matrix_complex_add_constant(dest: *mut gsl_matrix_complex, x: gsl_complex) -> enums::Value;
    pub fn gsl_matrix_complex_isnull(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_ispos(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_isneg(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_isnonneg(m: *const gsl_matrix_complex) -> c_int;
    pub fn gsl_matrix_complex_equal(u: *const gsl_matrix_complex, v: *const gsl_matrix_complex) -> c_int;

    // MatrixComplexFloat functions
    pub fn gsl_matrix_complex_float_alloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_complex_float;
    pub fn gsl_matrix_complex_float_calloc(size1: size_t, size2: size_t) -> *mut gsl_matrix_complex_float;
    pub fn gsl_matrix_complex_float_free(m: *mut gsl_matrix_complex_float);
    pub fn gsl_matrix_complex_float_get(m: *const gsl_matrix_complex_float, i: size_t, j: size_t) -> gsl_complex_float;
    pub fn gsl_matrix_complex_float_set(m: *mut gsl_matrix_complex_float, i: size_t, j: size_t, x: gsl_complex_float);
    pub fn gsl_matrix_complex_float_set_all(m: *mut gsl_matrix_complex_float, x: gsl_complex_float);
    pub fn gsl_matrix_complex_float_set_zero(m: *mut gsl_matrix_complex_float);
    pub fn gsl_matrix_complex_float_set_identity(m: *mut gsl_matrix_complex_float);
    pub fn gsl_matrix_complex_float_memcpy(dest: *mut gsl_matrix_complex_float, src: *const gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_swap(m: *mut gsl_matrix_complex_float, w: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_get_row(vector: *mut gsl_vector_complex_float, m: *const gsl_matrix_complex_float, i: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_float_get_col(vector: *mut gsl_vector_complex_float, m: *const gsl_matrix_complex_float, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_float_set_row(m: *mut gsl_matrix_complex_float, i: size_t, v: *const gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_set_col(m: *mut gsl_matrix_complex_float, j: size_t, v: *const gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_swap_rows(m: *mut gsl_matrix_complex_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_float_swap_columns(m: *mut gsl_matrix_complex_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_float_swap_rowcol(m: *mut gsl_matrix_complex_float, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_matrix_complex_float_transpose_memcpy(dest: *mut gsl_matrix_complex_float, src: *const gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_transpose(m: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_add(dest: *mut gsl_matrix_complex_float, src: *const gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_sub(dest: *mut gsl_matrix_complex_float, src: *const gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_mul_elements(dest: *mut gsl_matrix_complex_float, src: *const gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_div_elements(dest: *mut gsl_matrix_complex_float, src: *const gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_scale(dest: *mut gsl_matrix_complex_float, x: gsl_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_add_constant(dest: *mut gsl_matrix_complex_float, x: gsl_complex_float) -> enums::Value;
    pub fn gsl_matrix_complex_float_isnull(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_ispos(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_isneg(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_isnonneg(m: *const gsl_matrix_complex_float) -> c_int;
    pub fn gsl_matrix_complex_float_equal(u: *const gsl_matrix_complex_float, v: *const gsl_matrix_complex_float) -> c_int;

    // Mathieu functions
    // Mathieu functions Workspace
    pub fn gsl_sf_mathieu_alloc(n: size_t, qmax: c_double) -> *mut gsl_sf_mathieu_workspace;
    pub fn gsl_sf_mathieu_free(work: *mut gsl_sf_mathieu_workspace);
    // Mathieu functions Characteristic Values
    pub fn gsl_sf_mathieu_a(n: c_int, q: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_mathieu_b(n: c_int, q: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_mathieu_a_array(order_min: c_int, order_max: c_int, q: c_double, work: *mut gsl_sf_mathieu_workspace,
        result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_mathieu_b_array(order_min: c_int, order_max: c_int, q: c_double, work: *mut gsl_sf_mathieu_workspace,
        result_array: *mut c_double) -> enums::Value;
    // Angular Mathieu functions
    pub fn gsl_sf_mathieu_ce(n: c_int, q: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_mathieu_se(n: c_int, q: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_mathieu_ce_array(nmin: c_int, nmax: c_int, q: c_double, x: c_double, work: *mut gsl_sf_mathieu_workspace,
        result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_mathieu_se_array(nmin: c_int, nmax: c_int, q: c_double, x: c_double, work: *mut gsl_sf_mathieu_workspace,
        result_array: *mut c_double) -> enums::Value;
    // Radial Mathieu functions
    pub fn gsl_sf_mathieu_Mc(j: c_int, n: c_int, q: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_mathieu_Ms(j: c_int, n: c_int, q: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_mathieu_Mc_array(j: c_int, nmin: c_int, nmax: c_int, q: c_double, x: c_double, work: *mut gsl_sf_mathieu_workspace,
        result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_mathieu_Ms_array(j: c_int, nmin: c_int, nmax: c_int, q: c_double, x: c_double, work: *mut gsl_sf_mathieu_workspace,
        result_array: *mut c_double) -> enums::Value;

    // Complex number functions
    // https://www.gnu.org/software/gsl/manual/html_node/Representation-of-complex-numbers.html#Representation-of-complex-numbers
    pub fn gsl_complex_rect(x: c_double, y: c_double) -> gsl_complex;
    pub fn gsl_complex_polar(r: c_double, theta: c_double) -> gsl_complex;
    pub fn gsl_complex_arg(z: gsl_complex) -> c_double;
    pub fn gsl_complex_abs(z: gsl_complex) -> c_double;
    pub fn gsl_complex_abs2(z: gsl_complex) -> c_double;
    pub fn gsl_complex_logabs(z: gsl_complex) -> c_double;
    pub fn gsl_complex_add(a: gsl_complex, b: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_sub(a: gsl_complex, b: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_mul(a: gsl_complex, b: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_div(a: gsl_complex, b: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_add_real(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_sub_real(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_mul_real(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_div_real(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_add_imag(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_sub_imag(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_mul_imag(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_div_imag(a: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_conjugate(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_inverse(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_negative(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_sqrt(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_sqrt_real(x: c_double) -> gsl_complex;
    pub fn gsl_complex_pow(z: gsl_complex, a: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_pow_real(z: gsl_complex, x: c_double) -> gsl_complex;
    pub fn gsl_complex_exp(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_log(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_log10(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_log_b(z: gsl_complex, b: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_sin(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_cos(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_tan(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_sec(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_csc(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_cot(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arcsin(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arcsin_real(z: c_double) -> gsl_complex;
    pub fn gsl_complex_arccos(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arccos_real(z: c_double) -> gsl_complex;
    pub fn gsl_complex_arctan(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arcsec(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arcsec_real(z: c_double) -> gsl_complex;
    pub fn gsl_complex_arccsc(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arccsc_real(z: c_double) -> gsl_complex;
    pub fn gsl_complex_arccot(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_sinh(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_cosh(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_tanh(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_sech(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_csch(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_coth(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arcsinh(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arccosh(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arccosh_real(z: c_double) -> gsl_complex;
    pub fn gsl_complex_arctanh(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arctanh_real(z: c_double) -> gsl_complex;
    pub fn gsl_complex_arcsech(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arccsch(z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_arccoth(z: gsl_complex) -> gsl_complex;

    // ComplexFloat number functions
    /*pub fn gsl_complex_float_arg(z: gsl_complex_float) -> c_float;
    pub fn gsl_complex_float_abs(z: gsl_complex_float) -> c_float;
    pub fn gsl_complex_float_abs2(z: gsl_complex_float) -> c_float;
    pub fn gsl_complex_float_logabs(z: gsl_complex_float) -> c_float;
    pub fn gsl_complex_float_add(a: gsl_complex_float, b: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sub(a: gsl_complex_float, b: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_mul(a: gsl_complex_float, b: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_div(a: gsl_complex_float, b: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_add_real(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sub_real(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_mul_real(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_div_real(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_add_imag(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sub_imag(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_mul_imag(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_div_imag(a: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_conjugate(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_inverse(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_negative(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sqrt(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sqrt_real(x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_pow(z: gsl_complex_float, a: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_pow_real(z: gsl_complex_float, x: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_exp(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_log(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_log10(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_log_b(z: gsl_complex_float, b: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sin(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_cos(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_tan(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sec(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_csc(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_cot(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arcsin(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arcsin_real(z: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccos(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccos_real(z: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arctan(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arcsec(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arcsec_real(z: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccsc(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccsc_real(z: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccot(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sinh(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_cosh(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_tanh(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_sech(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_csch(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_coth(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arcsinh(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccosh(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccosh_real(z: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arctanh(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arctanh_real(z: c_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arcsech(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccsch(z: gsl_complex_float) -> gsl_complex_float;
    pub fn gsl_complex_float_arccoth(z: gsl_complex_float) -> gsl_complex_float;*/

    // Basis Splines
    pub fn gsl_bspline_alloc(k: size_t, nbreak: size_t) -> *mut gsl_bspline_workspace;
    pub fn gsl_bspline_free(w: *mut gsl_bspline_workspace);
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_alloc(k: size_t) -> *mut gsl_bspline_deriv_workspace;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_free(w: *mut gsl_bspline_deriv_workspace);
    pub fn gsl_bspline_knots(breakpts: *mut gsl_vector, w: *mut gsl_bspline_workspace) -> enums::Value;
    pub fn gsl_bspline_knots_uniform(a: c_double, b: c_double, w: *mut gsl_bspline_workspace) -> enums::Value;
    pub fn gsl_bspline_eval(x: c_double, B: *mut gsl_vector, w: *mut gsl_bspline_workspace) -> enums::Value;
    pub fn gsl_bspline_eval_nonzero(x: c_double, Bk: *mut gsl_vector, istart: *mut size_t, iend: *mut size_t,
        w: *mut gsl_bspline_workspace) -> enums::Value;
    pub fn gsl_bspline_ncoeffs(w: *mut gsl_bspline_workspace) -> size_t;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_eval(x: c_double, nderiv: size_t, dB: *mut gsl_matrix, w: *mut gsl_bspline_workspace,
        dw: *mut gsl_bspline_deriv_workspace) -> enums::Value;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_eval_nonzero(x: c_double, nderiv: size_t, Bk: *mut gsl_matrix, istart: *mut size_t, iend: *mut size_t,
        w: *mut gsl_bspline_workspace, dw: *mut gsl_bspline_deriv_workspace) -> enums::Value;
    pub fn gsl_bspline_greville_abscissa(i: size_t, w: *mut gsl_bspline_workspace) -> c_double;

    // Level 1 BLAS functions
    pub fn gsl_blas_sdsdot(alpha: c_float, x: *const gsl_vector_float, y: *const gsl_vector_float, result: *mut c_float) -> enums::Value;
    pub fn gsl_blas_sdot(x: *const gsl_vector_float, y: *const gsl_vector_float, result: *mut c_float) -> enums::Value;
    pub fn gsl_blas_dsdot(x: *const gsl_vector_float, y: *const gsl_vector_float, result: *mut c_double) -> enums::Value;
    pub fn gsl_blas_ddot(x: *const gsl_vector, y: *const gsl_vector, result: *mut c_double) -> enums::Value;
    pub fn gsl_blas_cdotu(x: *const gsl_vector_complex_float, y: *const gsl_vector_complex_float, dotu: *mut gsl_complex_float) -> enums::Value;
    pub fn gsl_blas_zdotu(x: *const gsl_vector_complex, y: *const gsl_vector_complex, dotu: *mut gsl_complex) -> enums::Value;
    pub fn gsl_blas_cdotc(x: *const gsl_vector_complex_float, y: *const gsl_vector_complex_float, dotc: *mut gsl_complex_float) -> enums::Value;
    pub fn gsl_blas_zdotc(x: *const gsl_vector_complex, y: *const gsl_vector_complex, dotc: *mut gsl_complex) -> enums::Value;
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
    pub fn gsl_blas_sswap(x: *mut gsl_vector_float, y: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_blas_dswap(x: *mut gsl_vector, y: *mut gsl_vector) -> enums::Value;
    pub fn gsl_blas_cswap(x: *mut gsl_vector_complex_float, y: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_blas_zswap(x: *mut gsl_vector_complex, y: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_blas_scopy(x: *const gsl_vector_float, y: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_blas_dcopy(x: *const gsl_vector, y: *mut gsl_vector) -> enums::Value;
    pub fn gsl_blas_ccopy(x: *const gsl_vector_complex_float, y: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_blas_zcopy(x: *const gsl_vector_complex, y: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_blas_saxpy(alpha: c_float, x: *const gsl_vector_float, y: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_blas_daxpy(alpha: f64, x: *const gsl_vector, y: *mut gsl_vector) -> enums::Value;
    pub fn gsl_blas_caxpy(alpha: gsl_complex_float, x: *const gsl_vector_complex_float, y: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_blas_zaxpy(alpha: gsl_complex, x: *const gsl_vector_complex, y: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_blas_sscal(alpha: c_float, x: *mut gsl_vector_float);
    pub fn gsl_blas_dscal(alpha: c_double, x: *mut gsl_vector);
    pub fn gsl_blas_cscal(alpha: gsl_complex_float, x: *mut gsl_vector_complex_float);
    pub fn gsl_blas_zscal(alpha: gsl_complex, x: *mut gsl_vector_complex);
    pub fn gsl_blas_csscal(alpha: c_float, x: *mut gsl_vector_complex_float);
    pub fn gsl_blas_zdscal(alpha: c_double, x: *mut gsl_vector_complex);
    pub fn gsl_blas_srotg(a: *mut c_float, b: *mut c_float, c: *mut c_float, d: *mut c_float) -> enums::Value;
    pub fn gsl_blas_drotg(a: *mut c_double, b: *mut c_double, c: *mut c_double, d: *mut c_double) -> enums::Value;
    pub fn gsl_blas_srot(a: *mut gsl_vector_float, b: *mut gsl_vector_float, c: c_float, d: c_float) -> enums::Value;
    pub fn gsl_blas_drot(a: *mut gsl_vector, b: *mut gsl_vector, c: c_double, d: c_double) -> enums::Value;
    pub fn gsl_blas_srotmg(d1: *mut c_float, d2: *mut c_float, b1: *mut c_float, b2: c_float, P: *mut c_float) -> enums::Value;
    pub fn gsl_blas_drotmg(d1: *mut c_double, d2: *mut c_double, b1: *mut c_double, b2: c_double, P: *mut c_double) -> enums::Value;
    pub fn gsl_blas_srotm(x: *mut gsl_vector_float, y: *mut gsl_vector_float, P: *mut c_float) -> enums::Value;
    pub fn gsl_blas_drotm(x: *mut gsl_vector, y: *mut gsl_vector, P: *mut c_double) -> enums::Value;
    // Level 2 BLAS functions
    pub fn gsl_blas_sgemv(transA: CBLAS_TRANSPOSE_t, alpha: c_float, A: *const gsl_matrix_float, x: *const gsl_vector_float, beta: c_float,
        y: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_blas_dgemv(transA: CBLAS_TRANSPOSE_t, alpha: c_double, A: *const gsl_matrix, x: *const gsl_vector, beta: c_double,
        y: *mut gsl_vector) -> enums::Value;
    pub fn gsl_blas_cgemv(transA: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float,
        x: *const gsl_vector_complex_float, beta: gsl_complex_float, y: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_blas_zgemv(transA: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: *const gsl_matrix_complex, x: *const gsl_vector_complex,
        beta: gsl_complex, y: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_blas_strmv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix_float,
        x: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_blas_dtrmv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix,
        x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_blas_ctrmv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix_complex_float,
        x: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_blas_ztrmv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix_complex,
        x: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_blas_strsv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix_float,
        x: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_blas_dtrsv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix,
        x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_blas_ctrsv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix_complex_float,
        x: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_blas_ztrsv(uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, A: *const gsl_matrix_complex,
        x: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_blas_ssymv(uplo: CBLAS_UPLO_t, alpha: c_float, A: *const gsl_matrix_float, x: *const gsl_vector_float, beta: c_float,
        y: *mut gsl_vector_float) -> enums::Value;
    pub fn gsl_blas_dsymv(uplo: CBLAS_UPLO_t, alpha: c_double, A: *const gsl_matrix, x: *const gsl_vector, beta: c_double,
        y: *mut gsl_vector) -> enums::Value;
    pub fn gsl_blas_chemv(uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float, x: *const gsl_vector_complex_float,
        beta: gsl_complex_float, y: *mut gsl_vector_complex_float) -> enums::Value;
    pub fn gsl_blas_zhemv(uplo: CBLAS_UPLO_t, alpha: gsl_complex, A: *const gsl_matrix_complex, x: *const gsl_vector_complex,
        beta: gsl_complex, y: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_blas_sger(alpha: c_float, x: *const gsl_vector_float, y: *const gsl_vector_float, A: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dger(alpha: c_double, x: *const gsl_vector, y: *const gsl_vector, A: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_cgeru(alpha: gsl_complex_float, x: *const gsl_vector_complex_float, y: *const gsl_vector_complex_float,
        A: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zgeru(alpha: gsl_complex, x: *const gsl_vector_complex, y: *const gsl_vector_complex, A: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_cgerc(alpha: gsl_complex_float, x: *const gsl_vector_complex_float, y: *const gsl_vector_complex_float,
        A: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zgerc(alpha: gsl_complex, x: *const gsl_vector_complex, y: *const gsl_vector_complex, A: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_ssyr(uplo: CBLAS_UPLO_t, alpha: c_float, x: *const gsl_vector_float, A: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dsyr(uplo: CBLAS_UPLO_t, alpha: c_double, x: *const gsl_vector, A: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_cher(uplo: CBLAS_UPLO_t, alpha: c_float, x: *const gsl_vector_complex_float, A: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zher(uplo: CBLAS_UPLO_t, alpha: c_double, x: *const gsl_vector_complex, A: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_ssyr2(uplo: CBLAS_UPLO_t, alpha: c_float, x: *const gsl_vector_float, y: *const gsl_vector_float,
        A: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dsyr2(uplo: CBLAS_UPLO_t, alpha: c_double, x: *const gsl_vector, y: *const gsl_vector,
        A: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_cher2(uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, x: *const gsl_vector_complex_float, y: *const gsl_vector_complex_float,
        A: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zher2(uplo: CBLAS_UPLO_t, alpha: gsl_complex, x: *const gsl_vector_complex, y: *const gsl_vector_complex,
        A: *mut gsl_matrix_complex) -> enums::Value;
    // Level 3 BLAS functions
    pub fn gsl_blas_sgemm(transA: CBLAS_TRANSPOSE_t, transB: CBLAS_TRANSPOSE_t, alpha: c_float, A: *const gsl_matrix_float,
        B: *const gsl_matrix_float, beta: c_float, C: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dgemm(transA: CBLAS_TRANSPOSE_t, transB: CBLAS_TRANSPOSE_t, alpha: c_double, A: *const gsl_matrix,
        B: *const gsl_matrix, beta: c_double, C: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_cgemm(transA: CBLAS_TRANSPOSE_t, transB: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float, beta: gsl_complex_float, C: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zgemm(transA: CBLAS_TRANSPOSE_t, transB: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex, beta: gsl_complex, C: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_ssymm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, alpha: c_float, A: *const gsl_matrix_float, B: *const gsl_matrix_float,
        beta: c_float, C: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dsymm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, alpha: c_double, A: *const gsl_matrix, B: *const gsl_matrix,
        beta: c_double, C: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_csymm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float, beta: gsl_complex_float, C: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zsymm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, alpha: gsl_complex, A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex, beta: gsl_complex, C: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_chemm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float, beta: gsl_complex_float, C: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zhemm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, alpha: gsl_complex, A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex, beta: gsl_complex, C: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_strmm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: c_float,
        A: *const gsl_matrix_float, B: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dtrmm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: c_double,
        A: *const gsl_matrix, B: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_ctrmm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float, B: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_ztrmm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: gsl_complex,
        A: *const gsl_matrix_complex, B: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_strsm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: c_float,
        A: *const gsl_matrix_float, B: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dtrsm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: c_double,
        A: *const gsl_matrix, B: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_ctrsm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: gsl_complex_float,
        A: *const gsl_matrix_complex_float, B: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_ztrsm(side: CBLAS_SIDE_t, uplo: CBLAS_UPLO_t, transA: CBLAS_TRANSPOSE_t, diag: CBLAS_DIAG_t, alpha: gsl_complex,
        A: *const gsl_matrix_complex, B: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_ssyrk(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: c_float, A: *const gsl_matrix_float, beta: c_float,
        C: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dsyrk(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: c_double, A: *const gsl_matrix, beta: c_double,
        C: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_csyrk(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float,
        beta: gsl_complex_float, C: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zsyrk(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: *const gsl_matrix_complex,
        beta: gsl_complex, C: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_cherk(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: c_float, A: *const gsl_matrix_complex_float,
        beta: c_float, C: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zherk(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: c_double, A: *const gsl_matrix_complex,
        beta: c_double, C: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_ssyr2k(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: c_float, A: *const gsl_matrix_float,
        B: *const gsl_matrix_float, beta: c_float, C: *mut gsl_matrix_float) -> enums::Value;
    pub fn gsl_blas_dsyr2k(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: c_double, A: *const gsl_matrix, B: *const gsl_matrix,
        beta: c_double, C: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_blas_csyr2k(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float, beta: gsl_complex_float, C: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zsyr2k(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex, beta: gsl_complex, C: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_blas_cher2k(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex_float, A: *const gsl_matrix_complex_float,
        B: *const gsl_matrix_complex_float, beta: c_float, C: *mut gsl_matrix_complex_float) -> enums::Value;
    pub fn gsl_blas_zher2k(uplo: CBLAS_UPLO_t, trans: CBLAS_TRANSPOSE_t, alpha: gsl_complex, A: *const gsl_matrix_complex,
        B: *const gsl_matrix_complex, beta: c_double, C: *mut gsl_matrix_complex) -> enums::Value;

    // Fit functions
    pub fn gsl_fit_linear(x: *const c_double, xstride: size_t, y: *const c_double, ystride: size_t, n: size_t, c0: *mut c_double, c1: *mut c_double,
        cov00: *mut c_double, cov01: *mut c_double, cov11: *mut c_double, sumsq: c_double) -> enums::Value;
    pub fn gsl_fit_wlinear(x: *const c_double, xstride: size_t, w: *const c_double, wstride: size_t, y: *const c_double, ystride: size_t,
        n: size_t, c0: *mut c_double, c1: *mut c_double, cov00: *mut c_double, cov01: *mut c_double, cov11: *mut c_double,
        chisq: *mut c_double) -> enums::Value;
    pub fn gsl_fit_linear_est(x: c_double, c0: c_double, c1: c_double, cov00: c_double, cov01: c_double, cov11: c_double, y: *mut c_double,
        y_err: *mut c_double) -> enums::Value;
    pub fn gsl_fit_mul(x: *const c_double, xstride: size_t, y: *const c_double, ystride: size_t, n: size_t, c1: *mut c_double,
        cov11: *mut c_double, sumsq: *mut c_double) -> enums::Value;
    pub fn gsl_fit_wmul(x: *const c_double, xstride: size_t, w: *const c_double, wstride: size_t, y: *const c_double, ystride: size_t,
        n: size_t, c1: *mut c_double, cov11: *mut c_double, sumsq: *mut c_double) -> enums::Value;
    pub fn gsl_fit_mul_est(x: c_double, c1: c_double, cov11: c_double, y: *mut c_double, y_err: *mut c_double) -> enums::Value;

    // Pow functions
    pub fn gsl_pow_int(x: c_double, n: c_int) -> c_double;
    pub fn gsl_pow_uint(x: c_double, n: c_uint) -> c_double;
    pub fn gsl_pow_2(x: c_double) -> c_double;
    pub fn gsl_pow_3(x: c_double) -> c_double;
    pub fn gsl_pow_4(x: c_double) -> c_double;
    pub fn gsl_pow_5(x: c_double) -> c_double;
    pub fn gsl_pow_6(x: c_double) -> c_double;
    pub fn gsl_pow_7(x: c_double) -> c_double;
    pub fn gsl_pow_8(x: c_double) -> c_double;
    pub fn gsl_pow_9(x: c_double) -> c_double;

    // Random Number Generation
    pub fn gsl_rng_alloc(T: *const gsl_rng_type) -> *mut gsl_rng;
    pub fn gsl_rng_set(r: *const gsl_rng, s: c_ulong);
    pub fn gsl_rng_free(r: *mut gsl_rng);
    pub fn gsl_rng_get(r: *const gsl_rng) -> c_ulong;
    pub fn gsl_rng_uniform(r: *const gsl_rng) -> c_double;
    pub fn gsl_rng_uniform_pos(r: *const gsl_rng) -> c_double;
    pub fn gsl_rng_uniform_int(r: *const gsl_rng, n: c_ulong) -> c_ulong;
    pub fn gsl_rng_name(r: *const gsl_rng) -> *const c_char;
    pub fn gsl_rng_max(r: *const gsl_rng) -> c_ulong;
    pub fn gsl_rng_min(r: *const gsl_rng) -> c_ulong;
    pub fn gsl_rng_state(r: *const gsl_rng) -> *mut c_void;
    pub fn gsl_rng_size(r: *const gsl_rng) -> size_t;
    pub fn gsl_rng_types_setup() -> *const *mut gsl_rng_type;
    pub fn gsl_rng_memcpy(dest: *mut gsl_rng, src: *const gsl_rng) -> enums::Value;
    pub fn gsl_rng_clone(r: *const gsl_rng) -> *mut gsl_rng;
    pub fn gsl_rng_env_setup() -> *const gsl_rng_type;

    // Random Number Distributions
    // The Gaussian Distribution
    pub fn gsl_ran_gaussian(r: *const gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_pdf(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_ziggurat(r: *const gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_ratio_method(r: *const gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_ugaussian(r: *const gsl_rng) -> c_double;
    pub fn gsl_ran_ugaussian_pdf(x: c_double) -> c_double;
    pub fn gsl_ran_ugaussian_ratio_method(r: *const gsl_rng) -> c_double;
    pub fn gsl_cdf_gaussian_P(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_gaussian_Q(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_gaussian_Pinv(P: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_gaussian_Qinv(Q: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_P(x: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_Q(x: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_Pinv(P: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_Qinv(Q: c_double) -> c_double;
    // The Gaussian Tail Distribution
    pub fn gsl_ran_gaussian_tail(r: *const gsl_rng, a: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_tail_pdf(x: c_double, a: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_ugaussian_tail(r: *const gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_ugaussian_tail_pdf(x: c_double, a: c_double) -> c_double;
    // The Bivariate Gaussian Distribution
    pub fn gsl_ran_bivariate_gaussian(r: *const gsl_rng, sigma_x: c_double, sigma_y: c_double, rho: c_double, x: *mut c_double,
        y: *mut c_double);
    pub fn gsl_ran_bivariate_gaussian_pdf(x: c_double, y: c_double, sigma_x: c_double, sigma_y: c_double, rho: c_double) -> c_double;
    // The Exponential Distribution
    pub fn gsl_ran_exponential(r: *const gsl_rng, mu: c_double) -> c_double;
    pub fn gsl_ran_exponential_pdf(x: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_P(x: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_Q(x: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_Pinv(P: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_Qinv(Q: c_double, mu: c_double) -> c_double;
    // The Laplace Distribution
    pub fn gsl_ran_laplace(r: *const gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_laplace_pdf(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_P(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_Q(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_Pinv(P: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_Qinv(Q: c_double, a: c_double) -> c_double;
    // The Exponential Power Distribution
    pub fn gsl_ran_exppow(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_exppow_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_exppow_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_exppow_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    // The Cauchy Distribution
    pub fn gsl_ran_cauchy(r: *const gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_cauchy_pdf(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_P(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_Q(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_Pinv(P: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_Qinv(Q: c_double, a: c_double) -> c_double;
    // The Rayleigh Distribution
    pub fn gsl_ran_rayleigh(r: *const gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_rayleigh_pdf(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_P(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_Q(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_Pinv(P: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_Qinv(Q: c_double, sigma: c_double) -> c_double;
    // The Rayleigh Tail Distribution
    pub fn gsl_ran_rayleigh_tail(r: *const gsl_rng, a: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_rayleigh_tail_pdf(x: c_double, a: c_double, sigma: c_double) -> c_double;
    // The Landau Distribution
    pub fn gsl_ran_landau(r: *const gsl_rng) -> c_double;
    pub fn gsl_ran_landau_pdf(x: c_double) -> c_double;
    // The Levy alpha-Stable Distributions
    pub fn gsl_ran_levy(r: *const gsl_rng, c: c_double, alpha: c_double) -> c_double;
    // The Levy skew alpha-Stable Distribution
    pub fn gsl_ran_levy_skew(r: *const gsl_rng, c: c_double, alpha: c_double, beta: c_double) -> c_double;
    // The Gamma Distribution
    pub fn gsl_ran_gamma(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gamma_knuth(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gamma_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Flat (Uniform) Distribution
    pub fn gsl_ran_flat(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_flat_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Lognormal Distribution
    pub fn gsl_ran_lognormal(r: *const gsl_rng, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_lognormal_pdf(x: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_P(x: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_Q(x: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_Pinv(P: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_Qinv(Q: c_double, zeta: c_double, sigma: c_double) -> c_double;
    // The Chi-squared Distribution
    pub fn gsl_ran_chisq(r: *const gsl_rng, nu: c_double) -> c_double;
    pub fn gsl_ran_chisq_pdf(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_P(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_Q(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_Pinv(P: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_Qinv(Q: c_double, nu: c_double) -> c_double;
    // The F-distribution
    pub fn gsl_ran_fdist(r: *const gsl_rng, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_ran_fdist_pdf(x: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_P(x: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_Q(x: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_Pinv(P: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_Qinv(Q: c_double, nu1: c_double, nu2: c_double) -> c_double;
    // The t-distribution
    pub fn gsl_ran_tdist(r: *const gsl_rng, nu: c_double) -> c_double;
    pub fn gsl_ran_tdist_pdf(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_P(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_Q(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_Pinv(P: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_Qinv(Q: c_double, nu: c_double) -> c_double;
    // The Beta Distribution
    pub fn gsl_ran_beta(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_beta_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Logistic Distribution
    pub fn gsl_ran_logistic(r: *const gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_logistic_pdf(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_P(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_Q(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_Pinv(P: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_Qinv(Q: c_double, a: c_double) -> c_double;
    // The Pareto Distribution
    pub fn gsl_ran_pareto(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_pareto_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // Spherical Vector Distributions
    pub fn gsl_ran_dir_2d(r: *const gsl_rng, x: *mut c_double, y: *mut c_double);
    pub fn gsl_ran_dir_2d_trig_method(r: *const gsl_rng, x: *mut c_double, y: *mut c_double);
    pub fn gsl_ran_dir_3d(r: *const gsl_rng, x: *mut c_double, y: *mut c_double, z: *mut c_double);
    pub fn gsl_ran_dir_nd(r: *const gsl_rng, n: size_t, x: *mut c_double);
    // The Weibull Distribution
    pub fn gsl_ran_weibull(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_weibull_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Type-1 Gumbel Distribution
    pub fn gsl_ran_gumbel1(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gumbel1_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Type-2 Gumbel Distribution
    pub fn gsl_ran_gumbel2(r: *const gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gumbel2_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Dirichlet Distribution
    pub fn gsl_ran_dirichlet(r: *const gsl_rng, K: size_t, alpha: *const c_double, theta: *mut c_double);
    pub fn gsl_ran_dirichlet_pdf(K: size_t, alpha: *const c_double, theta: *const c_double) -> c_double;
    pub fn gsl_ran_dirichlet_lnpdf(K: size_t, alpha: *const c_double, theta: *const c_double) -> c_double;
    // General Discrete Distributions
    pub fn gsl_ran_discrete_preproc(K: size_t, P: *const c_double) -> *mut gsl_ran_discrete_t;
    pub fn gsl_ran_discrete(r: *const gsl_rng, g: *const gsl_ran_discrete_t) -> size_t;
    pub fn gsl_ran_discrete_pdf(k: size_t, g: *const gsl_ran_discrete_t) -> c_double;
    pub fn gsl_ran_discrete_free(g: *mut gsl_ran_discrete_t);
    // The Poisson Distribution
    pub fn gsl_ran_poisson(r: *const gsl_rng, mu: c_double) -> c_uint;
    pub fn gsl_ran_poisson_pdf(k: c_uint, mu: c_double) -> c_double;
    pub fn gsl_cdf_poisson_P(k: c_uint, mu: c_double) -> c_double;
    pub fn gsl_cdf_poisson_Q(k: c_uint, mu: c_double) -> c_double;
    // The Bernoulli Distribution
    pub fn gsl_ran_bernoulli(r: *const gsl_rng, p: c_double) -> c_uint;
    pub fn gsl_ran_bernoulli_pdf(k: c_uint, p: c_double) -> c_double;
    // The Binomial Distribution
    pub fn gsl_ran_binomial(r: *const gsl_rng, p: c_double, n: c_uint) -> c_uint;
    pub fn gsl_ran_binomial_pdf(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_binomial_P(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_binomial_Q(k: c_uint, p: c_double, n: c_uint) -> c_double;
    // The Multinomial Distribution
    pub fn gsl_ran_multinomial(r: *const gsl_rng, K: size_t, N: c_uint, p: *const c_double, n: *mut c_uint);
    pub fn gsl_ran_multinomial_pdf(K: size_t, p: *const c_double, n: *const c_uint) -> c_double;
    pub fn gsl_ran_multinomial_lnpdf(K: size_t, p: *const c_double, n: *const c_uint) -> c_double;
    // The Negative Binomial Distribution
    pub fn gsl_ran_negative_binomial(r: *const gsl_rng, p: c_double, n: c_double) -> c_uint;
    pub fn gsl_ran_negative_binomial_pdf(k: c_uint, p: c_double, n: c_double) -> c_double;
    pub fn gsl_cdf_negative_binomial_P(k: c_uint, p: c_double, n: c_double) -> c_double;
    pub fn gsl_cdf_negative_binomial_Q(k: c_uint, p: c_double, n: c_double) -> c_double;
    // The Pascal Distribution
    pub fn gsl_ran_pascal(r: *const gsl_rng, p: c_double, n: c_uint) -> c_uint;
    pub fn gsl_ran_pascal_pdf(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_pascal_P(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_pascal_Q(k: c_uint, p: c_double, n: c_uint) -> c_double;
    // The Geometric Distribution
    pub fn gsl_ran_geometric(r: *const gsl_rng, p: c_double) -> c_uint;
    pub fn gsl_ran_geometric_pdf(k: c_uint, p: c_double) -> c_double;
    pub fn gsl_cdf_geometric_P(k: c_uint, p: c_double) -> c_double;
    pub fn gsl_cdf_geometric_Q(k: c_uint, p: c_double) -> c_double;
    // The Hypergeometric Distribution
    pub fn gsl_ran_hypergeometric(r: *const gsl_rng, n1: c_uint, n2: c_uint, t: c_uint) -> c_uint;
    pub fn gsl_ran_hypergeometric_pdf(k: c_uint, n1: c_uint, n2: c_uint, t: c_uint) -> c_double;
    pub fn gsl_cdf_hypergeometric_P(k: c_uint, n1: c_uint, n2: c_uint, t: c_uint) -> c_double;
    pub fn gsl_cdf_hypergeometric_Q(k: c_uint, n1: c_uint, n2: c_uint, t: c_uint) -> c_double;
    // The Logarithmic Distribution
    pub fn gsl_ran_logarithmic(r: *const gsl_rng, p: c_double) -> c_uint;
    pub fn gsl_ran_logarithmic_pdf(k: c_uint, p: c_double) -> c_double;
    // Shuffling and Sampling
    pub fn gsl_ran_shuffle(r: *const gsl_rng, base: *mut c_void, n: size_t, size: size_t);
    pub fn gsl_ran_choose(r: *const gsl_rng, dest: *mut c_void, k: size_t, src: *mut c_void, n: size_t, size: size_t) -> enums::Value;
    pub fn gsl_ran_sample(r: *const gsl_rng, dest: *mut c_void, k: size_t, src: *mut c_void, n: size_t, size: size_t) -> enums::Value;

    // Permutation struct
    pub fn gsl_permutation_alloc(size: size_t) -> *mut gsl_permutation;
    pub fn gsl_permutation_calloc(size: size_t) -> *mut gsl_permutation;
    pub fn gsl_permutation_init(p: *mut gsl_permutation);
    pub fn gsl_permutation_free(p: *mut gsl_permutation);
    pub fn gsl_permutation_memcpy(dest: *mut gsl_permutation, src: *const gsl_permutation) -> enums::Value;
    pub fn gsl_permutation_get(p: *const gsl_permutation, i: size_t) -> size_t;
    pub fn gsl_permutation_swap(p: *mut gsl_permutation, i: size_t, j: size_t) -> enums::Value;
    pub fn gsl_permutation_size(p: *const gsl_permutation) -> size_t;
    //pub fn gsl_permutation_data(p: *const gsl_permutation) -> *mut size_t;
    pub fn gsl_permutation_valid(p: *const gsl_permutation) -> enums::Value;
    pub fn gsl_permutation_reverse(p: *mut gsl_permutation);
    pub fn gsl_permutation_inverse(inv: *mut gsl_permutation, p: *const gsl_permutation) -> enums::Value;
    pub fn gsl_permutation_next(p: *mut gsl_permutation) -> enums::Value;
    pub fn gsl_permutation_prev(p: *mut gsl_permutation) -> enums::Value;
    pub fn gsl_permute(p: *const size_t, data: *mut c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_permute_inverse(p: *const size_t, data: *mut c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_permute_vector(p: *const gsl_permutation, v: *mut gsl_vector) -> enums::Value;
    pub fn gsl_permute_vector_inverse(p: *const gsl_permutation, v: *mut gsl_vector) -> enums::Value;
    pub fn gsl_permutation_mul(p: *mut gsl_permutation, pa: *const gsl_permutation, pb: *const gsl_permutation) -> enums::Value;
    pub fn gsl_permutation_linear_to_canonical(q: *mut gsl_permutation, p: *const gsl_permutation) -> enums::Value;
    pub fn gsl_permutation_canonical_to_linear(p: *mut gsl_permutation, q: *const gsl_permutation) -> enums::Value;
    pub fn gsl_permutation_inversions(p: *const gsl_permutation) -> size_t;
    pub fn gsl_permutation_linear_cycles(p: *const gsl_permutation) -> size_t;
    pub fn gsl_permutation_canonical_cycles(p: *const gsl_permutation) -> size_t;

    // Sorting functions
    // Sorting objects
    //pub fn gsl_heapsort(array: *mut c_void, count: size_t, size: size_t, compare: compare_fn);
    //pub fn gsl_heapsort_index(p: *mut size_t, array: *const c_void, count: size_t, size: size_t, compare: compare_fn) -> enums::Value;
    // Sorting vectors
    pub fn gsl_sort(data: *mut c_double, stride: size_t, n: size_t);
    pub fn gsl_sort2(data1: *mut c_double, stride1: size_t, data2: *mut c_double, stride2: size_t, n: size_t);
    pub fn gsl_sort_vector(v: *mut gsl_vector);
    pub fn gsl_sort_vector2(v1: *mut gsl_vector, v2: *mut gsl_vector);
    pub fn gsl_sort_index(p: *mut size_t, data: *const c_double, stride: size_t, n: size_t);
    pub fn gsl_sort_vector_index(p: *mut gsl_permutation, v: *const gsl_vector) -> enums::Value;
    // Selecting the k smallest or largest elements
    pub fn gsl_sort_smallest(dest: *mut c_double, k: size_t, src: *const c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_sort_largest(dest: *mut c_double, k: size_t, src: *const c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_sort_vector_smallest(dest: *mut c_double, k: size_t, v: *const gsl_vector) -> enums::Value;
    pub fn gsl_sort_vector_largest(dest: *mut c_double, k: size_t, v: *const gsl_vector) -> enums::Value;
    pub fn gsl_sort_smallest_index(p: *mut size_t, k: size_t, src: *const c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_sort_largest_index(p: *mut size_t, k: size_t, src: *const c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_sort_vector_smallest_index(p: *mut size_t, k: size_t, v: *const gsl_vector) -> enums::Value;
    pub fn gsl_sort_vector_largest_index(p: *mut size_t, k: size_t, v: *const gsl_vector) -> enums::Value;

    // Chebyshev Approximations
    // Creation and Calculation of Chebyshev Series
    pub fn gsl_cheb_alloc(n: size_t) -> *mut gsl_cheb_series;
    pub fn gsl_cheb_free(cs: *mut gsl_cheb_series);
    // Auxiliary functions
    pub fn gsl_cheb_order(cs: *const gsl_cheb_series) -> size_t;
    pub fn gsl_cheb_size(cs: *const gsl_cheb_series) -> size_t;
    // Chebyshev Series Evaluation
    pub fn gsl_cheb_eval(cs: *const gsl_cheb_series, x: c_double) -> c_double;
    pub fn gsl_cheb_eval_err(cs: *const gsl_cheb_series, x: c_double, result: *mut c_double, abs_err: *mut c_double) -> enums::Value;
    pub fn gsl_cheb_eval_n(cs: *const gsl_cheb_series, order: size_t, x: c_double) -> c_double;
    pub fn gsl_cheb_eval_n_err(cs: *const gsl_cheb_series, order: size_t, x: c_double, result: *mut c_double,
        abs_err: *mut c_double) -> enums::Value;
    // Derivatives and Integrals
    pub fn gsl_cheb_calc_deriv(cs: *mut gsl_cheb_series, deriv: *const gsl_cheb_series) -> enums::Value;
    pub fn gsl_cheb_calc_integ(cs: *mut gsl_cheb_series, integ: *const gsl_cheb_series) -> enums::Value;

    // Error function
    #[allow(dead_code)]
    pub fn gsl_error(reason: *const c_char, file: *const c_char, line: c_int, gsl_errno: c_int);

    // Combination
    // Combination allocation
    pub fn gsl_combination_alloc(n: size_t, k: size_t) -> *mut gsl_combination;
    pub fn gsl_combination_calloc(n: size_t, k: size_t) -> *mut gsl_combination;
    pub fn gsl_combination_init_first(c: *mut gsl_combination);
    pub fn gsl_combination_init_last(c: *mut gsl_combination);
    pub fn gsl_combination_free(c: *mut gsl_combination);
    pub fn gsl_combination_memcpy(dest: *mut gsl_combination, src: *const gsl_combination) -> enums::Value;
    // Accessing combination elements
    pub fn gsl_combination_get(c: *const gsl_combination, i: size_t) -> size_t;
    // Combination properties
    pub fn gsl_combination_n(c: *const gsl_combination) -> size_t;
    pub fn gsl_combination_k(c: *const gsl_combination) -> size_t;
    pub fn gsl_combination_valid(c: *mut gsl_combination) -> enums::Value;
    // Combination functions
    pub fn gsl_combination_next(c: *mut gsl_combination) -> enums::Value;
    pub fn gsl_combination_prev(c: *mut gsl_combination) -> enums::Value;

    // Polynomials
    // Polynomial Evaluation
    pub fn gsl_poly_eval(c: *const c_double, len: c_int, x: c_double) -> c_double;
    pub fn gsl_poly_complex_eval(c: *const c_double, len: c_int, z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_poly_complex_eval(c: *const gsl_complex, len: c_int, z: gsl_complex) -> gsl_complex;
    pub fn gsl_poly_eval_derivs(c: *const c_double, lenc: size_t, x: c_double, res: *mut c_double, lenres: size_t) -> enums::Value;
    // Divided Difference Representation of Polynomials
    pub fn gsl_poly_dd_init(dd: *mut c_double, xa: *const c_double, ya: *const c_double, size: size_t) -> enums::Value;
    pub fn gsl_poly_dd_eval(dd: *const c_double, xa: *const c_double, size: size_t, x: c_double) -> c_double;
    pub fn gsl_poly_dd_taylor(c: *mut c_double, xp: c_double, dd: *const c_double, xa: *const c_double, size: size_t,
        w: *mut c_double) -> enums::Value;
    pub fn gsl_poly_dd_hermite_init(dd: *mut c_double, za: *mut c_double, xa: *const c_double, ya: *const c_double, dya: *const c_double,
        size: size_t) -> enums::Value;
    // Quadratic Equations
    pub fn gsl_poly_solve_quadratic(a: c_double, b: c_double, c: c_double, x0: *mut c_double, x1: *mut c_double) -> c_int;
    pub fn gsl_poly_complex_solve_quadratic(a: c_double, b: c_double, c: c_double, x0: *mut gsl_complex, x1: *mut gsl_complex) -> c_int;
    // Cubic Equations
    pub fn gsl_poly_solve_cubic(a: c_double, b: c_double, c: c_double, x0: *mut c_double, x1: *mut c_double, x2: *mut c_double) -> c_int;
    pub fn gsl_poly_complex_solve_cubic(a: c_double, b: c_double, c: c_double, x0: *mut gsl_complex, x1: *mut gsl_complex,
        x2: *mut gsl_complex) -> c_int;
    // General Polynomial Equations
    pub fn gsl_poly_complex_workspace_alloc(n: size_t) -> *mut gsl_poly_complex_workspace;
    pub fn gsl_poly_complex_workspace_free(w: *mut gsl_poly_complex_workspace);
    pub fn gsl_poly_complex_solve(a: *const c_double, n: size_t, w: *mut gsl_poly_complex_workspace, z: gsl_complex_packed_ptr) -> enums::Value;

    // Discrete Hankel functions
    pub fn gsl_dht_alloc(size: size_t) -> *mut gsl_dht;
    pub fn gsl_dht_init(t: *mut gsl_dht, nu: c_double, xmax: c_double) -> enums::Value;
    pub fn gsl_dht_new(size: size_t, nu: c_double, xmax: c_double) -> *mut gsl_dht;
    pub fn gsl_dht_free(t: *mut gsl_dht);
    pub fn gsl_dht_apply(t: *const gsl_dht, f_in: *const c_double, f_out: *mut c_double) -> enums::Value;
    pub fn gsl_dht_x_sample(t: *const gsl_dht, n: c_int) -> c_double;
    pub fn gsl_dht_k_sample(t: *const gsl_dht, n: c_int) -> c_double;

    // Real Symmetric Matrices
    pub fn gsl_eigen_symm_alloc(n: size_t) -> *mut gsl_eigen_symm_workspace;
    pub fn gsl_eigen_symm_free(w: *mut gsl_eigen_symm_workspace);
    pub fn gsl_eigen_symm(A: *mut gsl_matrix, eval: *mut gsl_vector, w: *mut gsl_eigen_symm_workspace) -> enums::Value;
    pub fn gsl_eigen_symmv_alloc(n: size_t) -> *mut gsl_eigen_symmv_workspace;
    pub fn gsl_eigen_symmv_free(w: *mut gsl_eigen_symmv_workspace);
    pub fn gsl_eigen_symmv(A: *mut gsl_matrix, eval: *mut gsl_vector, evec: *mut gsl_matrix, w: *mut gsl_eigen_symmv_workspace) -> enums::Value;
    // Complex Hermitian Matrices
    pub fn gsl_eigen_herm_alloc(n: size_t) -> *mut gsl_eigen_herm_workspace;
    pub fn gsl_eigen_herm_free(w: *mut gsl_eigen_herm_workspace);
    pub fn gsl_eigen_herm(A: *mut gsl_matrix_complex, eval: *mut gsl_vector, w: *mut gsl_eigen_herm_workspace) -> enums::Value;
    pub fn gsl_eigen_hermv_alloc(n: size_t) -> *mut gsl_eigen_hermv_workspace;
    pub fn gsl_eigen_hermv_free(w: *mut gsl_eigen_hermv_workspace);
    pub fn gsl_eigen_hermv(A: *mut gsl_matrix_complex, eval: *mut gsl_vector, evec: *mut gsl_matrix_complex, w: *mut gsl_eigen_hermv_workspace) -> enums::Value;
    // Real Nonsymmetric Matrices
    pub fn gsl_eigen_nonsymm_alloc(n: size_t) -> *mut gsl_eigen_nonsymm_workspace;
    pub fn gsl_eigen_nonsymm_free(w: *mut gsl_eigen_nonsymm_workspace);
    pub fn gsl_eigen_nonsymm_params(compute_t: c_int, balance: c_int, w: *mut gsl_eigen_nonsymm_workspace);
    pub fn gsl_eigen_nonsymm(A: *mut gsl_matrix, eval: *mut gsl_vector_complex, w: *mut gsl_eigen_nonsymm_workspace) -> enums::Value;
    pub fn gsl_eigen_nonsymm_Z(A: *mut gsl_matrix, eval: *mut gsl_vector_complex, z: *mut gsl_matrix, w: *mut gsl_eigen_nonsymm_workspace) -> enums::Value;
    pub fn gsl_eigen_nonsymmv_alloc(n: size_t) -> *mut gsl_eigen_nonsymmv_workspace;
    pub fn gsl_eigen_nonsymmv_free(w: *mut gsl_eigen_nonsymmv_workspace);
    pub fn gsl_eigen_nonsymmv_params(balance: c_int, w: *mut gsl_eigen_nonsymmv_workspace);
    pub fn gsl_eigen_nonsymmv(A: *mut gsl_matrix, eval: *mut gsl_vector_complex, evec: *mut gsl_matrix_complex, w: *mut gsl_eigen_nonsymmv_workspace) -> enums::Value;
    pub fn gsl_eigen_nonsymmv_Z(A: *mut gsl_matrix, eval: *mut gsl_vector_complex, evec: *mut gsl_matrix_complex, z: *mut gsl_matrix,
        w: *mut gsl_eigen_nonsymmv_workspace) -> enums::Value;
    // Real Generalized Symmetric-Definite Eigensystems
    pub fn gsl_eigen_gensymm_alloc(n: size_t) -> *mut gsl_eigen_gensymm_workspace;
    pub fn gsl_eigen_gensymm_free(w: *mut gsl_eigen_gensymm_workspace);
    pub fn gsl_eigen_gensymm(A: *mut gsl_matrix, B: *mut gsl_matrix, eval: *mut gsl_vector, w: *mut gsl_eigen_gensymm_workspace) -> enums::Value;
    pub fn gsl_eigen_gensymmv_alloc(n: size_t) -> *mut gsl_eigen_gensymmv_workspace;
    pub fn gsl_eigen_gensymmv_free(w: *mut gsl_eigen_gensymmv_workspace);
    pub fn gsl_eigen_gensymmv(A: *mut gsl_matrix, B: *mut gsl_matrix, eval: *mut gsl_vector, evec: *mut gsl_matrix, w: *mut gsl_eigen_gensymmv_workspace) -> enums::Value;
    // Complex Generalized Hermitian-Definite Eigensystems
    pub fn gsl_eigen_genherm_alloc(n: size_t) -> *mut gsl_eigen_genherm_workspace;
    pub fn gsl_eigen_genherm_free(w: *mut gsl_eigen_genherm_workspace);
    pub fn gsl_eigen_genherm(A: *mut gsl_matrix_complex, B: *mut gsl_matrix_complex, eval: *mut gsl_vector,
        w: *mut gsl_eigen_genherm_workspace) -> enums::Value;
    pub fn gsl_eigen_genhermv_alloc(n: size_t) -> *mut gsl_eigen_genhermv_workspace;
    pub fn gsl_eigen_genhermv_free(w: *mut gsl_eigen_genhermv_workspace);
    pub fn gsl_eigen_genhermv(A: *mut gsl_matrix_complex, B: *mut gsl_matrix_complex, eval: *mut gsl_vector, evec: *mut gsl_matrix_complex,
        w: *mut gsl_eigen_genhermv_workspace) -> enums::Value;
    // Real Generalized Nonsymmetric Eigensystems
    pub fn gsl_eigen_gen_alloc(n: size_t) -> *mut gsl_eigen_gen_workspace;
    pub fn gsl_eigen_gen_free(w: *mut gsl_eigen_gen_workspace);
    pub fn gsl_eigen_gen_params(compute_s: c_int, compute_t: c_int, balance: c_int, w: *mut gsl_eigen_gen_workspace);
    pub fn gsl_eigen_gen(A: *mut gsl_matrix, B: *mut gsl_matrix, alpha: *mut gsl_vector_complex, beta: *mut gsl_vector,
        w: *mut gsl_eigen_gen_workspace) -> enums::Value;
    pub fn gsl_eigen_gen_QZ(A: *mut gsl_matrix, B: *mut gsl_matrix, alpha: *mut gsl_vector_complex, beta: *mut gsl_vector,
        Q: *mut gsl_matrix, Z: *mut gsl_matrix, w: *mut gsl_eigen_gen_workspace) -> enums::Value;
    pub fn gsl_eigen_genv_alloc(n: size_t) -> *mut gsl_eigen_genv_workspace;
    pub fn gsl_eigen_genv_free(w: *mut gsl_eigen_genv_workspace);
    pub fn gsl_eigen_genv(A: *mut gsl_matrix, B: *mut gsl_matrix, alpha: *mut gsl_vector_complex, beta: *mut gsl_vector, evec: *mut gsl_matrix_complex,
        w: *mut gsl_eigen_genv_workspace) -> enums::Value;
    pub fn gsl_eigen_genv_QZ(A: *mut gsl_matrix, B: *mut gsl_matrix, alpha: *mut gsl_vector_complex, beta: *mut gsl_vector, evec: *mut gsl_matrix_complex,
        Q: *mut gsl_matrix, Z: *mut gsl_matrix, w: *mut gsl_eigen_genv_workspace) -> enums::Value;
    // Sorting Eigenvalues and Eigenvectors
    pub fn gsl_eigen_symmv_sort(eval: *mut gsl_vector, evec: *mut gsl_matrix, sort_type: ::EigenSort) -> enums::Value;
    pub fn gsl_eigen_hermv_sort(eval: *mut gsl_vector, evec: *mut gsl_matrix_complex, sort_type: ::EigenSort) -> enums::Value;
    pub fn gsl_eigen_nonsymmv_sort(eval: *mut gsl_vector_complex, evec: *mut gsl_matrix_complex, sort_type: ::EigenSort) -> enums::Value;
    pub fn gsl_eigen_gensymmv_sort(eval: *mut gsl_vector, evec: *mut gsl_matrix, sort_type: ::EigenSort) -> enums::Value;
    pub fn gsl_eigen_genhermv_sort(eval: *mut gsl_vector, evec: *mut gsl_matrix_complex, sort_type: ::EigenSort) -> enums::Value;
    pub fn gsl_eigen_genv_sort(alpha: *mut gsl_vector_complex, beta: *mut gsl_vector, evec: *mut gsl_matrix_complex, sort_type: ::EigenSort) -> enums::Value;

    // Fast Fourier Transforms
    // Radix-2 FFT routines for complex data
    pub fn gsl_fft_complex_radix2_forward(data: gsl_complex_packed_array, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_complex_radix2_transform(data: gsl_complex_packed_array, stride: size_t, n: size_t, sign: ::FftDirection) -> enums::Value;
    pub fn gsl_fft_complex_radix2_backward(data: gsl_complex_packed_array, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_complex_radix2_inverse(data: gsl_complex_packed_array, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_complex_radix2_dif_forward(data: gsl_complex_packed_array, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_complex_radix2_dif_transform(data: gsl_complex_packed_array, stride: size_t, n: size_t, sign: ::FftDirection) -> enums::Value;
    pub fn gsl_fft_complex_radix2_dif_backward(data: gsl_complex_packed_array, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_complex_radix2_dif_inverse(data: gsl_complex_packed_array, stride: size_t, n: size_t) -> enums::Value;
    // Mixed-radix FFT routines for complex data
    pub fn gsl_fft_complex_wavetable_alloc(n: size_t) -> *mut gsl_fft_complex_wavetable;
    pub fn gsl_fft_complex_wavetable_free(w: *mut gsl_fft_complex_wavetable);
    pub fn gsl_fft_complex_workspace_alloc(n: size_t) -> *mut gsl_fft_complex_workspace;
    pub fn gsl_fft_complex_workspace_free(w: *mut gsl_fft_complex_workspace);
    pub fn gsl_fft_complex_forward(data: gsl_complex_packed_array, stride: size_t, n: size_t, wavetable: *const gsl_fft_complex_wavetable,
        work: *mut gsl_fft_complex_workspace) -> enums::Value;
    pub fn gsl_fft_complex_transform(data: gsl_complex_packed_array, stride: size_t, n: size_t, wavetable: *const gsl_fft_complex_wavetable,
        work: *mut gsl_fft_complex_workspace, sign: ::FftDirection) -> enums::Value;
    pub fn gsl_fft_complex_backward(data: gsl_complex_packed_array, stride: size_t, n: size_t, wavetable: *const gsl_fft_complex_wavetable,
        work: *mut gsl_fft_complex_workspace) -> enums::Value;
    pub fn gsl_fft_complex_inverse(data: gsl_complex_packed_array, stride: size_t, n: size_t, wavetable: *const gsl_fft_complex_wavetable,
        work: *mut gsl_fft_complex_workspace) -> enums::Value;
    // Radix-2 FFT routines for real data
    pub fn gsl_fft_real_radix2_transform(data: *mut c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_halfcomplex_radix2_inverse(data: *mut c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_halfcomplex_radix2_backward(data: *mut c_double, stride: size_t, n: size_t) -> enums::Value;
    pub fn gsl_fft_halfcomplex_radix2_unpack(halfcomplex_coefficient: *mut c_double, complex_coefficient: gsl_complex_packed_array,
        stride: size_t, n: size_t) -> enums::Value;

    // Histograms
    // Histogram allocation
    pub fn gsl_histogram_alloc(n: size_t) -> *mut gsl_histogram;
    pub fn gsl_histogram_set_ranges(h: *mut gsl_histogram, range: *const c_double, size: size_t) -> enums::Value;
    pub fn gsl_histogram_set_ranges_uniform(h: *mut gsl_histogram, xmin: c_double, xmax: c_double) -> enums::Value;
    pub fn gsl_histogram_free(h: *mut gsl_histogram);
    // Copying Histograms
    pub fn gsl_histogram_memcpy(dest: *mut gsl_histogram, src: *const gsl_histogram) -> enums::Value;
    pub fn gsl_histogram_clone(src: *const gsl_histogram) -> *mut gsl_histogram;
    // Updating and accessing histogram elements
    pub fn gsl_histogram_increment(h: *mut gsl_histogram, x: c_double) -> enums::Value;
    pub fn gsl_histogram_accumulate(h: *mut gsl_histogram, x: c_double, weigth: c_double) -> enums::Value;
    pub fn gsl_histogram_get(h: *const gsl_histogram, i: size_t) -> c_double;
    pub fn gsl_histogram_get_range(h: *const gsl_histogram, i: size_t, lower: *mut c_double, upper: *mut c_double) -> enums::Value;
    pub fn gsl_histogram_max(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_min(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_bins(h: *const gsl_histogram) -> size_t;
    pub fn gsl_histogram_reset(h: *mut gsl_histogram);
    pub fn gsl_histogram_find(h: *const gsl_histogram, x: c_double, i: *mut size_t) -> enums::Value;
    pub fn gsl_histogram_max_val(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_max_bin(h: *const gsl_histogram) -> size_t;
    pub fn gsl_histogram_min_val(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_min_bin(h: *const gsl_histogram) -> size_t;
    pub fn gsl_histogram_mean(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_sigma(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_sum(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_equal_bins_p(h1: *const gsl_histogram, h2: *const gsl_histogram) -> c_int;
    pub fn gsl_histogram_add(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> enums::Value;
    pub fn gsl_histogram_sub(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> enums::Value;
    pub fn gsl_histogram_mul(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> enums::Value;
    pub fn gsl_histogram_div(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> enums::Value;
    pub fn gsl_histogram_scale(h1: *mut gsl_histogram, scale: c_double) -> enums::Value;
    pub fn gsl_histogram_shift(h1: *mut gsl_histogram, offset: c_double) -> enums::Value;
    // The histogram probability distribution struct
    pub fn gsl_histogram_pdf_alloc(n: size_t) -> *mut gsl_histogram_pdf;
    pub fn gsl_histogram_pdf_init(p: *mut gsl_histogram_pdf, h: *const gsl_histogram) -> enums::Value;
    pub fn gsl_histogram_pdf_free(p: *mut gsl_histogram_pdf);
    pub fn gsl_histogram_pdf_sample(p: *const gsl_histogram_pdf, r: c_double) -> c_double;
    // 2D Histogram allocation
    pub fn gsl_histogram2d_alloc(nx: size_t, ny: size_t) -> *mut gsl_histogram2d;
    pub fn gsl_histogram2d_set_ranges(h: *mut gsl_histogram2d, xrange: *const c_double, xsize: size_t, yrange: *const c_double,
        ysize: size_t) -> enums::Value;
    pub fn gsl_histogram2d_set_ranges_uniform(h: *mut gsl_histogram2d, xmin: c_double, xmax: c_double, ymin: c_double,
        ymax: c_double) -> enums::Value;
    pub fn gsl_histogram2d_free(h: *mut gsl_histogram2d);
    pub fn gsl_histogram2d_memcpy(dest: *mut gsl_histogram2d, src: *const gsl_histogram2d) -> enums::Value;
    pub fn gsl_histogram2d_clone(src: *const gsl_histogram2d) -> *mut gsl_histogram2d;
    pub fn gsl_histogram2d_increment(h: *mut gsl_histogram2d, x: c_double, y: c_double) -> enums::Value;
    pub fn gsl_histogram2d_accumulate(h: *mut gsl_histogram2d, x: c_double, y: c_double, weight: c_double) -> enums::Value;
    pub fn gsl_histogram2d_get(h: *const gsl_histogram2d, i: size_t, j: size_t) -> c_double;
    pub fn gsl_histogram2d_get_xrange(h: *const gsl_histogram2d, i: size_t, xlower: *mut c_double, xupper: *mut c_double) -> enums::Value;
    pub fn gsl_histogram2d_get_yrange(h: *const gsl_histogram2d, j: size_t, ylower: *mut c_double, yupper: *mut c_double) -> enums::Value;
    pub fn gsl_histogram2d_xmax(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_xmin(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_nx(h: *const gsl_histogram2d) -> size_t;
    pub fn gsl_histogram2d_ymax(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_ymin(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_ny(h: *const gsl_histogram2d) -> size_t;
    pub fn gsl_histogram2d_reset(h: *mut gsl_histogram2d);
    pub fn gsl_histogram2d_find(h: *const gsl_histogram2d, x: c_double, y: c_double, i: *mut size_t, j: *mut size_t) -> enums::Value;
    // 2D Histogram Statistics
    pub fn gsl_histogram2d_max_val(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_max_bin(h: *const gsl_histogram2d, i: *mut size_t, j: *mut size_t);
    pub fn gsl_histogram2d_min_val(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_min_bin(h: *const gsl_histogram2d, i: *mut size_t, j: *mut size_t);
    pub fn gsl_histogram2d_xmean(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_ymean(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_xsigma(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_ysigma(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_cov(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_sum(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_equal_bins_p(h1: *const gsl_histogram2d, h2: *const gsl_histogram2d) -> c_int;
    pub fn gsl_histogram2d_add(h1: *mut gsl_histogram2d, h2: *const gsl_histogram2d) -> enums::Value;
    pub fn gsl_histogram2d_sub(h1: *mut gsl_histogram2d, h2: *const gsl_histogram2d) -> enums::Value;
    pub fn gsl_histogram2d_mul(h1: *mut gsl_histogram2d, h2: *const gsl_histogram2d) -> enums::Value;
    pub fn gsl_histogram2d_div(h1: *mut gsl_histogram2d, h2: *const gsl_histogram2d) -> enums::Value;
    pub fn gsl_histogram2d_scale(h1: *mut gsl_histogram2d, scale: c_double) -> enums::Value;
    pub fn gsl_histogram2d_shift(h1: *mut gsl_histogram2d, offset: c_double) -> enums::Value;
    // Resampling from 2D histograms
    pub fn gsl_histogram2d_pdf_alloc(nx: size_t, ny: size_t) -> *mut gsl_histogram2d_pdf;
    pub fn gsl_histogram2d_pdf_init(p: *mut gsl_histogram2d_pdf, h: *const gsl_histogram2d) -> enums::Value;
    pub fn gsl_histogram2d_pdf_free(p: *mut gsl_histogram2d_pdf);
    pub fn gsl_histogram2d_pdf_sample(p: *const gsl_histogram2d_pdf, r1: c_double, r2: c_double, x: *mut c_double, y: *mut c_double) -> enums::Value;

    // QAG adaptive integration
    pub fn gsl_integration_workspace_alloc(n: size_t) -> *mut gsl_integration_workspace;
    pub fn gsl_integration_workspace_free(w: *mut gsl_integration_workspace);
    // QAWS adaptive integration for singular functions
    pub fn gsl_integration_qaws_table_alloc(alpha: c_double, beta: c_double, mu: c_int, nu: c_int) -> *mut gsl_integration_qaws_table;
    pub fn gsl_integration_qaws_table_set(t: *mut gsl_integration_qaws_table, alpha: c_double, beta: c_double, mu: c_int, nu: c_int) -> enums::Value;
    pub fn gsl_integration_qaws_table_free(t: *mut gsl_integration_qaws_table);
    // QAWO adaptive integration for oscillatory functions
    pub fn gsl_integration_qawo_table_alloc(omega: c_double, l: c_double, sine: ::IntegrationQawo, n: size_t) -> *mut gsl_integration_qawo_table;
    pub fn gsl_integration_qawo_table_set(t: *mut gsl_integration_qawo_table, omega: c_double, l: c_double, sine: ::IntegrationQawo) -> enums::Value;
    pub fn gsl_integration_qawo_table_set_length(t: *mut gsl_integration_qawo_table, l: c_double) -> enums::Value;
    pub fn gsl_integration_qawo_table_free(t: *mut gsl_integration_qawo_table);
    // CQUAD doubly-adaptive integration
    pub fn gsl_integration_cquad_workspace_alloc(n: size_t) -> *mut gsl_integration_cquad_workspace;
    pub fn gsl_integration_cquad_workspace_free(w: *mut gsl_integration_cquad_workspace);
    // Gauss-Legendre integration
    pub fn gsl_integration_glfixed_table_alloc(n: size_t) -> *mut gsl_integration_glfixed_table;
    pub fn gsl_integration_glfixed_point(a: c_double, b: c_double, i: size_t, xi: *mut c_double, wi: *mut c_double, t: *const gsl_integration_glfixed_table) -> enums::Value;
    pub fn gsl_integration_glfixed_table_free(t: *mut gsl_integration_glfixed_table);

    // Interpolation Functions
    pub fn gsl_interp_alloc(t: *const gsl_interp_type, size: size_t) -> *mut gsl_interp;
    pub fn gsl_interp_init(interp: *mut gsl_interp, xa: *const c_double, ya: *const c_double, size: size_t) -> enums::Value;
    pub fn gsl_interp_free(interp: *mut gsl_interp);
    pub fn gsl_interp_min_size(interp: *const gsl_interp) -> c_uint;
    pub fn gsl_interp_name(interp: *const gsl_interp) -> *const c_char;
    // Interpolation Types
    pub fn gsl_interp_type_min_size(t: *const gsl_interp_type) -> c_uint;
    // Index Look-up and Acceleration
    pub fn gsl_interp_accel_find(a: *mut ::InterpAccel, x_array: *const c_double, size: size_t, x: c_double) -> size_t;
    pub fn gsl_interp_bsearch(x_array: *const c_double, x: c_double, index_lo: size_t, index_hi: size_t) -> size_t;
    // Evaluation of Interpolating Functions
    pub fn gsl_interp_eval(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, x: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_interp_eval_e(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, x: c_double, acc: *mut ::InterpAccel, y: *mut c_double) -> enums::Value;
    pub fn gsl_interp_eval_deriv(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, x: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_interp_eval_deriv_e(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, x: c_double, acc: *mut ::InterpAccel, d: *mut c_double) -> enums::Value;
    pub fn gsl_interp_eval_deriv2(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, x: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_interp_eval_deriv2_e(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, x: c_double, acc: *mut ::InterpAccel, d2: *mut c_double) -> enums::Value;
    pub fn gsl_interp_eval_integ(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, a: c_double, b: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_interp_eval_integ_e(interp: *const gsl_interp, xa: *const c_double, ya: *const c_double, a: c_double, b: c_double, acc: *mut ::InterpAccel, result: *mut c_double) -> enums::Value;
    // Higher-level Interface
    pub fn gsl_spline_alloc(t: *const gsl_interp_type, size: size_t) -> *mut gsl_spline;
    pub fn gsl_spline_init(spline: *mut gsl_spline, xa: *const c_double, ya: *const c_double, size: size_t) -> enums::Value;
    pub fn gsl_spline_free(spline: *mut gsl_spline);
    pub fn gsl_spline_min_size(spline: *const gsl_spline) -> c_uint;
    pub fn gsl_spline_name(spline: *const gsl_spline) -> *const c_char;
    pub fn gsl_spline_eval(spline: *const gsl_spline, x: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_spline_eval_e(spline: *const gsl_spline, x: c_double, acc: *mut ::InterpAccel, y: *mut c_double) -> enums::Value;
    pub fn gsl_spline_eval_deriv(spline: *const gsl_spline, x: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_spline_eval_deriv_e(spline: *const gsl_spline, x: c_double, acc: *mut ::InterpAccel, d: *mut c_double) -> enums::Value;
    pub fn gsl_spline_eval_deriv2(spline: *const gsl_spline, x: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_spline_eval_deriv2_e(spline: *const gsl_spline, x: c_double, acc: *mut ::InterpAccel, d2: *mut c_double) -> enums::Value;
    pub fn gsl_spline_eval_integ(spline: *const gsl_spline, a: c_double, b: c_double, acc: *mut ::InterpAccel) -> c_double;
    pub fn gsl_spline_eval_integ_e(spline: *const gsl_spline, a: c_double, b: c_double, acc: *mut ::InterpAccel, result: *mut c_double) -> enums::Value;

    // linear algebra
    // LU Decomposition
    pub fn gsl_linalg_LU_decomp(a: *mut gsl_matrix, p: *mut gsl_permutation, signum: *mut c_int) -> enums::Value;
    pub fn gsl_linalg_complex_LU_decomp(a: *mut gsl_matrix_complex, p: *mut gsl_permutation, signum: *mut c_int) -> enums::Value;
    pub fn gsl_linalg_LU_solve(lu: *const gsl_matrix, p: *const gsl_permutation, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_complex_LU_solve(lu: *const gsl_matrix_complex, p: *const gsl_permutation, b: *const gsl_vector_complex, x: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_linalg_LU_svx(lu: *const gsl_matrix, p: *const gsl_permutation, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_complex_LU_svx(lu: *const gsl_matrix_complex, p: *const gsl_permutation, x: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_linalg_LU_refine(a: *const gsl_matrix, lu: *const gsl_matrix, p: *const gsl_permutation, b: *const gsl_vector, x: *mut gsl_vector,
        residual: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_complex_LU_refine(a: *const gsl_matrix_complex, lu: *const gsl_matrix_complex, p: *const gsl_permutation, b: *const gsl_vector_complex,
        x: *mut gsl_vector_complex, residual: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_linalg_LU_invert(lu: *const gsl_matrix, p: *const gsl_permutation, inverse: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_complex_LU_invert(lu: *const gsl_matrix_complex, p: *const gsl_permutation, inverse: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_linalg_LU_det(lu: *mut gsl_matrix, signum: c_int) -> c_double;
    pub fn gsl_linalg_complex_LU_det(lu: *mut gsl_matrix_complex, signum: c_int) -> gsl_complex;
    pub fn gsl_linalg_LU_lndet(lu: *mut gsl_matrix) -> c_double;
    pub fn gsl_linalg_complex_LU_lndet(lu: *mut gsl_matrix_complex) -> c_double;
    pub fn gsl_linalg_LU_sgndet(lu: *mut gsl_matrix, signum: c_int) -> c_double;
    pub fn gsl_linalg_complex_LU_sgndet(lu: *mut gsl_matrix_complex, signum: c_int) -> gsl_complex;
    // QR Decomposition
    pub fn gsl_linalg_QR_decomp(a: *mut gsl_matrix, tau: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_solve(qr: *const gsl_matrix, tau: *const gsl_vector, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_svx(qr: *const gsl_matrix, tau: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_lssolve(qr: *const gsl_matrix, tau: *const gsl_vector, b: *const gsl_vector, x: *mut gsl_vector,
        residual: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_QTvec(qr: *const gsl_matrix, tau: *const gsl_vector, v: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_Qvec(qr: *const gsl_matrix, tau: *const gsl_vector, v: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_QTmat(qr: *const gsl_matrix, tau: *const gsl_vector, v: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_QR_Rsolve(qr: *const gsl_matrix, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_Rsvx(qr: *const gsl_matrix, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_unpack(qr: *const gsl_matrix, tau: *const gsl_vector, q: *mut gsl_matrix, r: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_QR_QRsolve(q: *mut gsl_matrix, r: *mut gsl_matrix, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QR_update(q: *mut gsl_matrix, r: *mut gsl_matrix, w: *mut gsl_vector, v: *const gsl_vector) -> enums::Value;
    pub fn gsl_linalg_R_solve(r: *const gsl_matrix, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_R_svx(r: *const gsl_matrix, x: *mut gsl_vector) -> enums::Value;
    // QR Decomposition with Column Pivoting
    pub fn gsl_linalg_QRPT_decomp(a: *mut gsl_matrix, tau: *mut gsl_vector, p: *mut gsl_permutation, signum: *mut c_int,
        norm: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QRPT_decomp2(a: *const gsl_matrix, q: *mut gsl_matrix, r: *mut gsl_matrix, tau: *mut gsl_vector, p: *mut gsl_permutation,
        signum: *mut c_int, norm: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QRPT_solve(qr: *const gsl_matrix, tau: *const gsl_vector, p: *const gsl_permutation, b: *const gsl_vector,
        x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QRPT_svx(qr: *const gsl_matrix, tau: *const gsl_vector, p: *const gsl_permutation, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QRPT_QRsolve(q: *const gsl_matrix, r: *const gsl_matrix, p: *const gsl_permutation, b: *const gsl_vector,
        x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QRPT_update(q: *const gsl_matrix, r: *const gsl_matrix, p: *const gsl_permutation, w: *mut gsl_vector,
        v: *const gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QRPT_Rsolve(qr: *const gsl_matrix, p: *const gsl_permutation, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_QRPT_Rsvx(qr: *const gsl_matrix, p: *const gsl_permutation, x: *mut gsl_vector) -> enums::Value;
    // Singular Value Decomposition
    pub fn gsl_linalg_SV_decomp(a: *mut gsl_matrix, v: *mut gsl_matrix, s: *mut gsl_vector, work: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_SV_decomp_mod(a: *mut gsl_matrix, x: *mut gsl_matrix, v: *mut gsl_matrix, s: *mut gsl_vector, work: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_SV_decomp_jacobi(a: *mut gsl_matrix, v: *mut gsl_matrix, s: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_SV_solve(u: *const gsl_matrix, v: *const gsl_matrix, s: *const gsl_vector, b: *const gsl_vector,
        x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_SV_leverage(u: *const gsl_matrix, h: *mut gsl_vector) -> enums::Value;
    // Cholesky Decomposition
    pub fn gsl_linalg_cholesky_decomp(a: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_complex_cholesky_decomp(a: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_linalg_cholesky_solve(cholesky: *const gsl_matrix, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_complex_cholesky_solve(cholesky: *const gsl_matrix_complex, b: *const gsl_vector_complex,
        x: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_linalg_cholesky_svx(cholesky: *const gsl_matrix, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_complex_cholesky_svx(cholesky: *const gsl_matrix_complex, x: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_linalg_cholesky_invert(cholesky: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_complex_cholesky_invert(cholesky: *mut gsl_matrix_complex) -> enums::Value;
    // Tridiagonal Decomposition of Real Symmetric Matrices
    pub fn gsl_linalg_symmtd_decomp(a: *mut gsl_matrix, tau: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_symmtd_unpack(a: *const gsl_matrix, tau: *const gsl_vector, q: *mut gsl_matrix, diag: *mut gsl_vector,
        subdiag: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_symmtd_unpack_T(a: *const gsl_matrix, diag: *mut gsl_vector, subdiag: *mut gsl_vector) -> enums::Value;
    // Tridiagonal Decomposition of Hermitian Matrices
    pub fn gsl_linalg_hermtd_decomp(a: *mut gsl_matrix_complex, tau: *mut gsl_vector_complex) -> enums::Value;
    pub fn gsl_linalg_hermtd_unpack(a: *const gsl_matrix_complex, tau: *const gsl_vector_complex, u: *mut gsl_matrix_complex, diag: *mut gsl_vector,
        subdiag: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_hermtd_unpack_T(a: *const gsl_matrix_complex, diag: *mut gsl_vector, subdiag: *mut gsl_vector) -> enums::Value;
    // Hessenberg Decomposition of Real Matrices
    pub fn gsl_linalg_hessenberg_decomp(a: *mut gsl_matrix, tau: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_hessenberg_unpack(h: *mut gsl_matrix, tau: *mut gsl_vector, u: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_hessenberg_unpack_accum(h: *mut gsl_matrix, tau: *mut gsl_vector, v: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_hessenberg_set_zero(a: *mut gsl_matrix) -> enums::Value;
    // Hessenberg-Triangular Decomposition of Real Matrices
    pub fn gsl_linalg_hesstri_decomp(a: *mut gsl_matrix, b: *mut gsl_matrix, u: *mut gsl_matrix, v: *mut gsl_matrix,
        work: *mut gsl_vector) -> enums::Value;
    // Bidiagonalization
    pub fn gsl_linalg_bidiag_decomp(a: *mut gsl_matrix, tau_u: *mut gsl_vector, tau_v: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_bidiag_unpack(a: *const gsl_matrix, tau_u: *const gsl_vector, u: *mut gsl_matrix, tau_v: *const gsl_vector,
        v: *mut gsl_matrix, diag: *mut gsl_vector, superdiag: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_bidiag_unpack2(a: *mut gsl_matrix, tau_u: *mut gsl_vector, tau_v: *mut gsl_vector, v: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_bidiag_unpack_B(a: *const gsl_matrix, diag: *mut gsl_vector, superdiag: *mut gsl_vector) -> enums::Value;
    // Householder Transformations
    pub fn gsl_linalg_householder_transform(v: *mut gsl_vector) -> c_double;
    pub fn gsl_linalg_complex_householder_transform(v: *mut gsl_vector_complex) -> gsl_complex;
    pub fn gsl_linalg_householder_hm(tau: c_double, v: *const gsl_vector, a: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_complex_householder_hm(tau: gsl_complex, v: *const gsl_vector_complex, a: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_linalg_householder_mh(tau: c_double, v: *const gsl_vector, a: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_complex_householder_mh(tau: gsl_complex, v: *const gsl_vector_complex, a: *mut gsl_matrix_complex) -> enums::Value;
    pub fn gsl_linalg_householder_hv(tau: c_double, v: *const gsl_vector, w: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_linalg_complex_householder_hv(tau: gsl_complex, v: *const gsl_vector_complex, w: *mut gsl_matrix_complex) -> enums::Value;
    // Householder solver for linear systems
    pub fn gsl_linalg_HH_solve(a: *mut gsl_matrix, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_HH_svx(a: *mut gsl_matrix, x: *mut gsl_vector) -> enums::Value;
    // Tridiagonal Systems
    pub fn gsl_linalg_solve_tridiag(diag: *const gsl_vector, e: *const gsl_vector, f: *const gsl_vector, b: *const gsl_vector,
        x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_solve_symm_tridiag(diag: *const gsl_vector, e: *const gsl_vector, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_solve_cyc_tridiag(diag: *const gsl_vector, e: *const gsl_vector, f: *const gsl_vector, b: *const gsl_vector,
        x: *mut gsl_vector) -> enums::Value;
    pub fn gsl_linalg_solve_symm_cyc_tridiag(diag: *const gsl_vector, e: *const gsl_vector, b: *const gsl_vector, x: *mut gsl_vector) -> enums::Value;
    // Balancing
    pub fn gsl_linalg_balance_matrix(a: *mut gsl_matrix, d: *mut gsl_vector) -> enums::Value;

    // Initializing the Minimizer
    /*pub fn gsl_min_fminimizer_alloc(t: *const gsl_min_fminimizer_type) -> *mut gsl_min_fminimizer;
    pub fn gsl_min_fminimizer_set(s: *mut gsl_min_fminimizer, f: *mut gsl_function, x_minimum: c_double, x_lower: c_double,
        x_upper: c_double) -> enums::Value;
    pub fn gsl_min_fminimizer_set_with_values(s: *mut gsl_min_fminimizer, f: *mut gsl_function, x_minimum: c_double, f_minimum: c_double,
        x_lower: c_double, f_lower: c_double, x_upper: c_double, f_upper: c_double) -> enums::Value;
    pub fn gsl_min_fminimizer_free(s: *mut gsl_min_fminimizer);
    pub fn gsl_min_fminimizer_name(s: *const gsl_min_fminimizer) -> *const c_char;*/
    pub fn gsl_min_test_interval(x_lower: c_double, x_upper: c_double, epsabs: c_double, epsrel: c_double) -> enums::Value;

    // Monte Carlo
    // PLAIN Monte Carlo
    pub fn gsl_monte_plain_alloc(dim: size_t) -> *mut gsl_monte_plain_state;
    pub fn gsl_monte_plain_init(s: *mut gsl_monte_plain_state) -> enums::Value;
    pub fn gsl_monte_plain_free(s: *mut gsl_monte_plain_state);
    pub fn gsl_monte_plain_integrate(f: *mut c_void, xl: *const c_double, xu: *const c_double, dim: size_t,
                                     calls: size_t, r: *mut gsl_rng, s: *mut gsl_monte_plain_state,
                                     result: *mut c_double, abserr: *mut c_double) -> enums::Value;
    // MISER
    pub fn gsl_monte_miser_alloc(dim: size_t) -> *mut gsl_monte_miser_state;
    pub fn gsl_monte_miser_init(s: *mut gsl_monte_miser_state) -> enums::Value;
    pub fn gsl_monte_miser_free(s: *mut gsl_monte_miser_state);
    pub fn gsl_monte_miser_integrate(f: *mut c_void, xl: *const c_double, xu: *const c_double, dim: size_t,
                                     calls: size_t, r: *mut gsl_rng, s: *mut gsl_monte_miser_state,
                                     result: *mut c_double, abserr: *mut c_double) -> enums::Value;
    pub fn gsl_monte_miser_params_get(s: *mut gsl_monte_miser_state, m: *mut ::MiserParams);
    pub fn gsl_monte_miser_params_set(s: *mut gsl_monte_miser_state, m: *const ::MiserParams);
    // VEGAS
    pub fn gsl_monte_vegas_alloc(dim: size_t) -> *mut gsl_monte_vegas_state;
    pub fn gsl_monte_vegas_init(s: *mut gsl_monte_vegas_state) -> enums::Value;
    pub fn gsl_monte_vegas_free(s: *mut gsl_monte_vegas_state);
    pub fn gsl_monte_vegas_integrate(f: *mut c_void, xl: *const c_double, xu: *const c_double, dim: size_t,
                                     calls: size_t, r: *mut gsl_rng, s: *mut gsl_monte_vegas_state,
                                     result: *mut c_double, abserr: *mut c_double) -> enums::Value;
    pub fn gsl_monte_vegas_chisq(s: *const gsl_monte_vegas_state) -> c_double;
    pub fn gsl_monte_vegas_runval(s: *const gsl_monte_vegas_state, result: *mut c_double, sigma: *mut c_double);
    //pub fn gsl_monte_vegas_params_get(s: *const gsl_monte_vegas_state, params: *mut ::VegasParams);
    //pub fn gsl_monte_vegas_params_set(s: *mut gsl_monte_vegas_state, params: *const ::VegasParams);

    // N-tuples
    // Creating ntuples
    pub fn gsl_ntuple_create(filename: *mut c_char, ntuple_data: *mut c_void, size: size_t) -> *mut gsl_ntuple;
    // Opening an existing ntuple file
    pub fn gsl_ntuple_open(filename: *mut c_char, ntuple_data: *mut c_void, size: size_t) -> *mut gsl_ntuple;
    // Writing ntuples
    pub fn gsl_ntuple_write(ntuple: *mut gsl_ntuple) -> enums::Value;
    pub fn gsl_ntuple_bookdata(ntuple: *mut gsl_ntuple) -> enums::Value;
    // Reading ntuples
    pub fn gsl_ntuple_read(ntuple: *mut gsl_ntuple) -> enums::Value;
    // Closing an ntuple file
    pub fn gsl_ntuple_close(ntuple: *mut gsl_ntuple) -> enums::Value;

    // Multisets
    // Multiset allocation
    pub fn gsl_multiset_alloc(n: size_t, k: size_t) -> *mut gsl_multiset;
    pub fn gsl_multiset_calloc(n: size_t, k: size_t) -> *mut gsl_multiset;
    pub fn gsl_multiset_init_first(c: *mut gsl_multiset);
    pub fn gsl_multiset_init_last(c: *mut gsl_multiset);
    pub fn gsl_multiset_free(c: *mut gsl_multiset);
    pub fn gsl_multiset_memcpy(dest: *mut gsl_multiset, src: *const gsl_multiset) -> enums::Value;
    // Accessing multiset elements
    pub fn gsl_multiset_get(c: *const gsl_multiset, i: size_t) -> size_t;
    // Multiset properties
    pub fn gsl_multiset_n(c: *const gsl_multiset) -> size_t;
    pub fn gsl_multiset_k(c: *const gsl_multiset) -> size_t;
    //pub fn gsl_multiset_data(c: *const gsl_multiset) -> *mut size_t;
    pub fn gsl_multiset_valid(c: *mut gsl_multiset) -> enums::Value;
    // Multiset functions
    pub fn gsl_multiset_next(c: *mut gsl_multiset) -> enums::Value;
    pub fn gsl_multiset_prev(c: *mut gsl_multiset) -> enums::Value;

    // Ordinary Differential Equations
    // Stepping Functions
    pub fn gsl_odeiv2_step_alloc(t: *const gsl_odeiv2_step_type, dim: size_t) -> *mut gsl_odeiv2_step;
    pub fn gsl_odeiv2_step_reset(s: *mut gsl_odeiv2_step) -> enums::Value;
    pub fn gsl_odeiv2_step_free(s: *mut gsl_odeiv2_step);
    pub fn gsl_odeiv2_step_name(s: *mut gsl_odeiv2_step) -> *const c_char;
    pub fn gsl_odeiv2_step_order(s: *const gsl_odeiv2_step) -> c_uint;
    pub fn gsl_odeiv2_step_set_driver(s: *mut gsl_odeiv2_step, d: *const gsl_odeiv2_driver) -> enums::Value;
    pub fn gsl_odeiv2_step_apply(s: *mut gsl_odeiv2_step, t: c_double, h: c_double, y: *mut c_double, yerr: *mut c_double,
        dydt_in: *const c_double, dydt_out: *mut c_double, sys: *const gsl_odeiv2_system) -> enums::Value;
    // Adaptive Step-size Control
    pub fn gsl_odeiv2_control_standard_new(eps_abs: c_double, eps_rel: c_double, a_y: c_double, a_dydt: c_double) -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_y_new(eps_abs: c_double, eps_rel: c_double) -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_yp_new(eps_abs: c_double, eps_rel: c_double) -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_scaled_new(eps_abs: c_double, eps_rel: c_double, a_y: c_double, a_dydt: c_double,
        scale_abs: *const c_double, dim: size_t) -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_alloc(t: *const gsl_odeiv2_control_type) -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_init(c: *mut gsl_odeiv2_control, eps_abs: c_double, eps_rel: c_double, a_y: c_double,
        a_dydt: c_double) -> enums::Value;
    pub fn gsl_odeiv2_control_free(c: *mut gsl_odeiv2_control);
    pub fn gsl_odeiv2_control_hadjust(c: *mut gsl_odeiv2_control, s: *mut gsl_odeiv2_step, y: *const c_double, yerr: *const c_double,
        dydt: *const c_double, h: *mut c_double) -> ::ODEiv;
    pub fn gsl_odeiv2_control_name(c: *const gsl_odeiv2_control) -> *const c_char;
    pub fn gsl_odeiv2_control_errlevel(c: *mut gsl_odeiv2_control, y: c_double, dydt: c_double, h: c_double, ind: size_t,
        errlev: *mut c_double) -> enums::Value;
    pub fn gsl_odeiv2_control_set_driver(c: *mut gsl_odeiv2_control, d: *const gsl_odeiv2_driver) -> enums::Value;
    // Evolution
    pub fn gsl_odeiv2_evolve_alloc(dim: size_t) -> *mut gsl_odeiv2_evolve;
    pub fn gsl_odeiv2_evolve_apply(e: *mut gsl_odeiv2_evolve, con: *mut gsl_odeiv2_control, step: *mut gsl_odeiv2_step, sys: *const gsl_odeiv2_system,
        t: *mut c_double, t1: c_double, h: *mut c_double, y: *mut c_double) -> enums::Value;
    pub fn gsl_odeiv2_evolve_apply_fixed_step(e: *mut gsl_odeiv2_evolve, con: *mut gsl_odeiv2_control, step: *mut gsl_odeiv2_step,
        sys: *const gsl_odeiv2_system, t: *mut c_double, h: c_double, y: *mut c_double) -> enums::Value;
    pub fn gsl_odeiv2_evolve_reset(e: *mut gsl_odeiv2_evolve) -> enums::Value;
    pub fn gsl_odeiv2_evolve_free(e: *mut gsl_odeiv2_evolve);
    pub fn gsl_odeiv2_evolve_set_driver(e: *mut gsl_odeiv2_evolve, d: *const gsl_odeiv2_driver) -> enums::Value;
    // Driver
    pub fn gsl_odeiv2_driver_alloc_y_new(sys: *const gsl_odeiv2_system, t: *const gsl_odeiv2_step_type, hstart: c_double, epsabs: c_double,
        epsrel: c_double) -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_alloc_yp_new(sys: *const gsl_odeiv2_system, t: *const gsl_odeiv2_step_type, hstart: c_double, epsabs: c_double,
        epsrel: c_double) -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_alloc_standard_new(sys: *const gsl_odeiv2_system, t: *const gsl_odeiv2_step_type, hstart: c_double, epsabs: c_double,
        epsrel: c_double, a_y: c_double, a_dydt: c_double) -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_alloc_scaled_new(sys: *const gsl_odeiv2_system, t: *const gsl_odeiv2_step_type, hstart: c_double, epsabs: c_double,
        epsrel: c_double, a_y: c_double, a_dydt: c_double, scale_abs: *const c_double) -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_set_hmin(d: *mut gsl_odeiv2_driver, hmin: c_double) -> enums::Value;
    pub fn gsl_odeiv2_driver_set_hmax(d: *mut gsl_odeiv2_driver, hmax: c_double) -> enums::Value;
    pub fn gsl_odeiv2_driver_set_nmax(d: *mut gsl_odeiv2_driver, nmax: usize) -> enums::Value;
    pub fn gsl_odeiv2_driver_apply(d: *mut gsl_odeiv2_driver, t: *mut c_double, t1: c_double, y: *mut c_double) -> enums::Value;
    pub fn gsl_odeiv2_driver_apply_fixed_step(d: *mut gsl_odeiv2_driver, t: *mut c_double, h: c_double, n: usize,
        y: *mut c_double) -> enums::Value;
    pub fn gsl_odeiv2_driver_reset(d: *mut gsl_odeiv2_driver) -> enums::Value;
    pub fn gsl_odeiv2_driver_reset_hstart(d: *mut gsl_odeiv2_driver, hstart: c_double) -> enums::Value;
    pub fn gsl_odeiv2_driver_free(d: *mut gsl_odeiv2_driver);

    // Quasi-Random Sequences
    // Quasi-random number generator initialization
    pub fn gsl_qrng_alloc(t: *const gsl_qrng_type, d: c_uint) -> *mut gsl_qrng;
    pub fn gsl_qrng_free(q: *mut gsl_qrng);
    pub fn gsl_qrng_init(q: *mut gsl_qrng);
    // Sampling from a quasi-random number generator
    pub fn gsl_qrng_get(q: *const gsl_qrng, x: *mut c_double) -> enums::Value;
    // Auxiliary quasi-random number generator functions
    pub fn gsl_qrng_name(q: *const gsl_qrng) -> *const c_char;
    pub fn gsl_qrng_size(q: *const gsl_qrng) -> size_t;
    pub fn gsl_qrng_state(q: *const gsl_qrng) -> *mut c_void;
    // Saving and resorting quasi-random number generator state
    pub fn gsl_qrng_memcpy(dest: *mut gsl_qrng, src: *const gsl_qrng) -> enums::Value;
    pub fn gsl_qrng_clone(q: *const gsl_qrng) -> *mut gsl_qrng;

    // Statistics
    // Mean, Standard Deviation and Variance
    pub fn gsl_stats_mean(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_variance(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_variance_m(data: *const c_double, stride: size_t, n: size_t, mean: c_double) -> c_double;
    pub fn gsl_stats_sd(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_sd_m(data: *const c_double, stride: size_t, n: size_t, mean: c_double) -> c_double;
    pub fn gsl_stats_tss(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_tss_m(data: *const c_double, stride: size_t, n: size_t, mean: c_double) -> c_double;
    pub fn gsl_stats_variance_with_fixed_mean(data: *const c_double, stride: size_t, n: size_t, mean: c_double) -> c_double;
    pub fn gsl_stats_sd_with_fixed_mean(data: *const c_double, stride: size_t, n: size_t, mean: c_double) -> c_double;
    // Absolute deviation
    pub fn gsl_stats_absdev(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_absdev_m(data: *const c_double, stride: size_t, n: size_t, mean: c_double) -> c_double;
    // Higher moments (skewness and kurtosis)
    pub fn gsl_stats_skew(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_skew_m_sd(data: *const c_double, stride: size_t, n: size_t, mean: c_double, sd: c_double) -> c_double;
    pub fn gsl_stats_kurtosis(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_kurtosis_m_sd(data: *const c_double, stride: size_t, n: size_t, mean: c_double, sd: c_double) -> c_double;
    // Autocorrelation
    pub fn gsl_stats_lag1_autocorrelation(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_lag1_autocorrelation_m(data: *const c_double, stride: size_t, n: size_t, mean: c_double) -> c_double;
    // Covariance
    pub fn gsl_stats_covariance(data1: *const c_double, stride1: size_t, data2: *const c_double, stride2: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_covariance_m(data1: *const c_double, stride1: size_t, data2: *const c_double, stride2: size_t, n: size_t,
        mean1: c_double, mean2: c_double) -> c_double;
    // Correlation
    pub fn gsl_stats_correlation(data1: *const c_double, stride1: size_t, data2: *const c_double, stride2: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_spearman(data1: *const c_double, stride1: size_t, data2: *const c_double, stride2: size_t, n: size_t,
        work: *mut c_double) -> c_double;
    // Weighted Samples
    pub fn gsl_stats_wmean(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_wvariance(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_wvariance_m(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double) -> c_double;
    pub fn gsl_stats_wsd(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_wsd_m(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double) -> c_double;
    pub fn gsl_stats_wvariance_with_fixed_mean(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double) -> c_double;
    pub fn gsl_stats_wsd_with_fixed_mean(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double) -> c_double;
    pub fn gsl_stats_wtss(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_wtss_m(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double) -> c_double;
    pub fn gsl_stats_wabsdev(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_wabsdev_m(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double) -> c_double;
    pub fn gsl_stats_wskew(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_wskew_m_sd(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double, wsd: c_double) -> c_double;
    pub fn gsl_stats_wkurtosis(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_wkurtosis_m_sd(w: *const c_double, wstride: size_t, data: *const c_double, stride: size_t, n: size_t,
        wmean: c_double, wsd: c_double) -> c_double;
    // Maximum and Minimum values
    pub fn gsl_stats_max(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_min(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_minmax(min: *mut c_double, max: *mut c_double, data: *const c_double, stride: size_t, n: size_t);
    pub fn gsl_stats_max_index(data: *const c_double, stride: size_t, n: size_t) -> size_t;
    pub fn gsl_stats_min_index(data: *const c_double, stride: size_t, n: size_t) -> size_t;
    pub fn gsl_stats_minmax_index(min: *mut size_t, max: *mut size_t, data: *const c_double, stride: size_t, n: size_t);
    // Median and Percentiles
    pub fn gsl_stats_median_from_sorted_data(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_quantile_from_sorted_data(data: *const c_double, stride: size_t, n: size_t, f: c_double) -> c_double;

    // Series Acceleration
    // Acceleration functions
    pub fn gsl_sum_levin_u_alloc(n: size_t) -> *mut gsl_sum_levin_u_workspace;
    pub fn gsl_sum_levin_u_free(w: *mut gsl_sum_levin_u_workspace);
    pub fn gsl_sum_levin_u_accel(array: *const c_double, array_size: size_t, w: *mut gsl_sum_levin_u_workspace, sum_accel: *mut c_double,
        abserr: *mut c_double) -> enums::Value;
    // Acceleration functions without error estimation
    pub fn gsl_sum_levin_utrunc_alloc(n: size_t) -> *mut gsl_sum_levin_utrunc_workspace;
    pub fn gsl_sum_levin_utrunc_free(w: *mut gsl_sum_levin_utrunc_workspace);
    pub fn gsl_sum_levin_utrunc_accel(array: *const c_double, array_size: size_t, w: *mut gsl_sum_levin_utrunc_workspace,
        sum_accel: *mut c_double, abserr_trunc: *mut c_double) -> enums::Value;

    // Wavelet Transforms
    // Initialization
    pub fn gsl_wavelet_alloc(t: *const gsl_wavelet_type, k: size_t) -> *mut gsl_wavelet;
    pub fn gsl_wavelet_name(w: *const gsl_wavelet) -> *const c_char;
    pub fn gsl_wavelet_free(w: *mut gsl_wavelet);
    pub fn gsl_wavelet_workspace_alloc(n: size_t) -> *mut gsl_wavelet_workspace;
    pub fn gsl_wavelet_workspace_free(w: *mut gsl_wavelet_workspace);
    // Wavelet transforms in one dimension
    pub fn gsl_wavelet_transform(w: *const gsl_wavelet, data: *mut c_double, stride: size_t, n: size_t, dir: ::WaveletDirection,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet_transform_forward(w: *const gsl_wavelet, data: *mut c_double, stride: size_t, n: size_t,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet_transform_inverse(w: *const gsl_wavelet, data: *mut c_double, stride: size_t, n: size_t,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    // Wavelet transforms in two dimension
    pub fn gsl_wavelet2d_transform(w: *const gsl_wavelet, data: *mut c_double, tda: size_t, size1: size_t, size2: size_t,
        dir: ::WaveletDirection, work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_transform_forward(w: *const gsl_wavelet, data: *mut c_double, tda: size_t, size1: size_t, size2: size_t,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_transform_inverse(w: *const gsl_wavelet, data: *mut c_double, tda: size_t, size1: size_t, size2: size_t,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_transform_matrix(w: *const gsl_wavelet, m: *mut gsl_matrix, dir: ::WaveletDirection,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_transform_matrix_forward(w: *const gsl_wavelet, m: *mut gsl_matrix, work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_transform_matrix_inverse(w: *const gsl_wavelet, m: *mut gsl_matrix, work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_nstransform(w: *const gsl_wavelet, data: *mut c_double, tda: size_t, size1: size_t, size2: size_t,
        dir: ::WaveletDirection, work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_nstransform_forward(w: *const gsl_wavelet, data: *mut c_double, tda: size_t, size1: size_t, size2: size_t,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_nstransform_inverse(w: *const gsl_wavelet, data: *mut c_double, tda: size_t, size1: size_t, size2: size_t,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_nstransform_matrix(w: *const gsl_wavelet, m: *mut gsl_matrix, dir: ::WaveletDirection,
        work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_nstransform_matrix_forward(w: *const gsl_wavelet, m: *mut gsl_matrix, work: *mut gsl_wavelet_workspace) -> enums::Value;
    pub fn gsl_wavelet2d_nstransform_matrix_inverse(w: *const gsl_wavelet, m: *mut gsl_matrix, work: *mut gsl_wavelet_workspace) -> enums::Value;

    // multifit
    pub fn gsl_multifit_covar(j: *const gsl_matrix, epsrel: c_double, covar: *mut gsl_matrix) -> enums::Value;
    pub fn gsl_multifit_test_delta(dx: *const gsl_vector, x: *const gsl_vector, epsabs: c_double, epsrel: c_double) -> enums::Value;
    pub fn gsl_multifit_gradient(j: *const gsl_matrix, f: *const gsl_vector, g: *mut gsl_vector) -> enums::Value;
}

#[repr(C)]
pub struct gsl_sf_result {
    pub val: c_double,
    pub err: c_double
}

#[repr(C)]
pub struct gsl_sf_result_e10 {
    pub val: c_double,
    pub err: c_double,
    pub e10: c_int
}

#[repr(C)]
pub struct gsl_vector_float {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_float,
    pub owner: c_int
}

#[repr(C)]
pub struct gsl_block_float {
    pub size: size_t,
    pub data: *mut c_float
}

#[repr(C)]
pub struct gsl_vector {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block,
    pub owner: c_int
}


#[repr(C)]
pub struct gsl_vector_view {
    pub vector: gsl_vector
}

#[repr(C)]
pub struct gsl_block {
    pub size: size_t,
    pub data: *mut c_double
}

#[repr(C)]
pub struct gsl_vector_complex_float {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_complex_float,
    pub owner: c_int
}

#[repr(C)]
pub struct gsl_block_complex_float {
    pub size: size_t,
    pub data: *mut c_float,
}

#[repr(C)]
pub struct gsl_vector_complex {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block_complex,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_block_complex {
    pub size: size_t,
    pub data: *mut c_double,
}

#[repr(C)]
pub struct gsl_complex {
    pub dat: [c_double; 2],
}

#[repr(C)]
pub struct gsl_complex_float {
    pub dat: [c_float; 2],
}

#[repr(C)]
pub struct gsl_matrix {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_matrix_view {
    pub mat: gsl_matrix,
}

#[repr(C)]
pub struct gsl_matrix_float {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_float,
    pub owner: c_int,
}

#[repr(C)]
pub struct gsl_matrix_complex {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block_complex,
    pub owner: c_int
}

#[repr(C)]
pub struct gsl_matrix_complex_float {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_complex_float,
    pub owner: c_int
}

#[repr(C)]
pub struct gsl_sf_mathieu_workspace {
    pub size: size_t,
    pub even_order: size_t,
    pub odd_order: size_t,
    pub extra_values: c_int,
    pub qa: c_double, // allow for caching of results: not implemented yet
    pub qb: c_double, // allow for caching of results: not implemented yet
    pub aa: *mut c_double,
    pub bb: *mut c_double,
    pub dd: *mut c_double,
    pub ee: *mut c_double,
    pub tt: *mut c_double,
    pub e2: *mut c_double,
    pub zz: *mut c_double,
    pub eval: *mut gsl_vector,
    pub evec: *mut gsl_matrix,
    pub wmat: *mut gsl_eigen_symmv_workspace
}

#[repr(C)]
pub struct gsl_eigen_symmv_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
    pub gc: *mut c_double,
    pub gs: *mut c_double
}

#[repr(C)]
pub struct gsl_bspline_workspace {
    pub k: size_t,      // spline order
    pub km1: size_t,    // k - 1 (polynomial order)
    pub l: size_t,      // number of polynomial pieces on interval
    pub nbreak: size_t, // number of breakpoints (l + 1)
    pub n: size_t,      // number of bspline basis functions (l + k - 1)
    pub knots: *mut gsl_vector,  // knots vector
    pub deltal: *mut gsl_vector, // left delta
    pub deltar: *mut gsl_vector, // right delta
    pub B: *mut gsl_vector       // temporary spline results
}

#[repr(C)]
#[cfg(not(feature = "v2"))]
pub struct gsl_bspline_deriv_workspace {
    pub k: size_t,          // spline order
    pub A: *mut gsl_matrix, // work matrix
    pub dB: *mut gsl_matrix // temporary derivative results
}

/*#[repr(C)]
pub struct gsl_multifit_fsolver {
    pub type_: *const gsl_multifit_fsolver_type,
    pub function: *mut gsl_multifit_function,
    pub x: *mut gsl_vector,
    pub f: *mut gsl_vector,
    pub dx: *mut gsl_vector,
    pub state: *mut c_void,
}

#[repr(C)]
pub struct gsl_multifit_fsolver_type {
    pub name: *const c_char,
    pub size: size_t,
    pub alloc: Option<extern "C" fn(state: *mut c_void, n: size_t, p: size_t) -> c_int>,
    pub set: Option<extern "C" fn(state: *mut c_void, function: *mut gsl_multifit_function,
                                  x: *mut gsl_vector, f: *mut gsl_vector,
                                  dx: *mut gsl_vector) -> c_int>,
    pub iterate: Option<extern "C" fn(state: *mut c_void, function: *mut gsl_multifit_function,
                                      x: *mut gsl_vector, f: *mut gsl_vector,
                                      dx: *mut gsl_vector) -> c_int>,
    pub free: Option<extern "C" fn(state: *mut c_void)>,
}

#[repr(C)]
pub struct gsl_multifit_function {
    f: Option<extern "C" fn(x: *mut gsl_vector, params: *mut c_void, f: *mut gsl_vector) -> c_int>,
    n: size_t,
    p: size_t,
    params: *mut c_void,
}*/

pub type rng_set = Option<extern "C" fn(state: *mut c_void, seed: c_ulong)>;
pub type rng_get = Option<extern "C" fn(state: *mut c_void) -> c_ulong>;
pub type rng_get_double = Option<extern "C" fn(state: *mut c_void) -> c_double>;

#[repr(C)]
pub struct gsl_rng_type {
    pub name: *const c_char,
    pub max: c_ulong,
    pub min: c_ulong,
    pub size: size_t,
    pub set: rng_set,
    pub get: rng_get,
    pub get_double: rng_get_double
}

#[repr(C)]
pub struct gsl_rng {
    pub _type: *const gsl_rng_type,
    pub state: *mut c_void
}

#[repr(C)]
pub struct gsl_ran_discrete_t {
    pub K: size_t,
    pub A: *mut size_t,
    pub F: *mut c_double
}

#[repr(C)]
pub struct gsl_permutation {
    pub size: size_t,
    pub data: *mut size_t
}

#[repr(C)]
pub struct gsl_cheb_series {
    pub c: *mut c_double, // coefficients
    pub order: c_int, // order of expansion
    pub a: c_double, // lower interval point
    pub b: c_double, // upper interval point
    pub order_sp: size_t,
    pub f: *mut c_double
}

#[repr(C)]
pub struct gsl_combination {
    pub n: size_t,
    pub k: size_t,
    pub data: *mut size_t
}

#[repr(C)]
pub struct gsl_poly_complex_workspace {
    pub nc: size_t,
    pub matrix: *mut c_double
}

#[repr(C)]
pub struct gsl_dht {
    pub size: size_t, // size of the sample arrays to be transformed
    pub nu: c_double, // Bessel function order
    pub xmax: c_double, // the upper limit to the x-sampling domain
    pub kmax: c_double, // the upper limit to the k-sampling domain
    pub j: *mut c_double, // array of computed J_nu zeros, j_{nu,s} = j[s]
    pub Jjj: *mut c_double, // transform numerator, J_nu(j_i j_m / j_N)
    pub J2: *mut c_double // transform denominator, J_{nu+1}^2(j_m)
}

#[repr(C)]
pub struct gsl_eigen_symm_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double
}

#[repr(C)]
pub struct gsl_eigen_herm_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
    pub tau: *mut c_double
}

#[repr(C)]
pub struct gsl_eigen_hermv_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
    pub gc: *mut c_double,
    pub gs: *mut c_double
}

#[repr(C)]
pub struct gsl_eigen_francis_workspace {
    pub size: size_t, // matrix size
    pub max_iterations: size_t, // max iterations since last eigenvalue found
    pub n_iter: size_t, // number of iterations since last eigenvalue found
    pub n_evals: size_t, // number of eigenvalues found so far
    pub compute_t: c_int, // compute Schur form T = Z^t A Z
    pub H: *mut gsl_matrix, // pointer to Hessenberg matrix
    pub Z: *mut gsl_matrix // pointer to Schur vector matrix
}

#[repr(C)]
pub struct gsl_eigen_nonsymm_workspace {
    pub size: size_t, // size of matrices
    pub diag: *mut gsl_vector, // diagonal matrix elements from balancing
    pub tau: *mut gsl_vector, // Householder coefficients
    pub Z: *mut gsl_matrix, // pointer to Z matrix
    pub do_balance: c_int, // perform balancing transformation?
    pub n_evals: size_t, // number of eigenvalues found
    pub francis_workspace_p: *mut gsl_eigen_francis_workspace
}

#[repr(C)]
pub struct gsl_eigen_nonsymmv_workspace {
    pub size: size_t, // size of matrices
    pub work: *mut gsl_vector, // scratch workspace
    pub work2: *mut gsl_vector, // scratch workspace
    pub work3: *mut gsl_vector, // scratch workspace
    pub Z: *mut gsl_matrix, // pointer to Schur vectors
    pub nonsymm_workspace_p: *mut gsl_eigen_nonsymm_workspace
}

#[repr(C)]
pub struct gsl_eigen_gensymm_workspace {
    pub size: size_t,
    pub symm_workspace_p: gsl_eigen_symm_workspace
}

#[repr(C)]
pub struct gsl_eigen_gensymmv_workspace {
    pub size: size_t,
    pub symmv_workspace_p: gsl_eigen_symmv_workspace
}

#[repr(C)]
pub struct gsl_eigen_genherm_workspace {
    pub size: size_t,
    pub herm_workspace_p: *mut gsl_eigen_herm_workspace
}

#[repr(C)]
pub struct gsl_eigen_genhermv_workspace {
    pub size: size_t,
    pub hermv_workspace_p: *mut gsl_eigen_hermv_workspace
}

#[repr(C)]
pub struct gsl_eigen_gen_workspace {
    pub size: size_t, // size of matrices
    pub work: *mut gsl_vector, // scratch workspace
    pub n_evals: size_t, // number of eigenvalues found
    pub max_iterations: size_t, // maximum QZ iterations allowed
    pub n_iter: size_t, // number of iterations since last eigenvalue found
    pub eshift: c_double, // exceptional shift counter
    pub need_top: c_int, // need to compute top index?
    pub atol: c_double, // tolerance for splitting A matrix
    pub btol: c_double, // tolerance for splitting B matrix
    pub ascale: c_double, // scaling factor for shifts
    pub bscale: c_double, // scaling factor for shifts
    pub H: *mut gsl_matrix, // pointer to hessenberg matrix
    pub R: *mut gsl_matrix, // pointer to upper triangular matrix
    pub compute_s: c_int, // compute generalized Schur form S
    pub compute_t: c_int, // compute generalized Schur form T
    pub Q: *mut gsl_matrix, // pointer to left Schur vectors
    pub Z: *mut gsl_matrix // pointer to right Schur vectors
}

#[repr(C)]
pub struct gsl_eigen_genv_workspace {
    pub size: size_t, // size of matrices
    pub work1: *mut gsl_vector, // 1-norm of columns of A
    pub work2: *mut gsl_vector, // 1-norm of columns of B
    pub work3: *mut gsl_vector, // real part of eigenvector
    pub work4: *mut gsl_vector, // imag part of eigenvector
    pub work5: *mut gsl_vector, // real part of back-transformed eigenvector
    pub work6: *mut gsl_vector, // imag part of back-transformed eigenvector
    pub Q: *mut gsl_matrix, // pointer to left Schur vectors
    pub Z: *mut gsl_matrix, // pointer to right Schur vectors
    pub gen_workspace_p: *mut gsl_eigen_gen_workspace
}

#[repr(C)]
pub struct gsl_fft_complex_wavetable {
    pub n: size_t,
    pub nf: size_t,
    pub factor: [size_t; 64],
    pub twiddle: [*mut gsl_complex; 64],
    pub trig: *mut gsl_complex
}

#[repr(C)]
pub struct gsl_fft_complex_workspace {
    pub n: size_t,
    pub scratch: *mut c_double
}

#[repr(C)]
pub struct gsl_histogram {
    pub n: size_t, // This is the number of histogram bins
    pub range: *mut c_double, // The ranges of the bins are stored in an array of n+1 elements pointed to by range.
    pub bin: *mut c_double /* The counts for each bin are stored in an array of n elements pointed to by bin. The bins are floating-point numbers, so you can increment them by non-integer values if necessary.

The range for bin[i] is given by range[i] to range[i+1]. For n bins there are n+1 entries in the array range. Each bin is inclusive at the lower end and exclusive at the upper end. Mathematically this means that the bins are defined by the following inequality,

bin[i] corresponds to range[i] <= x < range[i+1]
Here is a diagram of the correspondence between ranges and bins on the number-line for x,

     [ bin[0] )[ bin[1] )[ bin[2] )[ bin[3] )[ bin[4] )
  ---|---------|---------|---------|---------|---------|---  x
   r[0]      r[1]      r[2]      r[3]      r[4]      r[5]

In this picture the values of the range array are denoted by r. On the left-hand side of each bin the square bracket ‘[’ denotes an inclusive lower bound (r <= x), and the round parentheses ‘)’ on the right-hand side denote an exclusive upper bound (x < r). Thus any samples which fall on the upper end of the histogram are excluded. If you want to include this value for the last bin you will need to add an extra bin to your histogram.
*/
}

#[repr(C)]
pub struct gsl_histogram_pdf {
    pub n: size_t, // This is the number of bins used to approximate the probability distribution function.
    pub range: *mut c_double, // The ranges of the bins are stored in an array of n+1 elements pointed to by range.
    pub sum: *mut c_double // The cumulative probability for the bins is stored in an array of n elements pointed to by sum.
}

#[repr(C)]
pub struct gsl_histogram2d {
    pub nx: size_t, // This is the number of histogram bins in the x direction.
    pub ny: size_t, // This is the number of histogram bins in the y direction.
    pub xrange: *mut c_double, // The ranges of the bins in the x-direction are stored in an array of nx + 1 elements pointed to by xrange.
    pub yrange: *mut c_double, // The ranges of the bins in the y-direction are stored in an array of ny + 1 elements pointed to by yrange.
    pub bin: *mut c_double, /*The counts for each bin are stored in an array pointed to by bin. The bins are floating-point numbers, so you can increment them by non-integer values if necessary. The array bin stores the two dimensional array of bins in a single block of memory according to the mapping bin(i,j) = bin[i * ny + j].

The range for bin(i,j) is given by xrange[i] to xrange[i+1] in the x-direction and yrange[j] to yrange[j+1] in the y-direction. Each bin is inclusive at the lower end and exclusive at the upper end. Mathematically this means that the bins are defined by the following inequality,

bin(i,j) corresponds to xrange[i] <= x < xrange[i+1]
                    and yrange[j] <= y < yrange[j+1]
Note that any samples which fall on the upper sides of the histogram are excluded. If you want to include these values for the side bins you will need to add an extra row or column to your histogram.
*/
}

#[repr(C)]
pub struct gsl_histogram2d_pdf {
    pub nx: size_t, // This is the number of histogram bins used to approximate the probability distribution function in the x direction.
    pub ny: size_t, // This is the number of histogram bins used to approximate the probability distribution function in the y direction.
    pub xrange: *mut c_double, // The ranges of the bins in the x-direction are stored in an array of nx + 1 elements pointed to by xrange.
    pub yrange: *mut c_double, // The ranges of the bins in the y-direction are stored in an array of ny + 1 pointed to by yrange.
    pub sum: *mut c_double // The cumulative probability for the bins is stored in an array of nx*ny elements pointed to by sum.
}

#[repr(C)]
pub struct gsl_integration_workspace {
    pub limit: size_t,
    pub size: size_t,
    pub nrmax: size_t,
    pub i: size_t,
    pub maximum_level: size_t,
    pub alist: *mut c_double,
    pub blist: *mut c_double,
    pub rlist: *mut c_double,
    pub elist: *mut c_double,
    pub order: *mut size_t,
    pub level: *mut size_t
}

#[repr(C)]
pub struct extrapolation_table {
    pub n: size_t,
    pub rlist2: [c_double; 52],
    pub nres: size_t,
    pub res3la: [c_double; 3]
}

#[repr(C)]
pub struct gsl_integration_qaws_table {
    pub alpha: c_double,
    pub beta: c_double,
    pub mu: c_int,
    pub nu: c_int,
    pub ri: [c_double; 25],
    pub rj: [c_double; 25],
    pub rg: [c_double; 25],
    pub rh: [c_double; 25]
}

#[repr(C)]
pub struct gsl_integration_qawo_table {
    pub n: size_t,
    pub omega: c_double,
    pub L: c_double,
    pub par: c_double,
    pub sine: ::IntegrationQawo,
    pub chebmo: *mut c_double
}

/* Data of a single interval */
#[repr(C)]
pub struct gsl_integration_cquad_ival {
    pub a: c_double,
    pub b: c_double,
    pub c: [c_double; 64],
    pub fx: [c_double; 33],
    pub igral: c_double,
    pub err: c_double,
    pub depth: c_int,
    pub rdepth: c_int,
    pub ndiv: c_int
}


/* The workspace is just a collection of intervals */
#[repr(C)]
pub struct gsl_integration_cquad_workspace {
    pub size: size_t,
    pub ivals: *mut gsl_integration_cquad_ival,
    pub heap: *mut size_t
}

/* Workspace for fixed-order Gauss-Legendre integration */
#[repr(C)]
pub struct gsl_integration_glfixed_table {
    pub n: size_t, /* number of points */
    pub x: *mut c_double, /* Gauss abscissae/points */
    pub w: *mut c_double, /* Gauss weights for each abscissae */
    pub precomputed: c_int /* high precision abscissae/weights precomputed? */
}

/* interpolation object type */
#[repr(C)]
pub struct gsl_interp_type {
    pub name: *const c_char,
    pub min_size: c_uint,
    pub alloc: Option<extern "C" fn(size_t) -> *mut c_void>,
    pub init: Option<extern "C" fn(*mut c_void, *const c_double, *const c_double, size_t) -> enums::Value>,
    pub eval: Option<extern "C" fn(*const c_void, *const c_double, *const c_double, size_t, c_double, *mut ::InterpAccel, *mut c_double) -> enums::Value>,
    pub eval_deriv: Option<extern "C" fn(*const c_void, *const c_double, *const c_double, size_t, c_double, *mut ::InterpAccel, *mut c_double) -> enums::Value>,
    pub eval_deriv2: Option<extern "C" fn(*const c_void, *const c_double, *const c_double, size_t, c_double, *mut ::InterpAccel, *mut c_double) -> enums::Value>,
    pub eval_integ: Option<extern "C" fn(*const c_void, *const c_double, *const c_double, size_t, c_double, *mut ::InterpAccel, c_double, c_double, *mut c_double) -> enums::Value>,
    pub free: Option<extern "C" fn(*mut c_void)>
}

/* general interpolation object */
#[repr(C)]
pub struct gsl_interp {
    pub _type: *const gsl_interp_type,
    pub xmin: c_double,
    pub xmax: c_double,
    pub size: size_t,
    pub state: *mut c_void
}

/* general interpolation object */
#[repr(C)]
pub struct gsl_spline {
    pub interp: *mut gsl_interp,
    pub x: *mut c_double,
    pub y: *mut c_double,
    pub size: size_t
}

/*#[repr(C)]
pub struct gsl_min_fminimizer_type {
    pub name: *const c_char,
    pub size: size_t,
    pub set: Option<extern "C" fn(state: *mut c_void, f: *mut gsl_function, x_minimum: c_double, f_minimum: c_double, x_lower: c_double,
        f_lower: c_double, x_upper: c_double, f_upper: c_double) -> enums::Value>,
    pub iterate: Option<extern "C" fn(state: *mut c_void, f: *mut gsl_function, x_minimum: *mut c_double, f_minimum: *mut c_double,
        x_lower: *mut c_double, f_lower: *mut c_double, x_upper: *mut c_double, f_upper: *mut c_double) -> enums::Value>
}

#[repr(C)]
pub struct gsl_min_fminimizer {
    pub type_: *const gsl_min_fminimizer_type,
    pub function: *mut gsl_function,
    pub x_minimum: c_double,
    pub x_lower: c_double,
    pub x_upper: c_double,
    pub f_minimum: c_double,
    pub f_lower: c_double,
    pub f_upper: c_double,
    pub state: *mut c_void
}

#[repr(C)]
pub struct gsl_function {
    pub function: Option<extern "C" fn(x: c_double, params: *mut c_void) -> c_double>,
    pub params: *mut c_void
}*/

#[repr(C)]
pub struct gsl_monte_plain_state {
    pub dim: size_t,
    pub x: *mut c_double
}

#[repr(C)]
pub struct gsl_monte_function {
    pub f: *mut c_void,
    pub dim: size_t,
    pub params: *mut c_void,
}

#[repr(C)]
pub struct gsl_monte_miser_state {
    pub min_calls: size_t,
    pub min_calls_per_bisection: size_t,
    pub dither: c_double,
    pub estimate_frac: c_double,
    pub alpha: c_double,
    pub dim: size_t,
    pub estimate_style: c_int,
    pub depth: c_int,
    pub verbose: c_int,
    pub x: *mut c_double,
    pub xmid: *mut c_double,
    pub sigma_l: *mut c_double,
    pub sigma_r: *mut c_double,
    pub fmax_l: *mut c_double,
    pub fmax_r: *mut c_double,
    pub fmin_l: *mut c_double,
    pub fmin_r: *mut c_double,
    pub fsum_l: *mut c_double,
    pub fsum_r: *mut c_double,
    pub fsum2_l: *mut c_double,
    pub fsum2_r: *mut c_double,
    pub hits_l: *mut size_t,
    pub hits_r: *mut size_t
}

#[repr(C)]
pub struct gsl_monte_vegas_state {
    /* grid */
    pub dim: size_t,
    pub bins_max: size_t,
    pub bins: c_uint,
    pub boxes: c_uint, /* these are both counted along the axes */
    pub xi: *mut c_double,
    pub xin: *mut c_double,
    pub delx: *mut c_double,
    pub weight: *mut c_double,
    pub vol: c_double,

    pub x: *mut c_double,
    pub bin: *mut c_int,
    pub box_: *mut c_int,

    /* distribution */
    pub d: *mut c_double,

    /* control variables */
    pub alpha: c_double,
    pub mode: ::VegasMode,
    pub verbose: c_int,
    pub iterations: c_uint,
    pub stage: c_int,

    /* scratch variables preserved between calls to vegas1/2/3  */
    pub jac: c_double,
    pub wtd_int_sum: c_double,
    pub sum_wgts: c_double,
    pub chi_sum: c_double,
    pub chisq: c_double,

    pub result: c_double,
    pub sigma: c_double,

    pub it_start: c_uint,
    pub it_num: c_uint,
    pub samples: c_uint,
    pub calls_per_box: c_uint,

    pub ostream: *mut c_void // is supposed to be a FILE*
}

/*#[repr(C)]
pub struct gsl_monte_vegas_params {
    pub alpha: c_double,
    pub iterations: size_t,
    pub stage: c_int,
    pub mode: ::VegasMode,
    pub verbose: c_int,
    pub ostream: *mut c_void // is supposed to be a FILE*
}*/

#[repr(C)]
pub struct gsl_ntuple {
    pub file: *mut FILE,
    pub ntuple_data: *mut c_void,
    pub size: size_t,
}

#[repr(C)]
pub struct gsl_multiset {
    pub n: size_t,
    pub k: size_t,
    pub data: *mut size_t
}

#[repr(C)]
pub struct gsl_odeiv2_system {
    pub function: extern fn(t: f64, *const f64, *mut f64, *mut c_void) -> enums::Value,
    pub jacobian: Option<extern fn(t: f64, *const f64, *mut f64, *mut f64, *mut c_void) -> enums::Value>,
    pub dimension: usize,
    pub params: *mut c_void
}

#[repr(C)]
pub struct gsl_odeiv2_driver {
    // ODE system
    pub sys: *const gsl_odeiv2_system,
    // stepper object
    pub s: *mut gsl_odeiv2_step,
    // control object
    pub c: *mut gsl_odeiv2_control,
    // evolve object
    pub e: *mut gsl_odeiv2_evolve,
    // step size
    pub h: c_double,
    // minimum step size allowed
    pub hmin: c_double,
    // maximum step size allowed
    pub hmax: c_double,
    // number of steps taken
    pub n: usize,
    // Maximum number of steps allowed
    pub nmax: usize
}

#[repr(C)]
pub struct gsl_odeiv2_evolve {
    pub dimension: size_t,
    pub y0: *mut c_double,
    pub yerr: *mut c_double,
    pub dydt_in: *mut c_double,
    pub dydt_out: *mut c_double,
    pub last_step: c_double,
    pub count: usize,
    pub failed_steps: usize,
    pub driver: *const gsl_odeiv2_driver
}

#[repr(C)]
pub struct gsl_odeiv2_control {
    pub type_: *const gsl_odeiv2_control_type,
    pub state: *mut c_void
}

#[repr(C)]
pub struct gsl_odeiv2_control_type {
    pub name: *const c_char,
    pub alloc: fn() -> *mut c_void,
    pub init: fn(state: *mut c_void, eps_abs: c_double, eps_rel: c_double, a_y: c_double, a_dydt: c_double) -> enums::Value,
    pub hadjust: fn(state: *mut c_void, dim: size_t, ord: c_uint, y: *const c_double, yerr: *const c_double, yp: *const c_double,
        h: *mut c_double) -> enums::Value,
    pub errlevel: fn(state: *mut c_void, y: c_double, dydt: c_double, h: c_double, ind: size_t, errlev: *mut c_double) -> enums::Value,
    pub set_driver: fn(state: *mut c_void, d: *const gsl_odeiv2_driver) -> enums::Value,
    pub free: fn(state: *mut c_void)
}

#[repr(C)]
pub struct gsl_odeiv2_step {
    pub type_: *const gsl_odeiv2_step_type,
    pub dimension: size_t,
    pub state: *mut c_void
}

#[repr(C)]
pub struct gsl_odeiv2_step_type {
    pub name: *const c_char,
    pub can_use_dydt_in: c_int,
    pub gives_exact_dydt_out: c_int,
    pub alloc: fn(dim: size_t) -> *mut c_void,
    pub apply: fn(state: *mut c_void, dim: size_t, t: c_double, h: c_double, y: *mut c_double, yerr: *mut c_double, dydt_in: *const c_double,
        dydt_out: *mut c_double, dydt: *const gsl_odeiv2_system) -> enums::Value,
    pub set_driver: fn(state: *mut c_void, d: *const gsl_odeiv2_driver) -> enums::Value,
    pub reset: fn(state: *mut c_void, dim: size_t) -> enums::Value,
    pub order: fn(state: *mut c_void) -> c_uint,
    pub free: fn(state: *mut c_void)
}

// Structure describing a generator instance of a specified type, with generator-specific state info and dimension-specific info.
#[repr(C)]
pub struct gsl_qrng {
    pub type_: *const gsl_qrng_type,
    pub dimension: c_uint,
    pub state_size: size_t,
    pub state: *mut c_void
}

// Structure describing a type of generator.
#[repr(C)]
pub struct gsl_qrng_type {
    pub name: *const c_char,
    pub max_dimension: c_uint,
    pub state_size: Option<extern "C" fn(dimension: c_uint) -> size_t>,
    pub init_state: Option<extern "C" fn(state: *mut c_void, dimension: c_uint) -> enums::Value>,
    pub get: Option<extern "C" fn(state: *mut c_void, dimension: c_uint, x: *mut c_double) -> enums::Value>
}

/*
 * size        = number of terms the workspace can handle
 * sum_plain   = simple sum of series
 * q_num       = backward diagonal of numerator; length = size
 * q_den       = backward diagonal of denominator; length = size
 * dq_num      = table of numerator derivatives; length = size**2
 * dq_den      = table of denominator derivatives; length = size**2
 * dsum        = derivative of sum wrt term i; length = size
*/
#[repr(C)]
pub struct gsl_sum_levin_u_workspace {
    pub size: size_t,
    // position in array
    pub i: size_t,
    // number of calls
    pub terms_used: size_t,
    pub sum_plain: c_double,
    pub q_num: *mut c_double,
    pub q_den: *mut c_double,
    pub dq_num: *mut c_double,
    pub dq_den: *mut c_double,
    pub dsum: *mut c_double
}

#[repr(C)]
pub struct gsl_sum_levin_utrunc_workspace {
    pub size: size_t,
    // position in array
    pub i: size_t,
    // number of calls
    pub terms_used: size_t,
    pub sum_plain: c_double,
    pub q_num: *mut c_double,
    pub q_den: *mut c_double,
    pub dsum: *mut c_double
}

#[repr(C)]
pub struct gsl_wavelet_workspace {
    pub scratch: *mut c_double,
    pub n: size_t
}

#[repr(C)]
pub struct gsl_wavelet {
    pub type_: *const gsl_wavelet_type,
    pub h1: *const c_double,
    pub g1: *const c_double,
    pub h2: *const c_double,
    pub g2: *const c_double,
    pub nc: size_t,
    pub offset: size_t
}

#[repr(C)]
pub struct gsl_wavelet_type {
    pub name: *const c_char,
    pub init: Option<extern "C" fn(h1: *const *const c_double, g1: *const *const c_double, h2: *const *const c_double,
        g2: *const *const c_double, nc: *mut size_t, offset: *mut size_t, member: size_t) -> enums::Value>
}
