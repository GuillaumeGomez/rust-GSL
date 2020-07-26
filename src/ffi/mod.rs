//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![allow(improper_ctypes)]

use libc::{c_double, c_int, c_uint, c_float, c_void, size_t, c_ulong, c_char, FILE};

mod blas;
mod linalg;
mod monte_carlo;
mod randist;
mod solvers;

pub use self::blas::*;
pub use self::linalg::*;
pub use self::monte_carlo::*;
pub use self::randist::*;
pub use self::solvers::*;

pub type gsl_complex_packed_ptr = *mut c_double;
pub type gsl_complex_packed_array = *mut c_double;
#[allow(dead_code)]
pub type coord = c_int;

pub trait FFI<T> {
    fn wrap(r: *mut T) -> Self;
    fn soft_wrap(r: *mut T) -> Self;
    fn unwrap_shared(&Self) -> *const T;
    fn unwrap_unique(&mut Self) -> *mut T;
}

extern "C" {
    pub static gsl_rng_mt19937: *const gsl_rng_type;
    pub static gsl_rng_ranlxs0: *const gsl_rng_type;
    pub static gsl_rng_ranlxs1: *const gsl_rng_type;
    pub static gsl_rng_ranlxs2: *const gsl_rng_type;
    pub static gsl_rng_ranlxd1: *const gsl_rng_type;
    pub static gsl_rng_ranlxd2: *const gsl_rng_type;
    pub static gsl_rng_ranlux: *const gsl_rng_type;
    pub static gsl_rng_ranlux389: *const gsl_rng_type;
    pub static gsl_rng_cmrg: *const gsl_rng_type;
    pub static gsl_rng_mrg: *const gsl_rng_type;
    pub static gsl_rng_taus: *const gsl_rng_type;
    pub static gsl_rng_taus2: *const gsl_rng_type;
    pub static gsl_rng_gfsr4: *const gsl_rng_type;

    pub static gsl_rng_rand: *const gsl_rng_type;
    pub static gsl_rng_random_bsd: *const gsl_rng_type;
    pub static gsl_rng_random_libc5: *const gsl_rng_type;
    pub static gsl_rng_random_glibc2: *const gsl_rng_type;
    pub static gsl_rng_rand48: *const gsl_rng_type;

    pub static gsl_rng_default: *const gsl_rng_type;
    pub static gsl_rng_ranf: *const gsl_rng_type;
    pub static gsl_rng_ranmar: *const gsl_rng_type;
    pub static gsl_rng_r250: *const gsl_rng_type;
    pub static gsl_rng_tt800: *const gsl_rng_type;
    pub static gsl_rng_vax: *const gsl_rng_type;
    pub static gsl_rng_transputer: *const gsl_rng_type;
    pub static gsl_rng_randu: *const gsl_rng_type;
    pub static gsl_rng_minstd: *const gsl_rng_type;
    pub static gsl_rng_uni: *const gsl_rng_type;
    pub static gsl_rng_uni32: *const gsl_rng_type;
    pub static gsl_rng_slatec: *const gsl_rng_type;
    pub static gsl_rng_zuf: *const gsl_rng_type;
    pub static gsl_rng_knuthran2: *const gsl_rng_type;
    pub static gsl_rng_knuthran2002: *const gsl_rng_type;
    pub static gsl_rng_knuthran: *const gsl_rng_type;
    pub static gsl_rng_borosh13: *const gsl_rng_type;
    pub static gsl_rng_fishman18: *const gsl_rng_type;
    pub static gsl_rng_fishman20: *const gsl_rng_type;
    pub static gsl_rng_lecuyer21: *const gsl_rng_type;
    pub static gsl_rng_waterman14: *const gsl_rng_type;
    pub static gsl_rng_fishman2x: *const gsl_rng_type;
    pub static gsl_rng_coveyou: *const gsl_rng_type;

    pub static gsl_rng_default_seed: c_ulong;

    pub static gsl_interp_linear: *const gsl_interp_type;
    pub static gsl_interp_polynomial: *const gsl_interp_type;
    pub static gsl_interp_cspline: *const gsl_interp_type;
    pub static gsl_interp_cspline_periodic: *const gsl_interp_type;
    pub static gsl_interp_akima: *const gsl_interp_type;
    pub static gsl_interp_akima_periodic: *const gsl_interp_type;

    pub static gsl_odeiv2_step_rk2: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk4: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rkf45: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rkck: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk8pd: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk1imp: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk2imp: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_rk4imp: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_bsimp: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_msadams: *const gsl_odeiv2_step_type;
    pub static gsl_odeiv2_step_msbdf: *const gsl_odeiv2_step_type;

    pub static gsl_odeiv2_control_scaled: *const gsl_odeiv2_control_type;
    pub static gsl_odeiv2_control_standard: *const gsl_odeiv2_control_type;

    pub static gsl_qrng_niederreiter_2: *const gsl_qrng_type;
    pub static gsl_qrng_sobol: *const gsl_qrng_type;
    pub static gsl_qrng_halton: *const gsl_qrng_type;
    pub static gsl_qrng_reversehalton: *const gsl_qrng_type;

    pub static gsl_wavelet_daubechies: *const gsl_wavelet_type;
    pub static gsl_wavelet_daubechies_centered: *const gsl_wavelet_type;
    pub static gsl_wavelet_haar: *const gsl_wavelet_type;
    pub static gsl_wavelet_haar_centered: *const gsl_wavelet_type;
    pub static gsl_wavelet_bspline: *const gsl_wavelet_type;
    pub static gsl_wavelet_bspline_centered: *const gsl_wavelet_type;

    // Airy functions
    pub fn gsl_sf_airy_Ai(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Bi(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_e(x: c_double, mode: ::Mode, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_airy_Ai_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_scaled_e(x: c_double,
                                   mode: ::Mode,
                                   result: *mut gsl_sf_result)
                                   -> c_int;
    pub fn gsl_sf_airy_Bi_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_scaled_e(x: c_double,
                                   mode: ::Mode,
                                   result: *mut gsl_sf_result)
                                   -> c_int;
    // Derivatives of Airy Functions
    pub fn gsl_sf_airy_Ai_deriv(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_e(x: c_double,
                                  mode: ::Mode,
                                  result: *mut gsl_sf_result)
                                  -> c_int;
    pub fn gsl_sf_airy_Bi_deriv(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_e(x: c_double,
                                  mode: ::Mode,
                                  result: *mut gsl_sf_result)
                                  -> c_int;
    pub fn gsl_sf_airy_Ai_deriv_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_scaled_e(x: c_double,
                                         mode: ::Mode,
                                         result: *mut gsl_sf_result)
                                         -> c_int;
    pub fn gsl_sf_airy_Bi_deriv_scaled(x: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_scaled_e(x: c_double,
                                         mode: ::Mode,
                                         result: *mut gsl_sf_result)
                                         -> c_int;
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
    pub fn gsl_sf_bessel_In_array(nmin: c_int,
                                  nmax: c_int,
                                  x: c_double,
                                  result_array: *mut c_double)
                                  -> c_int;
    pub fn gsl_sf_bessel_I0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_I1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_I1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_In_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_In_scaled_e(n: c_int,
                                     x: c_double,
                                     result: *mut gsl_sf_result)
                                     -> c_int;
    pub fn gsl_sf_bessel_In_scaled_array(nmin: c_int,
                                         nmax: c_int,
                                         x: c_double,
                                         result_array: *mut c_double)
                                         -> c_int;
    // Regular Modified Spherical Bessel Functions
    // The regular modified spherical Bessel functions i_l(x) are related to the modified Bessel functions of fractional order, i_l(x) = \sqrt{\pi/(2x)} I_{l+1/2}(x)
    pub fn gsl_sf_bessel_i0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_i1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_i2_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_i2_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_il_scaled(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_il_scaled_e(l: c_int,
                                     x: c_double,
                                     result: *mut gsl_sf_result)
                                     -> c_int;
    pub fn gsl_sf_bessel_il_scaled_array(lmax: c_int,
                                         x: c_double,
                                         result_array: *mut c_double)
                                         -> c_int;
    pub fn gsl_sf_bessel_Inu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_e(nu: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_bessel_Inu_scaled(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Inu_scaled_e(nu: c_double,
                                      x: c_double,
                                      result: *mut gsl_sf_result)
                                      -> c_int;
    // Regular Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_J0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_J1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_J1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Jn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Jn_array(nmin: c_int,
                                  nmax: c_int,
                                  x: c_double,
                                  result_array: *mut c_double)
                                  -> c_int;
    // Regular Spherical Bessel Functions
    pub fn gsl_sf_bessel_j0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_j1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_j2(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_j2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_jl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_jl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_jl_array(lmax: c_int,
                                  x: c_double,
                                  result_array: *mut c_double)
                                  -> c_int;
    pub fn gsl_sf_bessel_jl_steed_array(lmax: c_int,
                                        x: c_double,
                                        result_array: *mut c_double)
                                        -> c_int;
    // Regular Bessel Function—Fractional Order
    pub fn gsl_sf_bessel_Jnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Jnu_e(nu: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_bessel_sequence_Jnu_e(nu: c_double,
                                        mode: ::Mode,
                                        size: i64,
                                        v: *mut c_double)
                                        -> c_int;
    // Irregular Modified Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_K0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_K1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_array(nmin: c_int,
                                  nmax: c_int,
                                  x: c_double,
                                  result_array: *mut c_double)
                                  -> c_int;
    pub fn gsl_sf_bessel_K0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_K1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_K1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Kn_scaled(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Kn_scaled_e(n: c_int,
                                     x: c_double,
                                     result: *mut gsl_sf_result)
                                     -> c_int;
    pub fn gsl_sf_bessel_Kn_scaled_array(nmin: c_int,
                                         nmax: c_int,
                                         x: c_double,
                                         result_array: *mut c_double)
                                         -> c_int;
    // Irregular Modified Spherical Bessel Functions
    pub fn gsl_sf_bessel_k0_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k0_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_k1_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k1_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_k2_scaled(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_k2_scaled_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_kl_scaled(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_kl_scaled_e(l: c_int,
                                     x: c_double,
                                     result: *mut gsl_sf_result)
                                     -> c_int;
    pub fn gsl_sf_bessel_kl_scaled_array(lmax: c_int,
                                         x: c_double,
                                         result_array: *mut c_double)
                                         -> c_int;
    // Irregular Modified Bessel Functions—Fractional Order
    pub fn gsl_sf_bessel_Knu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Knu_e(nu: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_bessel_lnKnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_lnKnu_e(nu: c_double,
                                 x: c_double,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    pub fn gsl_sf_bessel_Knu_scaled(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Knu_scaled_e(nu: c_double,
                                      x: c_double,
                                      result: *mut gsl_sf_result)
                                      -> c_int;
    // Irregular Cylindrical Bessel Functions
    pub fn gsl_sf_bessel_Y0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Y0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Y1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Y1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Yn(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Yn_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_Yn_array(nmin: c_int,
                                  nmax: c_int,
                                  x: c_double,
                                  result_array: *mut c_double)
                                  -> c_int;
    // Irregular Spherical Bessel Functions
    pub fn gsl_sf_bessel_y0(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_y1(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_y2(x: c_double) -> c_double;
    pub fn gsl_sf_bessel_y2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_yl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_yl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_yl_array(lmax: c_int,
                                  x: c_double,
                                  result_array: *mut c_double)
                                  -> c_int;
    // Irregular Bessel Functions—Fractional Order
    pub fn gsl_sf_bessel_Ynu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_Ynu_e(nu: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    // Zeros of Regular Bessel Functions
    pub fn gsl_sf_bessel_zero_J0(s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_J0_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_zero_J1(s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_J1_e(s: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_bessel_zero_Jnu(nu: c_double, s: c_uint) -> c_double;
    pub fn gsl_sf_bessel_zero_Jnu_e(nu: c_double,
                                    s: c_uint,
                                    result: *mut gsl_sf_result)
                                    -> c_int;

    // Trigonometric Functions
    pub fn gsl_sf_sin(x: c_double) -> c_double;
    pub fn gsl_sf_sin_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_cos(x: c_double) -> c_double;
    pub fn gsl_sf_cos_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_hypot(x: c_double) -> c_double;
    pub fn gsl_sf_hypot_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_sinc(x: c_double) -> c_double;
    pub fn gsl_sf_sinc_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_complex_sin_e(zr: c_double,
                                zi: c_double,
                                szr: *mut gsl_sf_result,
                                szi: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_complex_cos_e(zr: c_double,
                                zi: c_double,
                                czr: *mut gsl_sf_result,
                                czi: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_complex_logsin_e(zr: c_double,
                                   zi: c_double,
                                   lszr: *mut gsl_sf_result,
                                   lszi: *mut gsl_sf_result)
                                   -> c_int;
    pub fn gsl_sf_lnsinh(x: c_double) -> c_double;
    pub fn gsl_sf_lnsinh_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lncosh(x: c_double) -> c_double;
    pub fn gsl_sf_lncosh_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_polar_to_rect(r: c_double,
                                theta: c_double,
                                x: *mut gsl_sf_result,
                                y: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_rect_to_polar(x: c_double,
                                y: c_double,
                                r: *mut gsl_sf_result,
                                theta: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_angle_restrict_symm(theta: c_double) -> c_double;
    pub fn gsl_sf_angle_restrict_symm_e(theta: *mut c_double) -> c_int;
    pub fn gsl_sf_angle_restrict_pos(theta: c_double) -> c_double;
    pub fn gsl_sf_angle_restrict_pos_e(theta: *mut c_double) -> c_int;
    pub fn gsl_sf_sin_err_e(x: c_double, dx: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_cos_err_e(x: c_double, dx: c_double, result: *mut gsl_sf_result) -> c_int;

    // Exponential Integrals functions
    pub fn gsl_sf_expint_E1(x: c_double) -> c_double;
    pub fn gsl_sf_expint_E1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_expint_E2(x: c_double) -> c_double;
    pub fn gsl_sf_expint_E2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_expint_En(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_expint_En_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_expint_Ei(x: c_double) -> c_double;
    pub fn gsl_sf_expint_Ei_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_Shi(x: c_double) -> c_double;
    pub fn gsl_sf_Shi_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_Chi(x: c_double) -> c_double;
    pub fn gsl_sf_Chi_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_expint_3(x: c_double) -> c_double;
    pub fn gsl_sf_expint_3_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_Si(x: c_double) -> c_double;
    pub fn gsl_sf_Si_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_Ci(x: c_double) -> c_double;
    pub fn gsl_sf_Ci_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_atanint(x: c_double) -> c_double;
    pub fn gsl_sf_atanint_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Clausen functions
    pub fn gsl_sf_clausen(x: c_double) -> c_double;
    pub fn gsl_sf_clausen_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Coulomb functions
    // Normalized Hydrogenic Bound States
    pub fn gsl_sf_hydrogenicR_1(Z: c_double, r: c_double) -> c_double;
    pub fn gsl_sf_hydrogenicR_1_e(Z: c_double,
                                  r: c_double,
                                  result: *mut gsl_sf_result)
                                  -> c_int;
    pub fn gsl_sf_hydrogenicR(n: c_int, l: c_int, Z: c_double, r: c_double) -> c_double;
    pub fn gsl_sf_hydrogenicR_e(n: c_int,
                                l: c_int,
                                Z: c_double,
                                r: c_double,
                                result: *mut gsl_sf_result)
                                -> c_int;
    // Coulomb Wave Functions
    // The Coulomb wave functions F_L(\eta,x), G_L(\eta,x) are described in Abramowitz & Stegun, Chapter 14. Because there can be a large dynamic range of values for these functions, overflows are handled gracefully. If an overflow occurs, GSL_EOVRFLW is signalled and exponent(s) are returned through the modifiable parameters exp_F, exp_G. The full solution can be reconstructed from the following relations,
    //
    // F_L(eta,x)  =  fc[k_L] * exp(exp_F)
    // G_L(eta,x)  =  gc[k_L] * exp(exp_G)
    //
    // F_L'(eta,x) = fcp[k_L] * exp(exp_F)
    // G_L'(eta,x) = gcp[k_L] * exp(exp_G)
    pub fn gsl_sf_coulomb_wave_FG_e(eta: c_double,
                                    x: c_double,
                                    L_F: c_double,
                                    k: c_int,
                                    F: *mut gsl_sf_result,
                                    Fp: *mut gsl_sf_result,
                                    G: *mut gsl_sf_result,
                                    Gp: *mut gsl_sf_result,
                                    exp_F: *mut c_double,
                                    exp_G: *mut c_double)
                                    -> c_int;
    pub fn gsl_sf_coulomb_wave_F_array(L_min: c_double,
                                       kmax: c_int,
                                       eta: c_double,
                                       x: c_double,
                                       fc_array: *mut c_double,
                                       F_exponent: *mut c_double)
                                       -> c_int;
    pub fn gsl_sf_coulomb_wave_FG_array(L_min: c_double,
                                        kmax: c_int,
                                        eta: c_double,
                                        x: c_double,
                                        fc_array: *mut c_double,
                                        gc_array: *mut c_double,
                                        F_exponent: *mut c_double,
                                        G_exponent: *mut c_double)
                                        -> c_int;
    pub fn gsl_sf_coulomb_wave_FGp_array(L_min: c_double,
                                         kmax: c_int,
                                         eta: c_double,
                                         x: c_double,
                                         fc_array: *mut c_double,
                                         fcp_array: *mut c_double,
                                         gc_array: *mut c_double,
                                         gcp_array: *mut c_double,
                                         F_exponent: *mut c_double,
                                         G_exponent: *mut c_double)
                                         -> c_int;
    pub fn gsl_sf_coulomb_wave_sphF_array(L_min: c_double,
                                          kmax: c_int,
                                          eta: c_double,
                                          x: c_double,
                                          fc_array: *mut c_double,
                                          f_exponent: *mut c_double)
                                          -> c_int;
    // Coulomb Wave Function Normalization Constant
    pub fn gsl_sf_coulomb_CL_e(L: c_double,
                               eta: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_coulomb_CL_array(Lmin: c_double,
                                   kmax: c_int,
                                   eta: c_double,
                                   cl: *mut c_double)
                                   -> c_int;

    // Coupling Coefficients functions
    pub fn gsl_sf_coupling_3j(two_ja: c_int,
                              two_jb: c_int,
                              two_jc: c_int,
                              two_ma: c_int,
                              two_mc: c_int,
                              two_mc: c_int)
                              -> c_double;
    pub fn gsl_sf_coupling_3j_e(two_ja: c_int,
                                two_jb: c_int,
                                two_jc: c_int,
                                two_ma: c_int,
                                two_mc: c_int,
                                two_mc: c_int,
                                result: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_coupling_6j(two_ja: c_int,
                              two_jb: c_int,
                              two_jc: c_int,
                              two_jd: c_int,
                              two_je: c_int,
                              two_jf: c_int)
                              -> c_double;
    pub fn gsl_sf_coupling_6j_e(two_ja: c_int,
                                two_jb: c_int,
                                two_jc: c_int,
                                two_jd: c_int,
                                two_je: c_int,
                                two_jf: c_int,
                                result: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_coupling_9j(two_ja: c_int,
                              two_jb: c_int,
                              two_jc: c_int,
                              two_jd: c_int,
                              two_je: c_int,
                              two_jf: c_int,
                              two_jg: c_int,
                              two_jh: c_int,
                              two_ji: c_int)
                              -> c_double;
    pub fn gsl_sf_coupling_9j_e(two_ja: c_int,
                                two_jb: c_int,
                                two_jc: c_int,
                                two_jd: c_int,
                                two_je: c_int,
                                two_jf: c_int,
                                two_jg: c_int,
                                two_jh: c_int,
                                two_ji: c_int,
                                result: *mut gsl_sf_result)
                                -> c_int;

    // Dawson functions
    pub fn gsl_sf_dawson(x: c_double) -> c_double;
    pub fn gsl_sf_dawson_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Debye functions
    pub fn gsl_sf_debye_1(x: c_double) -> c_double;
    pub fn gsl_sf_debye_1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_debye_2(x: c_double) -> c_double;
    pub fn gsl_sf_debye_2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_debye_3(x: c_double) -> c_double;
    pub fn gsl_sf_debye_3_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_debye_4(x: c_double) -> c_double;
    pub fn gsl_sf_debye_4_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_debye_5(x: c_double) -> c_double;
    pub fn gsl_sf_debye_5_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_debye_6(x: c_double) -> c_double;
    pub fn gsl_sf_debye_6_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Dilogarithm functions
    // real argument
    pub fn gsl_sf_dilog(x: c_double) -> c_double;
    pub fn gsl_sf_dilog_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    // complex argument
    pub fn gsl_sf_complex_dilog_e(r: c_double,
                                  theta: c_double,
                                  result: *mut gsl_sf_result,
                                  result_im: *mut gsl_sf_result)
                                  -> c_int;

    // Elementary Operations functions
    pub fn gsl_sf_multiply_e(x: c_double, y: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_multiply_err_e(x: c_double,
                                 dx: c_double,
                                 y: c_double,
                                 dy: c_double,
                                 result: *mut gsl_sf_result)
                                 -> c_int;

    // Elliptic functions (Jacobi)
    pub fn gsl_sf_elljac_e(u: c_double,
                           m: c_double,
                           sn: *mut c_double,
                           cn: *mut c_double,
                           dn: *mut c_double)
                           -> c_int;

    // Error functions
    pub fn gsl_sf_erf(x: c_double) -> c_double;
    pub fn gsl_sf_erf_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_set_error_handler(x: Option<extern "C" fn(*const c_char, *const c_char, c_int, c_int)>);
    pub fn gsl_set_error_handler_off();
    // Complementary Error functions
    pub fn gsl_sf_erfc(x: c_double) -> c_double;
    pub fn gsl_sf_erfc_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Log Complementary Error functions
    pub fn gsl_sf_log_erfc(x: c_double) -> c_double;
    pub fn gsl_sf_log_erfc_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Probability functions
    // The probability functions for the Normal or Gaussian distribution are described in Abramowitz & Stegun, Section 26.2.
    pub fn gsl_sf_erf_Z(x: c_double) -> c_double;
    pub fn gsl_sf_erf_Z_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_erf_Q(x: c_double) -> c_double;
    pub fn gsl_sf_erf_Q_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_hazard(x: c_double) -> c_double;
    pub fn gsl_sf_hazard_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Exponential functions
    pub fn gsl_sf_exp(x: c_double) -> c_double;
    pub fn gsl_sf_exp_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_exp_e10_e(x: c_double, result: *mut gsl_sf_result_e10) -> c_int;
    pub fn gsl_sf_exp_mult(x: c_double, y: c_double) -> c_double;
    pub fn gsl_sf_exp_mult_e(x: c_double, y: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_exp_mult_e10_e(x: c_double,
                                 y: c_double,
                                 result: *mut gsl_sf_result_e10)
                                 -> c_int;
    // Relative Exponential functions
    pub fn gsl_sf_expm1(x: c_double) -> c_double;
    pub fn gsl_sf_expm1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_exprel(x: c_double) -> c_double;
    pub fn gsl_sf_exprel_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_exprel_2(x: c_double) -> c_double;
    pub fn gsl_sf_exprel_2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_exprel_n(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_exprel_n_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Exponentiation With Error Estimate
    pub fn gsl_sf_exp_err_e(x: c_double, dx: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_exp_err_e10_e(x: c_double,
                                dx: c_double,
                                result: *mut gsl_sf_result_e10)
                                -> c_int;
    pub fn gsl_sf_exp_mult_err_e(x: c_double,
                                 dx: c_double,
                                 y: c_double,
                                 dy: c_double,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    pub fn gsl_sf_exp_mult_err_e10_e(x: c_double,
                                     dx: c_double,
                                     y: c_double,
                                     dy: c_double,
                                     result: *mut gsl_sf_result_e10)
                                     -> c_int;

    // Gamma Beta functions
    // Gamma functions
    pub fn gsl_sf_gamma(x: c_double) -> c_double;
    pub fn gsl_sf_gamma_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lngamma(x: c_double) -> c_double;
    pub fn gsl_sf_lngamma_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lngamma_sgn_e(x: c_double,
                                result_lg: *mut gsl_sf_result,
                                sgn: *mut c_double)
                                -> c_int;
    pub fn gsl_sf_gammastar(x: c_double) -> c_double;
    pub fn gsl_sf_gammastar_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_gammainv(x: c_double) -> c_double;
    pub fn gsl_sf_gammainv_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lngamma_complex_e(zr: c_double,
                                    zi: c_double,
                                    lnr: *mut gsl_sf_result,
                                    arg: *mut gsl_sf_result)
                                    -> c_int;
    // Factorials
    pub fn gsl_sf_fact(n: c_uint) -> c_double;
    pub fn gsl_sf_fact_e(n: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_doublefact(n: c_uint) -> c_double;
    pub fn gsl_sf_doublefact_e(n: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lnfact(n: c_uint) -> c_double;
    pub fn gsl_sf_lnfact_e(n: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lndoublefact(n: c_uint) -> c_double;
    pub fn gsl_sf_lndoublefact_e(n: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_choose(n: c_uint, m: c_uint) -> c_double;
    pub fn gsl_sf_choose_e(n: c_uint, m: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lnchoose(n: c_uint, m: c_uint) -> c_double;
    pub fn gsl_sf_lnchoose_e(n: c_uint, m: c_uint, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_taylorcoeff(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_taylorcoeff_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Pochhammer Symbol
    pub fn gsl_sf_poch(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_poch_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lnpoch(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_lnpoch_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lnpoch_sgn_e(a: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result,
                               sgn: *mut c_double)
                               -> c_int;
    pub fn gsl_sf_pochrel(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_pochrel_e(a: c_double, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Beta functions
    pub fn gsl_sf_beta(a: c_double, b: c_double) -> c_double;
    pub fn gsl_sf_beta_e(a: c_double, b: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lnbeta(a: c_double, b: c_double) -> c_double;
    pub fn gsl_sf_lnbeta_e(a: c_double, b: c_double, result: *mut gsl_sf_result) -> c_int;
    // Incomplete Gamma functions
    pub fn gsl_sf_gamma_inc(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gamma_inc_e(a: c_double,
                              x: c_double,
                              result: *mut gsl_sf_result)
                              -> c_int;
    pub fn gsl_sf_gamma_inc_Q(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gamma_inc_Q_e(a: c_double,
                                x: c_double,
                                result: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_gamma_inc_P(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gamma_inc_P_e(a: c_double,
                                x: c_double,
                                result: *mut gsl_sf_result)
                                -> c_int;
    // Incomplete Beta functions
    pub fn gsl_sf_beta_inc(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_beta_inc_e(a: c_double,
                             b: c_double,
                             x: c_double,
                             result: *mut gsl_sf_result)
                             -> c_int;

    // Gegenbauer functions
    pub fn gsl_sf_gegenpoly_1(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_2(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_3(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_1_e(lambda: c_double,
                                x: c_double,
                                result: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_gegenpoly_2_e(lambda: c_double,
                                x: c_double,
                                result: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_gegenpoly_3_e(lambda: c_double,
                                x: c_double,
                                result: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_gegenpoly_n(n: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_gegenpoly_n_e(n: c_int,
                                lambda: c_double,
                                x: c_double,
                                result: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_gegenpoly_array(nmax: c_int,
                                  lambda: c_double,
                                  x: c_double,
                                  result_array: *mut c_double)
                                  -> c_int;

    // Hypergeometric functions
    pub fn gsl_sf_hyperg_0F1(c: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_0F1_e(c: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_hyperg_1F1_int(m: c_int, n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_1F1_int_e(m: c_int,
                                   n: c_int,
                                   x: c_double,
                                   result: *mut gsl_sf_result)
                                   -> c_int;
    pub fn gsl_sf_hyperg_1F1(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_1F1_e(a: c_double,
                               b: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_hyperg_U_int(m: c_int, n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_U_int_e(m: c_int,
                                 n: c_int,
                                 x: c_double,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    pub fn gsl_sf_hyperg_U_int_e10_e(m: c_int,
                                     n: c_int,
                                     x: c_double,
                                     result: *mut gsl_sf_result_e10)
                                     -> c_int;
    pub fn gsl_sf_hyperg_U(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_U_e(a: c_double,
                             b: c_double,
                             x: c_double,
                             result: *mut gsl_sf_result)
                             -> c_int;
    pub fn gsl_sf_hyperg_U_e10_e(a: c_double,
                                 b: c_double,
                                 x: c_double,
                                 result: *mut gsl_sf_result_e10)
                                 -> c_int;
    pub fn gsl_sf_hyperg_2F1(a: c_double, b: c_double, c: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_2F1_e(a: c_double,
                               b: c_double,
                               c: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_hyperg_2F1_conj(aR: c_double,
                                  aI: c_double,
                                  c: c_double,
                                  x: c_double)
                                  -> c_double;
    pub fn gsl_sf_hyperg_2F1_conj_e(aR: c_double,
                                    aI: c_double,
                                    c: c_double,
                                    x: c_double,
                                    result: *mut gsl_sf_result)
                                    -> c_int;
    pub fn gsl_sf_hyperg_2F1_renorm(a: c_double,
                                    b: c_double,
                                    c: c_double,
                                    x: c_double)
                                    -> c_double;
    pub fn gsl_sf_hyperg_2F1_renorm_e(a: c_double,
                                      b: c_double,
                                      c: c_double,
                                      x: c_double,
                                      result: *mut gsl_sf_result)
                                      -> c_int;
    pub fn gsl_sf_hyperg_2F1_conj_renorm(aR: c_double,
                                         aI: c_double,
                                         c: c_double,
                                         x: c_double)
                                         -> c_double;
    pub fn gsl_sf_hyperg_2F1_conj_renorm_e(aR: c_double,
                                           aI: c_double,
                                           c: c_double,
                                           x: c_double,
                                           result: *mut gsl_sf_result)
                                           -> c_int;
    pub fn gsl_sf_hyperg_2F0(a: c_double, b: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_hyperg_2F0_e(a: c_double,
                               b: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;

    /// Laguerre functions
    pub fn gsl_sf_laguerre_1(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_2(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_3(a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_1_e(a: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_laguerre_2_e(a: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_laguerre_3_e(a: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_laguerre_n(n: c_int, a: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_laguerre_n_e(n: c_int,
                               a: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;

    // Lambert W functions
    pub fn gsl_sf_lambert_W0(x: c_double) -> c_double;
    pub fn gsl_sf_lambert_W0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_lambert_Wm1(x: c_double) -> c_double;
    pub fn gsl_sf_lambert_Wm1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Legendre functions
    // Legendre Polynomials
    pub fn gsl_sf_legendre_P1(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_P2(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_P3(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_P1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_legendre_P2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_legendre_P3_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_legendre_Pl(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Pl_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_legendre_Pl_array(lmax: c_int,
                                    x: c_double,
                                    result_array: *mut c_double)
                                    -> c_int;
    pub fn gsl_sf_legendre_Pl_deriv_array(lmax: c_int,
                                          x: c_double,
                                          result_array: *mut c_double,
                                          result_deriv_array: *mut c_double)
                                          -> c_int;
    pub fn gsl_sf_legendre_Q0(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Q0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_legendre_Q1(x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Q1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_legendre_Ql(l: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Ql_e(l: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Associated Legendre Polynomials and Spherical Harmonics
    pub fn gsl_sf_legendre_Plm(l: c_int, m: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_Plm_e(l: c_int,
                                 m: c_int,
                                 x: c_double,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_Plm_array(lmax: c_int,
                                     m: c_int,
                                     x: c_double,
                                     result_array: *mut c_double)
                                     -> c_int;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_Plm_deriv_array(lmax: c_int,
                                           m: c_int,
                                           x: c_double,
                                           result_array: *mut c_double,
                                           result_deriv_array: *mut c_double)
                                           -> c_int;
    pub fn gsl_sf_legendre_sphPlm(l: c_int, m: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_sphPlm_e(l: c_int,
                                    m: c_int,
                                    x: c_double,
                                    result: *mut gsl_sf_result)
                                    -> c_int;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_sphPlm_array(lmax: c_int,
                                        m: c_int,
                                        x: c_double,
                                        result_array: *mut c_double)
                                        -> c_int;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_sphPlm_deriv_array(lmax: c_int,
                                              m: c_int,
                                              x: c_double,
                                              result_array: *mut c_double,
                                              result_deriv_array: *mut c_double)
                                              -> c_int;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_legendre_array_size(lmax: c_int, m: c_int) -> c_int;
    #[cfg(feature = "v2")]
    pub fn gsl_sf_legendre_array_n(lmax: size_t) -> size_t;
    #[cfg(feature = "v2")]
    pub fn gsl_sf_legendre_array_index(l: size_t, m: size_t) -> size_t;

    pub fn gsl_sf_legendre_array(norm: c_int,
                                       lmax: size_t, x: c_double,
                                       result_array: *mut c_double) -> c_int;
    #[cfg(feature = "v2")]
    pub fn gsl_sf_legendre_deriv_array(norm: c_int,
                                       lmax: size_t, x: c_double,
                                       result_array: *mut c_double,
                                       result_deriv_array: *mut c_double) -> c_int;

    // Conical functions
    pub fn gsl_sf_conicalP_half(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_half_e(lambda: c_double,
                                  x: c_double,
                                  result: *mut gsl_sf_result)
                                  -> c_int;
    pub fn gsl_sf_conicalP_mhalf(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_mhalf_e(lambda: c_double,
                                   x: c_double,
                                   result: *mut gsl_sf_result)
                                   -> c_int;
    pub fn gsl_sf_conicalP_0(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_0_e(lambda: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_conicalP_1(lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_1_e(lambda: c_double,
                               x: c_double,
                               result: *mut gsl_sf_result)
                               -> c_int;
    pub fn gsl_sf_conicalP_sph_reg(l: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_sph_reg_e(l: c_int,
                                     lambda: c_double,
                                     x: c_double,
                                     result: *mut gsl_sf_result)
                                     -> c_int;
    pub fn gsl_sf_conicalP_cyl_reg(m: c_int, lambda: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_conicalP_cyl_reg_e(m: c_int,
                                     lambda: c_double,
                                     x: c_double,
                                     result: *mut gsl_sf_result)
                                     -> c_int;
    // Radial Functions for Hyperbolic Space
    pub fn gsl_sf_legendre_H3d_0(lambda: c_double, eta: c_double) -> c_double;
    pub fn gsl_sf_legendre_H3d_0_e(lambda: c_double,
                                   eta: c_double,
                                   result: *mut gsl_sf_result)
                                   -> c_int;
    pub fn gsl_sf_legendre_H3d_1(lambda: c_double, eta: c_double) -> c_double;
    pub fn gsl_sf_legendre_H3d_1_e(lambda: c_double,
                                   eta: c_double,
                                   result: *mut gsl_sf_result)
                                   -> c_int;
    pub fn gsl_sf_legendre_H3d(l: c_int, lambda: c_double, eta: c_double) -> c_double;
    pub fn gsl_sf_legendre_H3d_e(l: c_int,
                                 lambda: c_double,
                                 eta: c_double,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    pub fn gsl_sf_legendre_H3d_array(lmax: c_int,
                                     lambda: c_double,
                                     eta: c_double,
                                     result_array: *mut c_double)
                                     -> c_int;

    // Logarithm and Related Functions
    pub fn gsl_sf_log(x: c_double) -> c_double;
    pub fn gsl_sf_log_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_log_abs(x: c_double) -> c_double;
    pub fn gsl_sf_log_abs_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_complex_log_e(zr: c_double,
                                zi: c_double,
                                lnr: *mut gsl_sf_result,
                                theta: *mut gsl_sf_result)
                                -> c_int;
    pub fn gsl_sf_log_1plusx(x: c_double) -> c_double;
    pub fn gsl_sf_log_1plusx_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_log_1plusx_mx(x: c_double) -> c_double;
    pub fn gsl_sf_log_1plusx_mx_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Power functions
    pub fn gsl_sf_pow_int(x: c_double, n: c_int) -> c_double;
    pub fn gsl_sf_pow_int_e(x: c_double, n: c_int, result: *mut gsl_sf_result) -> c_int;

    // Psi (Digamma) functions
    // Digamma functions
    pub fn gsl_sf_psi_int(n: c_int) -> c_double;
    pub fn gsl_sf_psi_int_e(n: c_int, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_psi(x: c_double) -> c_double;
    pub fn gsl_sf_psi_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_psi_1piy(y: c_double) -> c_double;
    pub fn gsl_sf_psi_1piy_e(y: c_double, result: *mut gsl_sf_result) -> c_int;
    // Trigamma functions
    pub fn gsl_sf_psi_1_int(n: c_int) -> c_double;
    pub fn gsl_sf_psi_1_int_e(n: c_int, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_psi_1(x: c_double) -> c_double;
    pub fn gsl_sf_psi_1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Polygamma functions
    pub fn gsl_sf_psi_n(n: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_psi_n_e(n: c_int, x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Synchrotron functions
    pub fn gsl_sf_synchrotron_1(x: c_double) -> c_double;
    pub fn gsl_sf_synchrotron_1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_synchrotron_2(x: c_double) -> c_double;
    pub fn gsl_sf_synchrotron_2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Transport functions
    pub fn gsl_sf_transport_2(x: c_double) -> c_double;
    pub fn gsl_sf_transport_2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_transport_3(x: c_double) -> c_double;
    pub fn gsl_sf_transport_3_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_transport_4(x: c_double) -> c_double;
    pub fn gsl_sf_transport_4_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_transport_5(x: c_double) -> c_double;
    pub fn gsl_sf_transport_5_e(x: c_double, result: *mut gsl_sf_result) -> c_int;

    // Zeta functions
    // Riemann Zeta functions
    pub fn gsl_sf_zeta_int(n: c_int) -> c_double;
    pub fn gsl_sf_zeta_int_e(n: c_int, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_zeta(s: c_double) -> c_double;
    pub fn gsl_sf_zeta_e(s: c_double, result: *mut gsl_sf_result) -> c_int;
    // Riemann Zeta functions Minus One
    pub fn gsl_sf_zetam1_int(n: c_int) -> c_double;
    pub fn gsl_sf_zetam1_int_e(n: c_int, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_zetam1(s: c_double) -> c_double;
    pub fn gsl_sf_zetam1_e(s: c_double, result: *mut gsl_sf_result) -> c_int;
    // Hurwitz Zeta functions
    pub fn gsl_sf_hzeta(s: c_double, q: c_double) -> c_double;
    pub fn gsl_sf_hzeta_e(s: c_double, q: c_double, result: *mut gsl_sf_result) -> c_int;
    // Eta functions
    pub fn gsl_sf_eta_int(n: c_int) -> c_double;
    pub fn gsl_sf_eta_int_e(n: c_int, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_eta(s: c_double) -> c_double;
    pub fn gsl_sf_eta_e(s: c_double, result: *mut gsl_sf_result) -> c_int;

    // Elliptic Integrals
    // Legendre Form of Complete Elliptic Integrals
    pub fn gsl_sf_ellint_Kcomp(k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_Kcomp_e(k: c_double,
                                 mode: ::Mode,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    pub fn gsl_sf_ellint_Ecomp(k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_Ecomp_e(k: c_double,
                                 mode: ::Mode,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    pub fn gsl_sf_ellint_Pcomp(k: c_double, n: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_Pcomp_e(k: c_double,
                                 n: c_double,
                                 mode: ::Mode,
                                 result: *mut gsl_sf_result)
                                 -> c_int;
    // Legendre Form of Incomplete Elliptic Integrals
    pub fn gsl_sf_ellint_F(phi: c_double, k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_F_e(phi: c_double,
                             k: c_double,
                             mode: ::Mode,
                             result: *mut gsl_sf_result)
                             -> c_int;
    pub fn gsl_sf_ellint_E(phi: c_double, k: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_E_e(phi: c_double,
                             k: c_double,
                             mode: ::Mode,
                             result: *mut gsl_sf_result)
                             -> c_int;
    pub fn gsl_sf_ellint_P(phi: c_double, k: c_double, n: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_P_e(phi: c_double,
                             k: c_double,
                             n: c_double,
                             mode: ::Mode,
                             result: *mut gsl_sf_result)
                             -> c_int;
    #[cfg(feature = "v2")]
    pub fn gsl_sf_ellint_D(phi: c_double, k: c_double, mode: ::Mode) -> c_double;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_sf_ellint_D(phi: c_double, k: c_double, n: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_D_e(phi: c_double,
                             k: c_double,
                             n: c_double,
                             mode: ::Mode,
                             result: *mut gsl_sf_result)
                             -> c_int;
    // Carlson Forms
    pub fn gsl_sf_ellint_RC(x: c_double, y: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_RC_e(x: c_double,
                              y: c_double,
                              mode: ::Mode,
                              result: *mut gsl_sf_result)
                              -> c_int;
    pub fn gsl_sf_ellint_RD(x: c_double, y: c_double, z: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_RD_e(x: c_double,
                              y: c_double,
                              z: c_double,
                              mode: ::Mode,
                              result: *mut gsl_sf_result)
                              -> c_int;
    pub fn gsl_sf_ellint_RF(x: c_double, y: c_double, z: c_double, mode: ::Mode) -> c_double;
    pub fn gsl_sf_ellint_RF_e(x: c_double,
                              y: c_double,
                              z: c_double,
                              mode: ::Mode,
                              result: *mut gsl_sf_result)
                              -> c_int;
    pub fn gsl_sf_ellint_RJ(x: c_double,
                            y: c_double,
                            z: c_double,
                            p: c_double,
                            mode: ::Mode)
                            -> c_double;
    pub fn gsl_sf_ellint_RJ_e(x: c_double,
                              y: c_double,
                              z: c_double,
                              p: c_double,
                              mode: ::Mode,
                              result: *mut gsl_sf_result)
                              -> c_int;

    // Fermi-Dirac functions
    // Complete Fermi-Dirac Integrals
    pub fn gsl_sf_fermi_dirac_m1(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_m1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_fermi_dirac_0(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_0_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_fermi_dirac_1(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_1_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_fermi_dirac_2(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_2_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_fermi_dirac_int(j: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_int_e(j: c_int,
                                    x: c_double,
                                    result: *mut gsl_sf_result)
                                    -> c_int;
    pub fn gsl_sf_fermi_dirac_mhalf(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_mhalf_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_fermi_dirac_half(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_half_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_fermi_dirac_3half(x: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_3half_e(x: c_double, result: *mut gsl_sf_result) -> c_int;
    // Incomplete Fermi-Dirac Integrals
    pub fn gsl_sf_fermi_dirac_inc_0(x: c_double, b: c_double) -> c_double;
    pub fn gsl_sf_fermi_dirac_inc_0_e(x: c_double,
                                      b: c_double,
                                      result: *mut gsl_sf_result)
                                      -> c_int;

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

    // Mathieu functions
    // Mathieu functions Workspace
    pub fn gsl_sf_mathieu_alloc(n: size_t, qmax: c_double) -> *mut gsl_sf_mathieu_workspace;
    pub fn gsl_sf_mathieu_free(work: *mut gsl_sf_mathieu_workspace);
    // Mathieu functions Characteristic Values
    pub fn gsl_sf_mathieu_a_e(n: c_int, q: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_mathieu_b_e(n: c_int, q: c_double, result: *mut gsl_sf_result) -> c_int;
    pub fn gsl_sf_mathieu_a_array(order_min: c_int,
                                  order_max: c_int,
                                  q: c_double,
                                  work: *mut gsl_sf_mathieu_workspace,
                                  result_array: *mut c_double)
                                  -> c_int;
    pub fn gsl_sf_mathieu_b_array(order_min: c_int,
                                  order_max: c_int,
                                  q: c_double,
                                  work: *mut gsl_sf_mathieu_workspace,
                                  result_array: *mut c_double)
                                  -> c_int;
    // Angular Mathieu functions
    pub fn gsl_sf_mathieu_ce_e(n: c_int,
                             q: c_double,
                             x: c_double,
                             result: *mut gsl_sf_result)
                             -> c_int;
    pub fn gsl_sf_mathieu_se_e(n: c_int,
                             q: c_double,
                             x: c_double,
                             result: *mut gsl_sf_result)
                             -> c_int;
    pub fn gsl_sf_mathieu_ce_array(nmin: c_int,
                                   nmax: c_int,
                                   q: c_double,
                                   x: c_double,
                                   work: *mut gsl_sf_mathieu_workspace,
                                   result_array: *mut c_double)
                                   -> c_int;
    pub fn gsl_sf_mathieu_se_array(nmin: c_int,
                                   nmax: c_int,
                                   q: c_double,
                                   x: c_double,
                                   work: *mut gsl_sf_mathieu_workspace,
                                   result_array: *mut c_double)
                                   -> c_int;
    // Radial Mathieu functions
    pub fn gsl_sf_mathieu_Mc_e(j: c_int,
                             n: c_int,
                             q: c_double,
                             x: c_double,
                             result: *mut gsl_sf_result)
                             -> c_int;
    pub fn gsl_sf_mathieu_Ms_e(j: c_int,
                             n: c_int,
                             q: c_double,
                             x: c_double,
                             result: *mut gsl_sf_result)
                             -> c_int;
    pub fn gsl_sf_mathieu_Mc_array(j: c_int,
                                   nmin: c_int,
                                   nmax: c_int,
                                   q: c_double,
                                   x: c_double,
                                   work: *mut gsl_sf_mathieu_workspace,
                                   result_array: *mut c_double)
                                   -> c_int;
    pub fn gsl_sf_mathieu_Ms_array(j: c_int,
                                   nmin: c_int,
                                   nmax: c_int,
                                   q: c_double,
                                   x: c_double,
                                   work: *mut gsl_sf_mathieu_workspace,
                                   result_array: *mut c_double)
                                   -> c_int;

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
    /*
    pub fn gsl_complex_float_arg(z: gsl_complex_float) -> c_float;
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
    pub fn gsl_complex_float_arccoth(z: gsl_complex_float) -> gsl_complex_float;
    */



    // Basis Splines
    pub fn gsl_bspline_alloc(k: size_t, nbreak: size_t) -> *mut gsl_bspline_workspace;
    pub fn gsl_bspline_free(w: *mut gsl_bspline_workspace);
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_alloc(k: size_t) -> *mut gsl_bspline_deriv_workspace;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_free(w: *mut gsl_bspline_deriv_workspace);
    pub fn gsl_bspline_knots(breakpts: *const gsl_vector,
                             w: *mut gsl_bspline_workspace)
                             -> c_int;
    pub fn gsl_bspline_knots_uniform(a: c_double,
                                     b: c_double,
                                     w: *mut gsl_bspline_workspace)
                                     -> c_int;
    pub fn gsl_bspline_eval(x: c_double,
                            B: *mut gsl_vector,
                            w: *mut gsl_bspline_workspace)
                            -> c_int;
    pub fn gsl_bspline_eval_nonzero(x: c_double,
                                    Bk: *mut gsl_vector,
                                    istart: *mut size_t,
                                    iend: *mut size_t,
                                    w: *mut gsl_bspline_workspace)
                                    -> c_int;
    pub fn gsl_bspline_ncoeffs(w: *mut gsl_bspline_workspace) -> size_t;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_eval(x: c_double,
                                  nderiv: size_t,
                                  dB: *mut gsl_matrix,
                                  w: *mut gsl_bspline_workspace,
                                  dw: *mut gsl_bspline_deriv_workspace)
                                  -> c_int;
    #[cfg(not(feature = "v2"))]
    pub fn gsl_bspline_deriv_eval_nonzero(x: c_double,
                                          nderiv: size_t,
                                          Bk: *mut gsl_matrix,
                                          istart: *mut size_t,
                                          iend: *mut size_t,
                                          w: *mut gsl_bspline_workspace,
                                          dw: *mut gsl_bspline_deriv_workspace)
                                          -> c_int;
    pub fn gsl_bspline_greville_abscissa(i: size_t, w: *mut gsl_bspline_workspace) -> c_double;

    // Fit functions
    pub fn gsl_fit_linear(x: *const c_double,
                          xstride: size_t,
                          y: *const c_double,
                          ystride: size_t,
                          n: size_t,
                          c0: *mut c_double,
                          c1: *mut c_double,
                          cov00: *mut c_double,
                          cov01: *mut c_double,
                          cov11: *mut c_double,
                          sumsq: *mut c_double)
                          -> c_int;
    pub fn gsl_fit_wlinear(x: *const c_double,
                           xstride: size_t,
                           w: *const c_double,
                           wstride: size_t,
                           y: *const c_double,
                           ystride: size_t,
                           n: size_t,
                           c0: *mut c_double,
                           c1: *mut c_double,
                           cov00: *mut c_double,
                           cov01: *mut c_double,
                           cov11: *mut c_double,
                           chisq: *mut c_double)
                           -> c_int;
    pub fn gsl_fit_linear_est(x: c_double,
                              c0: c_double,
                              c1: c_double,
                              cov00: c_double,
                              cov01: c_double,
                              cov11: c_double,
                              y: *mut c_double,
                              y_err: *mut c_double)
                              -> c_int;
    pub fn gsl_fit_mul(x: *const c_double,
                       xstride: size_t,
                       y: *const c_double,
                       ystride: size_t,
                       n: size_t,
                       c1: *mut c_double,
                       cov11: *mut c_double,
                       sumsq: *mut c_double)
                       -> c_int;
    pub fn gsl_fit_wmul(x: *const c_double,
                        xstride: size_t,
                        w: *const c_double,
                        wstride: size_t,
                        y: *const c_double,
                        ystride: size_t,
                        n: size_t,
                        c1: *mut c_double,
                        cov11: *mut c_double,
                        sumsq: *mut c_double)
                        -> c_int;
    pub fn gsl_fit_mul_est(x: c_double,
                           c1: c_double,
                           cov11: c_double,
                           y: *mut c_double,
                           y_err: *mut c_double)
                           -> c_int;

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
    pub fn gsl_rng_set(r: *mut gsl_rng, s: c_ulong);
    pub fn gsl_rng_free(r: *mut gsl_rng);
    pub fn gsl_rng_get(r: *mut gsl_rng) -> c_ulong;
    pub fn gsl_rng_uniform(r: *mut gsl_rng) -> c_double;
    pub fn gsl_rng_uniform_pos(r: *mut gsl_rng) -> c_double;
    pub fn gsl_rng_uniform_int(r: *mut gsl_rng, n: c_ulong) -> c_ulong;
    pub fn gsl_rng_name(r: *const gsl_rng) -> *const c_char;
    pub fn gsl_rng_max(r: *const gsl_rng) -> c_ulong;
    pub fn gsl_rng_min(r: *const gsl_rng) -> c_ulong;
    pub fn gsl_rng_state(r: *const gsl_rng) -> *mut c_void;
    pub fn gsl_rng_size(r: *const gsl_rng) -> size_t;
    pub fn gsl_rng_types_setup() -> *const *mut gsl_rng_type;
    pub fn gsl_rng_memcpy(dest: *mut gsl_rng, src: *const gsl_rng) -> c_int;
    pub fn gsl_rng_clone(r: *const gsl_rng) -> *mut gsl_rng;
    pub fn gsl_rng_env_setup() -> *const gsl_rng_type;

    // Permutation struct
    pub fn gsl_permutation_alloc(size: size_t) -> *mut gsl_permutation;
    pub fn gsl_permutation_calloc(size: size_t) -> *mut gsl_permutation;
    pub fn gsl_permutation_init(p: *mut gsl_permutation);
    pub fn gsl_permutation_free(p: *mut gsl_permutation);
    pub fn gsl_permutation_memcpy(dest: *mut gsl_permutation,
                                  src: *const gsl_permutation)
                                  -> c_int;
    pub fn gsl_permutation_get(p: *const gsl_permutation, i: size_t) -> size_t;
    pub fn gsl_permutation_swap(p: *mut gsl_permutation, i: size_t, j: size_t) -> c_int;
    pub fn gsl_permutation_size(p: *const gsl_permutation) -> size_t;
    //pub fn gsl_permutation_data(p: *const gsl_permutation) -> *mut size_t;
    pub fn gsl_permutation_valid(p: *const gsl_permutation) -> c_int;
    pub fn gsl_permutation_reverse(p: *mut gsl_permutation);
    pub fn gsl_permutation_inverse(inv: *mut gsl_permutation,
                                   p: *const gsl_permutation)
                                   -> c_int;
    pub fn gsl_permutation_next(p: *mut gsl_permutation) -> c_int;
    pub fn gsl_permutation_prev(p: *mut gsl_permutation) -> c_int;
    pub fn gsl_permute(p: *const size_t,
                       data: *mut c_double,
                       stride: size_t,
                       n: size_t)
                       -> c_int;
    pub fn gsl_permute_inverse(p: *const size_t,
                               data: *mut c_double,
                               stride: size_t,
                               n: size_t)
                               -> c_int;
    pub fn gsl_permute_vector(p: *const gsl_permutation, v: *mut gsl_vector) -> c_int;
    pub fn gsl_permute_vector_inverse(p: *const gsl_permutation,
                                      v: *mut gsl_vector)
                                      -> c_int;
    pub fn gsl_permutation_mul(p: *mut gsl_permutation,
                               pa: *const gsl_permutation,
                               pb: *const gsl_permutation)
                               -> c_int;
    pub fn gsl_permutation_linear_to_canonical(q: *mut gsl_permutation,
                                               p: *const gsl_permutation)
                                               -> c_int;
    pub fn gsl_permutation_canonical_to_linear(p: *mut gsl_permutation,
                                               q: *const gsl_permutation)
                                               -> c_int;
    pub fn gsl_permutation_inversions(p: *const gsl_permutation) -> size_t;
    pub fn gsl_permutation_linear_cycles(p: *const gsl_permutation) -> size_t;
    pub fn gsl_permutation_canonical_cycles(p: *const gsl_permutation) -> size_t;

    // Sorting functions
    // Sorting objects
    //pub fn gsl_heapsort(array: *mut c_void, count: size_t, size: size_t, compare: compare_fn);
    //pub fn gsl_heapsort_index(p: *mut size_t, array: *const c_void, count: size_t, size: size_t, compare: compare_fn) -> c_int;
    // Sorting vectors
    pub fn gsl_sort(data: *mut c_double, stride: size_t, n: size_t);
    pub fn gsl_sort2(data1: *mut c_double,
                     stride1: size_t,
                     data2: *mut c_double,
                     stride2: size_t,
                     n: size_t);
    pub fn gsl_sort_vector(v: *mut gsl_vector);
    pub fn gsl_sort_vector2(v1: *mut gsl_vector, v2: *mut gsl_vector);
    pub fn gsl_sort_index(p: *mut size_t, data: *const c_double, stride: size_t, n: size_t);
    pub fn gsl_sort_vector_index(p: *mut gsl_permutation, v: *const gsl_vector) -> c_int;
    // Selecting the k smallest or largest elements
    pub fn gsl_sort_smallest(dest: *mut c_double,
                             k: size_t,
                             src: *const c_double,
                             stride: size_t,
                             n: size_t)
                             -> c_int;
    pub fn gsl_sort_largest(dest: *mut c_double,
                            k: size_t,
                            src: *const c_double,
                            stride: size_t,
                            n: size_t)
                            -> c_int;
    pub fn gsl_sort_vector_smallest(dest: *mut c_double,
                                    k: size_t,
                                    v: *const gsl_vector)
                                    -> c_int;
    pub fn gsl_sort_vector_largest(dest: *mut c_double,
                                   k: size_t,
                                   v: *const gsl_vector)
                                   -> c_int;
    pub fn gsl_sort_smallest_index(p: *mut size_t,
                                   k: size_t,
                                   src: *const c_double,
                                   stride: size_t,
                                   n: size_t)
                                   -> c_int;
    pub fn gsl_sort_largest_index(p: *mut size_t,
                                  k: size_t,
                                  src: *const c_double,
                                  stride: size_t,
                                  n: size_t)
                                  -> c_int;
    pub fn gsl_sort_vector_smallest_index(p: *mut size_t,
                                          k: size_t,
                                          v: *const gsl_vector)
                                          -> c_int;
    pub fn gsl_sort_vector_largest_index(p: *mut size_t,
                                         k: size_t,
                                         v: *const gsl_vector)
                                         -> c_int;

    // Chebyshev Approximations
    // Creation and Calculation of Chebyshev Series
    pub fn gsl_cheb_alloc(n: size_t) -> *mut gsl_cheb_series;
    pub fn gsl_cheb_free(cs: *mut gsl_cheb_series);
    // Auxiliary functions
    pub fn gsl_cheb_order(cs: *const gsl_cheb_series) -> size_t;
    pub fn gsl_cheb_size(cs: *const gsl_cheb_series) -> size_t;
    // Chebyshev Series Evaluation
    pub fn gsl_cheb_eval(cs: *const gsl_cheb_series, x: c_double) -> c_double;
    pub fn gsl_cheb_eval_err(cs: *const gsl_cheb_series,
                             x: c_double,
                             result: *mut c_double,
                             abs_err: *mut c_double)
                             -> c_int;
    pub fn gsl_cheb_eval_n(cs: *const gsl_cheb_series, order: size_t, x: c_double) -> c_double;
    pub fn gsl_cheb_eval_n_err(cs: *const gsl_cheb_series,
                               order: size_t,
                               x: c_double,
                               result: *mut c_double,
                               abs_err: *mut c_double)
                               -> c_int;
    // Derivatives and Integrals
    pub fn gsl_cheb_calc_deriv(cs: *mut gsl_cheb_series,
                               deriv: *const gsl_cheb_series)
                               -> c_int;
    pub fn gsl_cheb_calc_integ(cs: *mut gsl_cheb_series,
                               integ: *const gsl_cheb_series)
                               -> c_int;

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
    pub fn gsl_combination_memcpy(dest: *mut gsl_combination,
                                  src: *const gsl_combination)
                                  -> c_int;
    // Accessing combination elements
    pub fn gsl_combination_get(c: *const gsl_combination, i: size_t) -> size_t;
    // Combination properties
    pub fn gsl_combination_n(c: *const gsl_combination) -> size_t;
    pub fn gsl_combination_k(c: *const gsl_combination) -> size_t;
    pub fn gsl_combination_valid(c: *mut gsl_combination) -> c_int;
    // Combination functions
    pub fn gsl_combination_next(c: *mut gsl_combination) -> c_int;
    pub fn gsl_combination_prev(c: *mut gsl_combination) -> c_int;

    // Polynomials
    // Polynomial Evaluation
    pub fn gsl_poly_eval(c: *const c_double, len: c_int, x: c_double) -> c_double;
    pub fn gsl_poly_complex_eval(c: *const c_double, len: c_int, z: gsl_complex) -> gsl_complex;
    pub fn gsl_complex_poly_complex_eval(c: *const gsl_complex,
                                         len: c_int,
                                         z: gsl_complex)
                                         -> gsl_complex;
    pub fn gsl_poly_eval_derivs(c: *const c_double,
                                lenc: size_t,
                                x: c_double,
                                res: *mut c_double,
                                lenres: size_t)
                                -> c_int;
    // Divided Difference Representation of Polynomials
    pub fn gsl_poly_dd_init(dd: *mut c_double,
                            xa: *const c_double,
                            ya: *const c_double,
                            size: size_t)
                            -> c_int;
    pub fn gsl_poly_dd_eval(dd: *const c_double,
                            xa: *const c_double,
                            size: size_t,
                            x: c_double)
                            -> c_double;
    pub fn gsl_poly_dd_taylor(c: *mut c_double,
                              xp: c_double,
                              dd: *const c_double,
                              xa: *const c_double,
                              size: size_t,
                              w: *mut c_double)
                              -> c_int;
    pub fn gsl_poly_dd_hermite_init(dd: *mut c_double,
                                    za: *mut c_double,
                                    xa: *const c_double,
                                    ya: *const c_double,
                                    dya: *const c_double,
                                    size: size_t)
                                    -> c_int;
    // Quadratic Equations
    pub fn gsl_poly_solve_quadratic(a: c_double,
                                    b: c_double,
                                    c: c_double,
                                    x0: *mut c_double,
                                    x1: *mut c_double)
                                    -> c_int;
    pub fn gsl_poly_complex_solve_quadratic(a: c_double,
                                            b: c_double,
                                            c: c_double,
                                            x0: *mut gsl_complex,
                                            x1: *mut gsl_complex)
                                            -> c_int;
    // Cubic Equations
    pub fn gsl_poly_solve_cubic(a: c_double,
                                b: c_double,
                                c: c_double,
                                x0: *mut c_double,
                                x1: *mut c_double,
                                x2: *mut c_double)
                                -> c_int;
    pub fn gsl_poly_complex_solve_cubic(a: c_double,
                                        b: c_double,
                                        c: c_double,
                                        x0: *mut gsl_complex,
                                        x1: *mut gsl_complex,
                                        x2: *mut gsl_complex)
                                        -> c_int;
    // General Polynomial Equations
    pub fn gsl_poly_complex_workspace_alloc(n: size_t) -> *mut gsl_poly_complex_workspace;
    pub fn gsl_poly_complex_workspace_free(w: *mut gsl_poly_complex_workspace);
    pub fn gsl_poly_complex_solve(a: *const c_double,
                                  n: size_t,
                                  w: *mut gsl_poly_complex_workspace,
                                  z: gsl_complex_packed_ptr)
                                  -> c_int;

    // Discrete Hankel functions
    pub fn gsl_dht_alloc(size: size_t) -> *mut gsl_dht;
    pub fn gsl_dht_init(t: *mut gsl_dht, nu: c_double, xmax: c_double) -> c_int;
    pub fn gsl_dht_new(size: size_t, nu: c_double, xmax: c_double) -> *mut gsl_dht;
    pub fn gsl_dht_free(t: *mut gsl_dht);
    pub fn gsl_dht_apply(t: *const gsl_dht,
                         f_in: *const c_double,
                         f_out: *mut c_double)
                         -> c_int;
    pub fn gsl_dht_x_sample(t: *const gsl_dht, n: c_int) -> c_double;
    pub fn gsl_dht_k_sample(t: *const gsl_dht, n: c_int) -> c_double;

    // Fast Fourier Transforms
    // Radix-2 FFT routines for complex data
    pub fn gsl_fft_complex_radix2_forward(data: gsl_complex_packed_array,
                                          stride: size_t,
                                          n: size_t)
                                          -> c_int;
    pub fn gsl_fft_complex_radix2_transform(data: gsl_complex_packed_array,
                                            stride: size_t,
                                            n: size_t,
                                            sign: c_int)
                                            -> c_int;
    pub fn gsl_fft_complex_radix2_backward(data: gsl_complex_packed_array,
                                           stride: size_t,
                                           n: size_t)
                                           -> c_int;
    pub fn gsl_fft_complex_radix2_inverse(data: gsl_complex_packed_array,
                                          stride: size_t,
                                          n: size_t)
                                          -> c_int;
    pub fn gsl_fft_complex_radix2_dif_forward(data: gsl_complex_packed_array,
                                              stride: size_t,
                                              n: size_t)
                                              -> c_int;
    pub fn gsl_fft_complex_radix2_dif_transform(data: gsl_complex_packed_array,
                                                stride: size_t,
                                                n: size_t,
                                                sign: c_int)
                                                -> c_int;
    pub fn gsl_fft_complex_radix2_dif_backward(data: gsl_complex_packed_array,
                                               stride: size_t,
                                               n: size_t)
                                               -> c_int;
    pub fn gsl_fft_complex_radix2_dif_inverse(data: gsl_complex_packed_array,
                                              stride: size_t,
                                              n: size_t)
                                              -> c_int;
    // Mixed-radix FFT routines for complex data
    pub fn gsl_fft_complex_wavetable_alloc(n: size_t) -> *mut gsl_fft_complex_wavetable;
    pub fn gsl_fft_complex_wavetable_free(w: *mut gsl_fft_complex_wavetable);
    pub fn gsl_fft_complex_workspace_alloc(n: size_t) -> *mut gsl_fft_complex_workspace;
    pub fn gsl_fft_complex_workspace_free(w: *mut gsl_fft_complex_workspace);
    pub fn gsl_fft_complex_forward(data: gsl_complex_packed_array,
                                   stride: size_t,
                                   n: size_t,
                                   wavetable: *const gsl_fft_complex_wavetable,
                                   work: *mut gsl_fft_complex_workspace)
                                   -> c_int;
    pub fn gsl_fft_complex_transform(data: gsl_complex_packed_array,
                                     stride: size_t,
                                     n: size_t,
                                     wavetable: *const gsl_fft_complex_wavetable,
                                     work: *mut gsl_fft_complex_workspace,
                                     sign: c_int)
                                     -> c_int;
    pub fn gsl_fft_complex_backward(data: gsl_complex_packed_array,
                                    stride: size_t,
                                    n: size_t,
                                    wavetable: *const gsl_fft_complex_wavetable,
                                    work: *mut gsl_fft_complex_workspace)
                                    -> c_int;
    pub fn gsl_fft_complex_inverse(data: gsl_complex_packed_array,
                                   stride: size_t,
                                   n: size_t,
                                   wavetable: *const gsl_fft_complex_wavetable,
                                   work: *mut gsl_fft_complex_workspace)
                                   -> c_int;
    // Radix-2 FFT routines for real data
    pub fn gsl_fft_real_radix2_transform(data: *mut c_double,
                                         stride: size_t,
                                         n: size_t)
                                         -> c_int;
    pub fn gsl_fft_halfcomplex_radix2_inverse(data: *mut c_double,
                                              stride: size_t,
                                              n: size_t)
                                              -> c_int;
    pub fn gsl_fft_halfcomplex_radix2_backward(data: *mut c_double,
                                               stride: size_t,
                                               n: size_t)
                                               -> c_int;
    pub fn gsl_fft_halfcomplex_radix2_unpack(halfcomplex_coefficient: *mut c_double,
                                             complex_coefficient: gsl_complex_packed_array,
                                             stride: size_t,
                                             n: size_t)
                                             -> c_int;

    // Histograms
    // Histogram allocation
    pub fn gsl_histogram_alloc(n: size_t) -> *mut gsl_histogram;
    pub fn gsl_histogram_set_ranges(h: *mut gsl_histogram,
                                    range: *const c_double,
                                    size: size_t)
                                    -> c_int;
    pub fn gsl_histogram_set_ranges_uniform(h: *mut gsl_histogram,
                                            xmin: c_double,
                                            xmax: c_double)
                                            -> c_int;
    pub fn gsl_histogram_free(h: *mut gsl_histogram);
    // Copying Histograms
    pub fn gsl_histogram_memcpy(dest: *mut gsl_histogram,
                                src: *const gsl_histogram)
                                -> c_int;
    pub fn gsl_histogram_clone(src: *const gsl_histogram) -> *mut gsl_histogram;
    // Updating and accessing histogram elements
    pub fn gsl_histogram_increment(h: *mut gsl_histogram, x: c_double) -> c_int;
    pub fn gsl_histogram_accumulate(h: *mut gsl_histogram,
                                    x: c_double,
                                    weigth: c_double)
                                    -> c_int;
    pub fn gsl_histogram_get(h: *const gsl_histogram, i: size_t) -> c_double;
    pub fn gsl_histogram_get_range(h: *const gsl_histogram,
                                   i: size_t,
                                   lower: *mut c_double,
                                   upper: *mut c_double)
                                   -> c_int;
    pub fn gsl_histogram_max(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_min(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_bins(h: *const gsl_histogram) -> size_t;
    pub fn gsl_histogram_reset(h: *mut gsl_histogram);
    pub fn gsl_histogram_find(h: *const gsl_histogram,
                              x: c_double,
                              i: *mut size_t)
                              -> c_int;
    pub fn gsl_histogram_max_val(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_max_bin(h: *const gsl_histogram) -> size_t;
    pub fn gsl_histogram_min_val(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_min_bin(h: *const gsl_histogram) -> size_t;
    pub fn gsl_histogram_mean(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_sigma(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_sum(h: *const gsl_histogram) -> c_double;
    pub fn gsl_histogram_equal_bins_p(h1: *const gsl_histogram, h2: *const gsl_histogram) -> c_int;
    pub fn gsl_histogram_add(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> c_int;
    pub fn gsl_histogram_sub(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> c_int;
    pub fn gsl_histogram_mul(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> c_int;
    pub fn gsl_histogram_div(h1: *mut gsl_histogram, h2: *const gsl_histogram) -> c_int;
    pub fn gsl_histogram_scale(h1: *mut gsl_histogram, scale: c_double) -> c_int;
    pub fn gsl_histogram_shift(h1: *mut gsl_histogram, offset: c_double) -> c_int;
    // The histogram probability distribution struct
    pub fn gsl_histogram_pdf_alloc(n: size_t) -> *mut gsl_histogram_pdf;
    pub fn gsl_histogram_pdf_init(p: *mut gsl_histogram_pdf,
                                  h: *const gsl_histogram)
                                  -> c_int;
    pub fn gsl_histogram_pdf_free(p: *mut gsl_histogram_pdf);
    pub fn gsl_histogram_pdf_sample(p: *const gsl_histogram_pdf, r: c_double) -> c_double;
    // 2D Histogram allocation
    pub fn gsl_histogram2d_alloc(nx: size_t, ny: size_t) -> *mut gsl_histogram2d;
    pub fn gsl_histogram2d_set_ranges(h: *mut gsl_histogram2d,
                                      xrange: *const c_double,
                                      xsize: size_t,
                                      yrange: *const c_double,
                                      ysize: size_t)
                                      -> c_int;
    pub fn gsl_histogram2d_set_ranges_uniform(h: *mut gsl_histogram2d,
                                              xmin: c_double,
                                              xmax: c_double,
                                              ymin: c_double,
                                              ymax: c_double)
                                              -> c_int;
    pub fn gsl_histogram2d_free(h: *mut gsl_histogram2d);
    pub fn gsl_histogram2d_memcpy(dest: *mut gsl_histogram2d,
                                  src: *const gsl_histogram2d)
                                  -> c_int;
    pub fn gsl_histogram2d_clone(src: *const gsl_histogram2d) -> *mut gsl_histogram2d;
    pub fn gsl_histogram2d_increment(h: *mut gsl_histogram2d,
                                     x: c_double,
                                     y: c_double)
                                     -> c_int;
    pub fn gsl_histogram2d_accumulate(h: *mut gsl_histogram2d,
                                      x: c_double,
                                      y: c_double,
                                      weight: c_double)
                                      -> c_int;
    pub fn gsl_histogram2d_get(h: *const gsl_histogram2d, i: size_t, j: size_t) -> c_double;
    pub fn gsl_histogram2d_get_xrange(h: *const gsl_histogram2d,
                                      i: size_t,
                                      xlower: *mut c_double,
                                      xupper: *mut c_double)
                                      -> c_int;
    pub fn gsl_histogram2d_get_yrange(h: *const gsl_histogram2d,
                                      j: size_t,
                                      ylower: *mut c_double,
                                      yupper: *mut c_double)
                                      -> c_int;
    pub fn gsl_histogram2d_xmax(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_xmin(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_nx(h: *const gsl_histogram2d) -> size_t;
    pub fn gsl_histogram2d_ymax(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_ymin(h: *const gsl_histogram2d) -> c_double;
    pub fn gsl_histogram2d_ny(h: *const gsl_histogram2d) -> size_t;
    pub fn gsl_histogram2d_reset(h: *mut gsl_histogram2d);
    pub fn gsl_histogram2d_find(h: *const gsl_histogram2d,
                                x: c_double,
                                y: c_double,
                                i: *mut size_t,
                                j: *mut size_t)
                                -> c_int;
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
    pub fn gsl_histogram2d_equal_bins_p(h1: *const gsl_histogram2d,
                                        h2: *const gsl_histogram2d)
                                        -> c_int;
    pub fn gsl_histogram2d_add(h1: *mut gsl_histogram2d,
                               h2: *const gsl_histogram2d)
                               -> c_int;
    pub fn gsl_histogram2d_sub(h1: *mut gsl_histogram2d,
                               h2: *const gsl_histogram2d)
                               -> c_int;
    pub fn gsl_histogram2d_mul(h1: *mut gsl_histogram2d,
                               h2: *const gsl_histogram2d)
                               -> c_int;
    pub fn gsl_histogram2d_div(h1: *mut gsl_histogram2d,
                               h2: *const gsl_histogram2d)
                               -> c_int;
    pub fn gsl_histogram2d_scale(h1: *mut gsl_histogram2d, scale: c_double) -> c_int;
    pub fn gsl_histogram2d_shift(h1: *mut gsl_histogram2d, offset: c_double) -> c_int;
    // Resampling from 2D histograms
    pub fn gsl_histogram2d_pdf_alloc(nx: size_t, ny: size_t) -> *mut gsl_histogram2d_pdf;
    pub fn gsl_histogram2d_pdf_init(p: *mut gsl_histogram2d_pdf,
                                    h: *const gsl_histogram2d)
                                    -> c_int;
    pub fn gsl_histogram2d_pdf_free(p: *mut gsl_histogram2d_pdf);
    pub fn gsl_histogram2d_pdf_sample(p: *const gsl_histogram2d_pdf,
                                      r1: c_double,
                                      r2: c_double,
                                      x: *mut c_double,
                                      y: *mut c_double)
                                      -> c_int;

    // QAG adaptive integration
    pub fn gsl_integration_workspace_alloc(n: size_t) -> *mut gsl_integration_workspace;
    pub fn gsl_integration_workspace_free(w: *mut gsl_integration_workspace);
    // QAWS adaptive integration for singular functions
    pub fn gsl_integration_qaws_table_alloc(alpha: c_double,
                                            beta: c_double,
                                            mu: c_int,
                                            nu: c_int)
                                            -> *mut gsl_integration_qaws_table;
    pub fn gsl_integration_qaws_table_set(t: *mut gsl_integration_qaws_table,
                                          alpha: c_double,
                                          beta: c_double,
                                          mu: c_int,
                                          nu: c_int)
                                          -> c_int;
    pub fn gsl_integration_qaws_table_free(t: *mut gsl_integration_qaws_table);
    // QAWO adaptive integration for oscillatory functions
    pub fn gsl_integration_qawo_table_alloc(omega: c_double,
                                            l: c_double,
                                            sine: c_int,
                                            n: size_t)
                                            -> *mut gsl_integration_qawo_table;
    pub fn gsl_integration_qawo_table_set(t: *mut gsl_integration_qawo_table,
                                          omega: c_double,
                                          l: c_double,
                                          sine: c_int)
                                          -> c_int;
    pub fn gsl_integration_qawo_table_set_length(t: *mut gsl_integration_qawo_table,
                                                 l: c_double)
                                                 -> c_int;
    pub fn gsl_integration_qawo_table_free(t: *mut gsl_integration_qawo_table);
    // CQUAD doubly-adaptive integration
    pub fn gsl_integration_cquad_workspace_alloc(n: size_t)
                                                 -> *mut gsl_integration_cquad_workspace;
    pub fn gsl_integration_cquad_workspace_free(w: *mut gsl_integration_cquad_workspace);
    // Gauss-Legendre integration
    pub fn gsl_integration_glfixed_table_alloc(n: size_t) -> *mut gsl_integration_glfixed_table;
    pub fn gsl_integration_glfixed_point(a: c_double,
                                         b: c_double,
                                         i: size_t,
                                         xi: *mut c_double,
                                         wi: *mut c_double,
                                         t: *const gsl_integration_glfixed_table)
                                         -> c_int;
    pub fn gsl_integration_glfixed_table_free(t: *mut gsl_integration_glfixed_table);

    // Interpolation Functions
    pub fn gsl_interp_alloc(t: *const gsl_interp_type, size: size_t) -> *mut gsl_interp;
    pub fn gsl_interp_init(interp: *mut gsl_interp,
                           xa: *const c_double,
                           ya: *const c_double,
                           size: size_t)
                           -> c_int;
    pub fn gsl_interp_free(interp: *mut gsl_interp);
    pub fn gsl_interp_min_size(interp: *const gsl_interp) -> c_uint;
    pub fn gsl_interp_name(interp: *const gsl_interp) -> *const c_char;
    // Interpolation Types
    pub fn gsl_interp_type_min_size(t: *const gsl_interp_type) -> c_uint;
    // Index Look-up and Acceleration
    pub fn gsl_interp_accel_find(a: *mut ::InterpAccel,
                                 x_array: *const c_double,
                                 size: size_t,
                                 x: c_double)
                                 -> size_t;
    pub fn gsl_interp_bsearch(x_array: *const c_double,
                              x: c_double,
                              index_lo: size_t,
                              index_hi: size_t)
                              -> size_t;
    // Evaluation of Interpolating Functions
    pub fn gsl_interp_eval(interp: *const gsl_interp,
                           xa: *const c_double,
                           ya: *const c_double,
                           x: c_double,
                           acc: *mut ::InterpAccel)
                           -> c_double;
    pub fn gsl_interp_eval_e(interp: *const gsl_interp,
                             xa: *const c_double,
                             ya: *const c_double,
                             x: c_double,
                             acc: *mut ::InterpAccel,
                             y: *mut c_double)
                             -> c_int;
    pub fn gsl_interp_eval_deriv(interp: *const gsl_interp,
                                 xa: *const c_double,
                                 ya: *const c_double,
                                 x: c_double,
                                 acc: *mut ::InterpAccel)
                                 -> c_double;
    pub fn gsl_interp_eval_deriv_e(interp: *const gsl_interp,
                                   xa: *const c_double,
                                   ya: *const c_double,
                                   x: c_double,
                                   acc: *mut ::InterpAccel,
                                   d: *mut c_double)
                                   -> c_int;
    pub fn gsl_interp_eval_deriv2(interp: *const gsl_interp,
                                  xa: *const c_double,
                                  ya: *const c_double,
                                  x: c_double,
                                  acc: *mut ::InterpAccel)
                                  -> c_double;
    pub fn gsl_interp_eval_deriv2_e(interp: *const gsl_interp,
                                    xa: *const c_double,
                                    ya: *const c_double,
                                    x: c_double,
                                    acc: *mut ::InterpAccel,
                                    d2: *mut c_double)
                                    -> c_int;
    pub fn gsl_interp_eval_integ(interp: *const gsl_interp,
                                 xa: *const c_double,
                                 ya: *const c_double,
                                 a: c_double,
                                 b: c_double,
                                 acc: *mut ::InterpAccel)
                                 -> c_double;
    pub fn gsl_interp_eval_integ_e(interp: *const gsl_interp,
                                   xa: *const c_double,
                                   ya: *const c_double,
                                   a: c_double,
                                   b: c_double,
                                   acc: *mut ::InterpAccel,
                                   result: *mut c_double)
                                   -> c_int;
    // Higher-level Interface
    pub fn gsl_spline_alloc(t: *const gsl_interp_type, size: size_t) -> *mut gsl_spline;
    pub fn gsl_spline_init(spline: *mut gsl_spline,
                           xa: *const c_double,
                           ya: *const c_double,
                           size: size_t)
                           -> c_int;
    pub fn gsl_spline_free(spline: *mut gsl_spline);
    pub fn gsl_spline_min_size(spline: *const gsl_spline) -> c_uint;
    pub fn gsl_spline_name(spline: *const gsl_spline) -> *const c_char;
    pub fn gsl_spline_eval(spline: *const gsl_spline,
                           x: c_double,
                           acc: *mut ::InterpAccel)
                           -> c_double;
    pub fn gsl_spline_eval_e(spline: *const gsl_spline,
                             x: c_double,
                             acc: *mut ::InterpAccel,
                             y: *mut c_double)
                             -> c_int;
    pub fn gsl_spline_eval_deriv(spline: *const gsl_spline,
                                 x: c_double,
                                 acc: *mut ::InterpAccel)
                                 -> c_double;
    pub fn gsl_spline_eval_deriv_e(spline: *const gsl_spline,
                                   x: c_double,
                                   acc: *mut ::InterpAccel,
                                   d: *mut c_double)
                                   -> c_int;
    pub fn gsl_spline_eval_deriv2(spline: *const gsl_spline,
                                  x: c_double,
                                  acc: *mut ::InterpAccel)
                                  -> c_double;
    pub fn gsl_spline_eval_deriv2_e(spline: *const gsl_spline,
                                    x: c_double,
                                    acc: *mut ::InterpAccel,
                                    d2: *mut c_double)
                                    -> c_int;
    pub fn gsl_spline_eval_integ(spline: *const gsl_spline,
                                 a: c_double,
                                 b: c_double,
                                 acc: *mut ::InterpAccel)
                                 -> c_double;
    pub fn gsl_spline_eval_integ_e(spline: *const gsl_spline,
                                   a: c_double,
                                   b: c_double,
                                   acc: *mut ::InterpAccel,
                                   result: *mut c_double)
                                   -> c_int;
    pub fn gsl_min_test_interval(x_lower: c_double,
                                 x_upper: c_double,
                                 epsabs: c_double,
                                 epsrel: c_double)
                                 -> c_int;
    
    // N-tuples
    // Creating ntuples
    pub fn gsl_ntuple_create(filename: *mut c_char,
                             ntuple_data: *mut c_void,
                             size: size_t)
                             -> *mut gsl_ntuple;
    // Opening an existing ntuple file
    pub fn gsl_ntuple_open(filename: *mut c_char,
                           ntuple_data: *mut c_void,
                           size: size_t)
                           -> *mut gsl_ntuple;
    // Writing ntuples
    pub fn gsl_ntuple_write(ntuple: *mut gsl_ntuple) -> c_int;
    pub fn gsl_ntuple_bookdata(ntuple: *mut gsl_ntuple) -> c_int;
    // Reading ntuples
    pub fn gsl_ntuple_read(ntuple: *mut gsl_ntuple) -> c_int;
    // Closing an ntuple file
    pub fn gsl_ntuple_close(ntuple: *mut gsl_ntuple) -> c_int;

    // Multisets
    // Multiset allocation
    pub fn gsl_multiset_alloc(n: size_t, k: size_t) -> *mut gsl_multiset;
    pub fn gsl_multiset_calloc(n: size_t, k: size_t) -> *mut gsl_multiset;
    pub fn gsl_multiset_init_first(c: *mut gsl_multiset);
    pub fn gsl_multiset_init_last(c: *mut gsl_multiset);
    pub fn gsl_multiset_free(c: *mut gsl_multiset);
    pub fn gsl_multiset_memcpy(dest: *mut gsl_multiset, src: *const gsl_multiset) -> c_int;
    // Accessing multiset elements
    pub fn gsl_multiset_get(c: *const gsl_multiset, i: size_t) -> size_t;
    // Multiset properties
    pub fn gsl_multiset_n(c: *const gsl_multiset) -> size_t;
    pub fn gsl_multiset_k(c: *const gsl_multiset) -> size_t;
    //pub fn gsl_multiset_data(c: *const gsl_multiset) -> *mut size_t;
    pub fn gsl_multiset_valid(c: *mut gsl_multiset) -> c_int;
    // Multiset functions
    pub fn gsl_multiset_next(c: *mut gsl_multiset) -> c_int;
    pub fn gsl_multiset_prev(c: *mut gsl_multiset) -> c_int;

    // Ordinary Differential Equations
    // Stepping Functions
    pub fn gsl_odeiv2_step_alloc(t: *const gsl_odeiv2_step_type,
                                 dim: size_t)
                                 -> *mut gsl_odeiv2_step;
    pub fn gsl_odeiv2_step_reset(s: *mut gsl_odeiv2_step) -> c_int;
    pub fn gsl_odeiv2_step_free(s: *mut gsl_odeiv2_step);
    pub fn gsl_odeiv2_step_name(s: *mut gsl_odeiv2_step) -> *const c_char;
    pub fn gsl_odeiv2_step_order(s: *const gsl_odeiv2_step) -> c_uint;
    pub fn gsl_odeiv2_step_set_driver(s: *mut gsl_odeiv2_step,
                                      d: *const gsl_odeiv2_driver)
                                      -> c_int;
    pub fn gsl_odeiv2_step_apply(s: *mut gsl_odeiv2_step,
                                 t: c_double,
                                 h: c_double,
                                 y: *mut c_double,
                                 yerr: *mut c_double,
                                 dydt_in: *const c_double,
                                 dydt_out: *mut c_double,
                                 sys: *const gsl_odeiv2_system)
                                 -> c_int;
    // Adaptive Step-size Control
    pub fn gsl_odeiv2_control_standard_new(eps_abs: c_double,
                                           eps_rel: c_double,
                                           a_y: c_double,
                                           a_dydt: c_double)
                                           -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_y_new(eps_abs: c_double,
                                    eps_rel: c_double)
                                    -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_yp_new(eps_abs: c_double,
                                     eps_rel: c_double)
                                     -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_scaled_new(eps_abs: c_double,
                                         eps_rel: c_double,
                                         a_y: c_double,
                                         a_dydt: c_double,
                                         scale_abs: *const c_double,
                                         dim: size_t)
                                         -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_alloc(t: *const gsl_odeiv2_control_type) -> *mut gsl_odeiv2_control;
    pub fn gsl_odeiv2_control_init(c: *mut gsl_odeiv2_control,
                                   eps_abs: c_double,
                                   eps_rel: c_double,
                                   a_y: c_double,
                                   a_dydt: c_double)
                                   -> c_int;
    pub fn gsl_odeiv2_control_free(c: *mut gsl_odeiv2_control);
    pub fn gsl_odeiv2_control_hadjust(c: *mut gsl_odeiv2_control,
                                      s: *mut gsl_odeiv2_step,
                                      y: *const c_double,
                                      yerr: *const c_double,
                                      dydt: *const c_double,
                                      h: *mut c_double)
                                      -> c_int;
    pub fn gsl_odeiv2_control_name(c: *const gsl_odeiv2_control) -> *const c_char;
    pub fn gsl_odeiv2_control_errlevel(c: *mut gsl_odeiv2_control,
                                       y: c_double,
                                       dydt: c_double,
                                       h: c_double,
                                       ind: size_t,
                                       errlev: *mut c_double)
                                       -> c_int;
    pub fn gsl_odeiv2_control_set_driver(c: *mut gsl_odeiv2_control,
                                         d: *const gsl_odeiv2_driver)
                                         -> c_int;
    // Evolution
    pub fn gsl_odeiv2_evolve_alloc(dim: size_t) -> *mut gsl_odeiv2_evolve;
    pub fn gsl_odeiv2_evolve_apply(e: *mut gsl_odeiv2_evolve,
                                   con: *mut gsl_odeiv2_control,
                                   step: *mut gsl_odeiv2_step,
                                   sys: *const gsl_odeiv2_system,
                                   t: *mut c_double,
                                   t1: c_double,
                                   h: *mut c_double,
                                   y: *mut c_double)
                                   -> c_int;
    pub fn gsl_odeiv2_evolve_apply_fixed_step(e: *mut gsl_odeiv2_evolve,
                                              con: *mut gsl_odeiv2_control,
                                              step: *mut gsl_odeiv2_step,
                                              sys: *const gsl_odeiv2_system,
                                              t: *mut c_double,
                                              h: c_double,
                                              y: *mut c_double)
                                              -> c_int;
    pub fn gsl_odeiv2_evolve_reset(e: *mut gsl_odeiv2_evolve) -> c_int;
    pub fn gsl_odeiv2_evolve_free(e: *mut gsl_odeiv2_evolve);
    pub fn gsl_odeiv2_evolve_set_driver(e: *mut gsl_odeiv2_evolve,
                                        d: *const gsl_odeiv2_driver)
                                        -> c_int;
    // Driver
    pub fn gsl_odeiv2_driver_alloc_y_new(sys: *const gsl_odeiv2_system,
                                         t: *const gsl_odeiv2_step_type,
                                         hstart: c_double,
                                         epsabs: c_double,
                                         epsrel: c_double)
                                         -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_alloc_yp_new(sys: *const gsl_odeiv2_system,
                                          t: *const gsl_odeiv2_step_type,
                                          hstart: c_double,
                                          epsabs: c_double,
                                          epsrel: c_double)
                                          -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_alloc_standard_new(sys: *const gsl_odeiv2_system,
                                                t: *const gsl_odeiv2_step_type,
                                                hstart: c_double,
                                                epsabs: c_double,
                                                epsrel: c_double,
                                                a_y: c_double,
                                                a_dydt: c_double)
                                                -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_alloc_scaled_new(sys: *const gsl_odeiv2_system,
                                              t: *const gsl_odeiv2_step_type,
                                              hstart: c_double,
                                              epsabs: c_double,
                                              epsrel: c_double,
                                              a_y: c_double,
                                              a_dydt: c_double,
                                              scale_abs: *const c_double)
                                              -> *mut gsl_odeiv2_driver;
    pub fn gsl_odeiv2_driver_set_hmin(d: *mut gsl_odeiv2_driver, hmin: c_double) -> c_int;
    pub fn gsl_odeiv2_driver_set_hmax(d: *mut gsl_odeiv2_driver, hmax: c_double) -> c_int;
    pub fn gsl_odeiv2_driver_set_nmax(d: *mut gsl_odeiv2_driver, nmax: usize) -> c_int;
    pub fn gsl_odeiv2_driver_apply(d: *mut gsl_odeiv2_driver,
                                   t: *mut c_double,
                                   t1: c_double,
                                   y: *mut c_double)
                                   -> c_int;
    pub fn gsl_odeiv2_driver_apply_fixed_step(d: *mut gsl_odeiv2_driver,
                                              t: *mut c_double,
                                              h: c_double,
                                              n: usize,
                                              y: *mut c_double)
                                              -> c_int;
    pub fn gsl_odeiv2_driver_reset(d: *mut gsl_odeiv2_driver) -> c_int;
    pub fn gsl_odeiv2_driver_reset_hstart(d: *mut gsl_odeiv2_driver,
                                          hstart: c_double)
                                          -> c_int;
    pub fn gsl_odeiv2_driver_free(d: *mut gsl_odeiv2_driver);

    // Quasi-Random Sequences
    // Quasi-random number generator initialization
    pub fn gsl_qrng_alloc(t: *const gsl_qrng_type, d: c_uint) -> *mut gsl_qrng;
    pub fn gsl_qrng_free(q: *mut gsl_qrng);
    pub fn gsl_qrng_init(q: *mut gsl_qrng);
    // Sampling from a quasi-random number generator
    pub fn gsl_qrng_get(q: *const gsl_qrng, x: *mut c_double) -> c_int;
    // Auxiliary quasi-random number generator functions
    pub fn gsl_qrng_name(q: *const gsl_qrng) -> *const c_char;
    pub fn gsl_qrng_size(q: *const gsl_qrng) -> size_t;
    pub fn gsl_qrng_state(q: *const gsl_qrng) -> *mut c_void;
    // Saving and resorting quasi-random number generator state
    pub fn gsl_qrng_memcpy(dest: *mut gsl_qrng, src: *const gsl_qrng) -> c_int;
    pub fn gsl_qrng_clone(q: *const gsl_qrng) -> *mut gsl_qrng;

    // Statistics
    // Mean, Standard Deviation and Variance
    pub fn gsl_stats_mean(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_variance(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_variance_m(data: *const c_double,
                                stride: size_t,
                                n: size_t,
                                mean: c_double)
                                -> c_double;
    pub fn gsl_stats_sd(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_sd_m(data: *const c_double,
                          stride: size_t,
                          n: size_t,
                          mean: c_double)
                          -> c_double;
    pub fn gsl_stats_tss(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_tss_m(data: *const c_double,
                           stride: size_t,
                           n: size_t,
                           mean: c_double)
                           -> c_double;
    pub fn gsl_stats_variance_with_fixed_mean(data: *const c_double,
                                              stride: size_t,
                                              n: size_t,
                                              mean: c_double)
                                              -> c_double;
    pub fn gsl_stats_sd_with_fixed_mean(data: *const c_double,
                                        stride: size_t,
                                        n: size_t,
                                        mean: c_double)
                                        -> c_double;
    // Absolute deviation
    pub fn gsl_stats_absdev(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_absdev_m(data: *const c_double,
                              stride: size_t,
                              n: size_t,
                              mean: c_double)
                              -> c_double;
    // Higher moments (skewness and kurtosis)
    pub fn gsl_stats_skew(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_skew_m_sd(data: *const c_double,
                               stride: size_t,
                               n: size_t,
                               mean: c_double,
                               sd: c_double)
                               -> c_double;
    pub fn gsl_stats_kurtosis(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_kurtosis_m_sd(data: *const c_double,
                                   stride: size_t,
                                   n: size_t,
                                   mean: c_double,
                                   sd: c_double)
                                   -> c_double;
    // Autocorrelation
    pub fn gsl_stats_lag1_autocorrelation(data: *const c_double,
                                          stride: size_t,
                                          n: size_t)
                                          -> c_double;
    pub fn gsl_stats_lag1_autocorrelation_m(data: *const c_double,
                                            stride: size_t,
                                            n: size_t,
                                            mean: c_double)
                                            -> c_double;
    // Covariance
    pub fn gsl_stats_covariance(data1: *const c_double,
                                stride1: size_t,
                                data2: *const c_double,
                                stride2: size_t,
                                n: size_t)
                                -> c_double;
    pub fn gsl_stats_covariance_m(data1: *const c_double,
                                  stride1: size_t,
                                  data2: *const c_double,
                                  stride2: size_t,
                                  n: size_t,
                                  mean1: c_double,
                                  mean2: c_double)
                                  -> c_double;
    // Correlation
    pub fn gsl_stats_correlation(data1: *const c_double,
                                 stride1: size_t,
                                 data2: *const c_double,
                                 stride2: size_t,
                                 n: size_t)
                                 -> c_double;
    pub fn gsl_stats_spearman(data1: *const c_double,
                              stride1: size_t,
                              data2: *const c_double,
                              stride2: size_t,
                              n: size_t,
                              work: *mut c_double)
                              -> c_double;
    // Weighted Samples
    pub fn gsl_stats_wmean(w: *const c_double,
                           wstride: size_t,
                           data: *const c_double,
                           stride: size_t,
                           n: size_t)
                           -> c_double;
    pub fn gsl_stats_wvariance(w: *const c_double,
                               wstride: size_t,
                               data: *const c_double,
                               stride: size_t,
                               n: size_t)
                               -> c_double;
    pub fn gsl_stats_wvariance_m(w: *const c_double,
                                 wstride: size_t,
                                 data: *const c_double,
                                 stride: size_t,
                                 n: size_t,
                                 wmean: c_double)
                                 -> c_double;
    pub fn gsl_stats_wsd(w: *const c_double,
                         wstride: size_t,
                         data: *const c_double,
                         stride: size_t,
                         n: size_t)
                         -> c_double;
    pub fn gsl_stats_wsd_m(w: *const c_double,
                           wstride: size_t,
                           data: *const c_double,
                           stride: size_t,
                           n: size_t,
                           wmean: c_double)
                           -> c_double;
    pub fn gsl_stats_wvariance_with_fixed_mean(w: *const c_double,
                                               wstride: size_t,
                                               data: *const c_double,
                                               stride: size_t,
                                               n: size_t,
                                               wmean: c_double)
                                               -> c_double;
    pub fn gsl_stats_wsd_with_fixed_mean(w: *const c_double,
                                         wstride: size_t,
                                         data: *const c_double,
                                         stride: size_t,
                                         n: size_t,
                                         wmean: c_double)
                                         -> c_double;
    pub fn gsl_stats_wtss(w: *const c_double,
                          wstride: size_t,
                          data: *const c_double,
                          stride: size_t,
                          n: size_t)
                          -> c_double;
    pub fn gsl_stats_wtss_m(w: *const c_double,
                            wstride: size_t,
                            data: *const c_double,
                            stride: size_t,
                            n: size_t,
                            wmean: c_double)
                            -> c_double;
    pub fn gsl_stats_wabsdev(w: *const c_double,
                             wstride: size_t,
                             data: *const c_double,
                             stride: size_t,
                             n: size_t)
                             -> c_double;
    pub fn gsl_stats_wabsdev_m(w: *const c_double,
                               wstride: size_t,
                               data: *const c_double,
                               stride: size_t,
                               n: size_t,
                               wmean: c_double)
                               -> c_double;
    pub fn gsl_stats_wskew(w: *const c_double,
                           wstride: size_t,
                           data: *const c_double,
                           stride: size_t,
                           n: size_t)
                           -> c_double;
    pub fn gsl_stats_wskew_m_sd(w: *const c_double,
                                wstride: size_t,
                                data: *const c_double,
                                stride: size_t,
                                n: size_t,
                                wmean: c_double,
                                wsd: c_double)
                                -> c_double;
    pub fn gsl_stats_wkurtosis(w: *const c_double,
                               wstride: size_t,
                               data: *const c_double,
                               stride: size_t,
                               n: size_t)
                               -> c_double;
    pub fn gsl_stats_wkurtosis_m_sd(w: *const c_double,
                                    wstride: size_t,
                                    data: *const c_double,
                                    stride: size_t,
                                    n: size_t,
                                    wmean: c_double,
                                    wsd: c_double)
                                    -> c_double;
    // Maximum and Minimum values
    pub fn gsl_stats_max(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_min(data: *const c_double, stride: size_t, n: size_t) -> c_double;
    pub fn gsl_stats_minmax(min: *mut c_double,
                            max: *mut c_double,
                            data: *const c_double,
                            stride: size_t,
                            n: size_t);
    pub fn gsl_stats_max_index(data: *const c_double, stride: size_t, n: size_t) -> size_t;
    pub fn gsl_stats_min_index(data: *const c_double, stride: size_t, n: size_t) -> size_t;
    pub fn gsl_stats_minmax_index(min: *mut size_t,
                                  max: *mut size_t,
                                  data: *const c_double,
                                  stride: size_t,
                                  n: size_t);
    // Median and Percentiles
    pub fn gsl_stats_median_from_sorted_data(data: *const c_double,
                                             stride: size_t,
                                             n: size_t)
                                             -> c_double;
    pub fn gsl_stats_quantile_from_sorted_data(data: *const c_double,
                                               stride: size_t,
                                               n: size_t,
                                               f: c_double)
                                               -> c_double;

    // Series Acceleration
    // Acceleration functions
    pub fn gsl_sum_levin_u_alloc(n: size_t) -> *mut gsl_sum_levin_u_workspace;
    pub fn gsl_sum_levin_u_free(w: *mut gsl_sum_levin_u_workspace);
    pub fn gsl_sum_levin_u_accel(array: *const c_double,
                                 array_size: size_t,
                                 w: *mut gsl_sum_levin_u_workspace,
                                 sum_accel: *mut c_double,
                                 abserr: *mut c_double)
                                 -> c_int;
    // Acceleration functions without error estimation
    pub fn gsl_sum_levin_utrunc_alloc(n: size_t) -> *mut gsl_sum_levin_utrunc_workspace;
    pub fn gsl_sum_levin_utrunc_free(w: *mut gsl_sum_levin_utrunc_workspace);
    pub fn gsl_sum_levin_utrunc_accel(array: *const c_double,
                                      array_size: size_t,
                                      w: *mut gsl_sum_levin_utrunc_workspace,
                                      sum_accel: *mut c_double,
                                      abserr_trunc: *mut c_double)
                                      -> c_int;

    // Wavelet Transforms
    // Initialization
    pub fn gsl_wavelet_alloc(t: *const gsl_wavelet_type, k: size_t) -> *mut gsl_wavelet;
    pub fn gsl_wavelet_name(w: *const gsl_wavelet) -> *const c_char;
    pub fn gsl_wavelet_free(w: *mut gsl_wavelet);
    pub fn gsl_wavelet_workspace_alloc(n: size_t) -> *mut gsl_wavelet_workspace;
    pub fn gsl_wavelet_workspace_free(w: *mut gsl_wavelet_workspace);
    // Wavelet transforms in one dimension
    pub fn gsl_wavelet_transform(w: *const gsl_wavelet,
                                 data: *mut c_double,
                                 stride: size_t,
                                 n: size_t,
                                 dir: c_int,
                                 work: *mut gsl_wavelet_workspace)
                                 -> c_int;
    pub fn gsl_wavelet_transform_forward(w: *const gsl_wavelet,
                                         data: *mut c_double,
                                         stride: size_t,
                                         n: size_t,
                                         work: *mut gsl_wavelet_workspace)
                                         -> c_int;
    pub fn gsl_wavelet_transform_inverse(w: *const gsl_wavelet,
                                         data: *mut c_double,
                                         stride: size_t,
                                         n: size_t,
                                         work: *mut gsl_wavelet_workspace)
                                         -> c_int;
    // Wavelet transforms in two dimension
    pub fn gsl_wavelet2d_transform(w: *const gsl_wavelet,
                                   data: *mut c_double,
                                   tda: size_t,
                                   size1: size_t,
                                   size2: size_t,
                                   dir: c_int,
                                   work: *mut gsl_wavelet_workspace)
                                   -> c_int;
    pub fn gsl_wavelet2d_transform_forward(w: *const gsl_wavelet,
                                           data: *mut c_double,
                                           tda: size_t,
                                           size1: size_t,
                                           size2: size_t,
                                           work: *mut gsl_wavelet_workspace)
                                           -> c_int;
    pub fn gsl_wavelet2d_transform_inverse(w: *const gsl_wavelet,
                                           data: *mut c_double,
                                           tda: size_t,
                                           size1: size_t,
                                           size2: size_t,
                                           work: *mut gsl_wavelet_workspace)
                                           -> c_int;
    pub fn gsl_wavelet2d_transform_matrix(w: *const gsl_wavelet,
                                          m: *mut gsl_matrix,
                                          dir: c_int,
                                          work: *mut gsl_wavelet_workspace)
                                          -> c_int;
    pub fn gsl_wavelet2d_transform_matrix_forward(w: *const gsl_wavelet,
                                                  m: *mut gsl_matrix,
                                                  work: *mut gsl_wavelet_workspace)
                                                  -> c_int;
    pub fn gsl_wavelet2d_transform_matrix_inverse(w: *const gsl_wavelet,
                                                  m: *mut gsl_matrix,
                                                  work: *mut gsl_wavelet_workspace)
                                                  -> c_int;
    pub fn gsl_wavelet2d_nstransform(w: *const gsl_wavelet,
                                     data: *mut c_double,
                                     tda: size_t,
                                     size1: size_t,
                                     size2: size_t,
                                     dir: c_int,
                                     work: *mut gsl_wavelet_workspace)
                                     -> c_int;
    pub fn gsl_wavelet2d_nstransform_forward(w: *const gsl_wavelet,
                                             data: *mut c_double,
                                             tda: size_t,
                                             size1: size_t,
                                             size2: size_t,
                                             work: *mut gsl_wavelet_workspace)
                                             -> c_int;
    pub fn gsl_wavelet2d_nstransform_inverse(w: *const gsl_wavelet,
                                             data: *mut c_double,
                                             tda: size_t,
                                             size1: size_t,
                                             size2: size_t,
                                             work: *mut gsl_wavelet_workspace)
                                             -> c_int;
    pub fn gsl_wavelet2d_nstransform_matrix(w: *const gsl_wavelet,
                                            m: *mut gsl_matrix,
                                            dir: c_int,
                                            work: *mut gsl_wavelet_workspace)
                                            -> c_int;
    pub fn gsl_wavelet2d_nstransform_matrix_forward(w: *const gsl_wavelet,
                                                    m: *mut gsl_matrix,
                                                    work: *mut gsl_wavelet_workspace)
                                                    -> c_int;
    pub fn gsl_wavelet2d_nstransform_matrix_inverse(w: *const gsl_wavelet,
                                                    m: *mut gsl_matrix,
                                                    work: *mut gsl_wavelet_workspace)
                                                    -> c_int;
}

#[repr(C)]
pub struct gsl_sf_result {
    pub val: c_double,
    pub err: c_double,
}

#[repr(C)]
pub struct gsl_sf_result_e10 {
    pub val: c_double,
    pub err: c_double,
    pub e10: c_int,
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
    pub wmat: *mut gsl_eigen_symmv_workspace,
}

#[repr(C)]
pub struct gsl_bspline_workspace {
    pub k: size_t, // spline order
    pub km1: size_t, // k - 1 (polynomial order)
    pub l: size_t, // number of polynomial pieces on interval
    pub nbreak: size_t, // number of breakpoints (l + 1)
    pub n: size_t, // number of bspline basis functions (l + k - 1)
    pub knots: *mut gsl_vector, // knots vector
    pub deltal: *mut gsl_vector, // left delta
    pub deltar: *mut gsl_vector, // right delta
    pub B: *mut gsl_vector, // temporary spline results
}

#[repr(C)]
#[cfg(not(feature = "v2"))]
pub struct gsl_bspline_deriv_workspace {
    pub k: size_t, // spline order
    pub A: *mut gsl_matrix, // work matrix
    pub dB: *mut gsl_matrix, // temporary derivative results
}

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
    pub get_double: rng_get_double,
}

#[repr(C)]
pub struct gsl_rng {
    pub _type: *const gsl_rng_type,
    pub state: *mut c_void,
}

#[repr(C)]
pub struct gsl_ran_discrete_t {
    pub K: size_t,
    pub A: *mut size_t,
    pub F: *mut c_double,
}

#[repr(C)]
pub struct gsl_permutation {
    pub size: size_t,
    pub data: *mut size_t,
}

#[repr(C)]
pub struct gsl_cheb_series {
    pub c: *mut c_double, // coefficients
    pub order: c_int, // order of expansion
    pub a: c_double, // lower interval point
    pub b: c_double, // upper interval point
    pub order_sp: size_t,
    pub f: *mut c_double,
}

#[repr(C)]
pub struct gsl_combination {
    pub n: size_t,
    pub k: size_t,
    pub data: *mut size_t,
}

#[repr(C)]
pub struct gsl_poly_complex_workspace {
    pub nc: size_t,
    pub matrix: *mut c_double,
}

#[repr(C)]
pub struct gsl_dht {
    pub size: size_t, // size of the sample arrays to be transformed
    pub nu: c_double, // Bessel function order
    pub xmax: c_double, // the upper limit to the x-sampling domain
    pub kmax: c_double, // the upper limit to the k-sampling domain
    pub j: *mut c_double, // array of computed J_nu zeros, j_{nu,s} = j[s]
    pub Jjj: *mut c_double, // transform numerator, J_nu(j_i j_m / j_N)
    pub J2: *mut c_double, // transform denominator, J_{nu+1}^2(j_m)
}

#[repr(C)]
pub struct gsl_fft_complex_wavetable {
    pub n: size_t,
    pub nf: size_t,
    pub factor: [size_t; 64],
    pub twiddle: [*mut gsl_complex; 64],
    pub trig: *mut gsl_complex,
}

#[repr(C)]
pub struct gsl_fft_complex_workspace {
    pub n: size_t,
    pub scratch: *mut c_double,
}

#[repr(C)]
pub struct gsl_histogram {
    pub n: size_t, // This is the number of histogram bins
    pub range: *mut c_double, // The ranges of the bins are stored in an array of n+1 elements pointed to by range.
    pub bin: *mut c_double, /* The counts for each bin are stored in an array of n elements pointed to by bin. The bins are floating-point numbers, so you can increment them by non-integer values if necessary.

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
    pub sum: *mut c_double, // The cumulative probability for the bins is stored in an array of n elements pointed to by sum.
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
    pub sum: *mut c_double, // The cumulative probability for the bins is stored in an array of nx*ny elements pointed to by sum.
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
    pub level: *mut size_t,
}

#[repr(C)]
pub struct extrapolation_table {
    pub n: size_t,
    pub rlist2: [c_double; 52],
    pub nres: size_t,
    pub res3la: [c_double; 3],
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
    pub rh: [c_double; 25],
}

#[repr(C)]
pub struct gsl_integration_qawo_table {
    pub n: size_t,
    pub omega: c_double,
    pub L: c_double,
    pub par: c_double,
    pub sine: c_int,
    pub chebmo: *mut c_double,
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
    pub ndiv: c_int,
}


/* The workspace is just a collection of intervals */
#[repr(C)]
pub struct gsl_integration_cquad_workspace {
    pub size: size_t,
    pub ivals: *mut gsl_integration_cquad_ival,
    pub heap: *mut size_t,
}

/* Workspace for fixed-order Gauss-Legendre integration */
#[repr(C)]
pub struct gsl_integration_glfixed_table {
    pub n: size_t, /* number of points */
    pub x: *mut c_double, /* Gauss abscissae/points */
    pub w: *mut c_double, /* Gauss weights for each abscissae */
    pub precomputed: c_int, /* high precision abscissae/weights precomputed? */
}

/* interpolation object type */
#[repr(C)]
pub struct gsl_interp_type {
    pub name: *const c_char,
    pub min_size: c_uint,
    pub alloc: Option<extern "C" fn(size_t) -> *mut c_void>,
    pub init: Option<extern "C" fn(*mut c_void, *const c_double, *const c_double, size_t)
                                   -> c_int>,
    pub eval: Option<extern "C" fn(*const c_void,
                                   *const c_double,
                                   *const c_double,
                                   size_t,
                                   c_double,
                                   *mut ::InterpAccel,
                                   *mut c_double)
                                   -> c_int>,
    pub eval_deriv: Option<extern "C" fn(*const c_void,
                                         *const c_double,
                                         *const c_double,
                                         size_t,
                                         c_double,
                                         *mut ::InterpAccel,
                                         *mut c_double)
                                         -> c_int>,
    pub eval_deriv2: Option<extern "C" fn(*const c_void,
                                          *const c_double,
                                          *const c_double,
                                          size_t,
                                          c_double,
                                          *mut ::InterpAccel,
                                          *mut c_double)
                                          -> c_int>,
    pub eval_integ: Option<extern "C" fn(*const c_void,
                                         *const c_double,
                                         *const c_double,
                                         size_t,
                                         c_double,
                                         *mut ::InterpAccel,
                                         c_double,
                                         c_double,
                                         *mut c_double)
                                         -> c_int>,
    pub free: Option<extern "C" fn(*mut c_void)>,
}

/* general interpolation object */
#[repr(C)]
pub struct gsl_interp {
    pub _type: *const gsl_interp_type,
    pub xmin: c_double,
    pub xmax: c_double,
    pub size: size_t,
    pub state: *mut c_void,
}

/* general interpolation object */
#[repr(C)]
pub struct gsl_spline {
    pub interp: *mut gsl_interp,
    pub x: *mut c_double,
    pub y: *mut c_double,
    pub size: size_t,
}

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
    pub data: *mut size_t,
}

#[repr(C)]
pub struct gsl_odeiv2_system {
    pub function: extern "C" fn(t: c_double, *const c_double, *mut c_double, *mut c_void)
                                -> c_int,
    pub jacobian: Option<extern "C" fn(t: c_double,
                                       *const c_double,
                                       *mut c_double,
                                       *mut c_double,
                                       *mut c_void)
                                       -> c_int>,
    pub dimension: usize,
    pub params: *mut c_void,
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
    pub nmax: usize,
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
    pub driver: *const gsl_odeiv2_driver,
}

#[repr(C)]
pub struct gsl_odeiv2_control {
    pub type_: *const gsl_odeiv2_control_type,
    pub state: *mut c_void,
}

#[repr(C)]
pub struct gsl_odeiv2_control_type {
    pub name: *const c_char,
    pub alloc: fn() -> *mut c_void,
    pub init: fn(state: *mut c_void,
                 eps_abs: c_double,
                 eps_rel: c_double,
                 a_y: c_double,
                 a_dydt: c_double)
                 -> c_int,
    pub hadjust: fn(state: *mut c_void,
                    dim: size_t,
                    ord: c_uint,
                    y: *const c_double,
                    yerr: *const c_double,
                    yp: *const c_double,
                    h: *mut c_double)
                    -> c_int,
    pub errlevel: fn(state: *mut c_void,
                     y: c_double,
                     dydt: c_double,
                     h: c_double,
                     ind: size_t,
                     errlev: *mut c_double)
                     -> c_int,
    pub set_driver: fn(state: *mut c_void, d: *const gsl_odeiv2_driver) -> c_int,
    pub free: fn(state: *mut c_void),
}

#[repr(C)]
pub struct gsl_odeiv2_step {
    pub type_: *const gsl_odeiv2_step_type,
    pub dimension: size_t,
    pub state: *mut c_void,
}

#[repr(C)]
pub struct gsl_odeiv2_step_type {
    pub name: *const c_char,
    pub can_use_dydt_in: c_int,
    pub gives_exact_dydt_out: c_int,
    pub alloc: fn(dim: size_t) -> *mut c_void,
    pub apply: fn(state: *mut c_void,
                  dim: size_t,
                  t: c_double,
                  h: c_double,
                  y: *mut c_double,
                  yerr: *mut c_double,
                  dydt_in: *const c_double,
                  dydt_out: *mut c_double,
                  dydt: *const gsl_odeiv2_system)
                  -> c_int,
    pub set_driver: fn(state: *mut c_void, d: *const gsl_odeiv2_driver) -> c_int,
    pub reset: fn(state: *mut c_void, dim: size_t) -> c_int,
    pub order: fn(state: *mut c_void) -> c_uint,
    pub free: fn(state: *mut c_void),
}

// Structure describing a generator instance of a specified type, with generator-specific state info and dimension-specific info.
#[repr(C)]
pub struct gsl_qrng {
    pub type_: *const gsl_qrng_type,
    pub dimension: c_uint,
    pub state_size: size_t,
    pub state: *mut c_void,
}

// Structure describing a type of generator.
#[repr(C)]
pub struct gsl_qrng_type {
    pub name: *const c_char,
    pub max_dimension: c_uint,
    pub state_size: Option<extern "C" fn(dimension: c_uint) -> size_t>,
    pub init_state: Option<extern "C" fn(state: *mut c_void, dimension: c_uint) -> c_int>,
    pub get: Option<extern "C" fn(state: *mut c_void, dimension: c_uint, x: *mut c_double)
                                  -> c_int>,
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
    pub dsum: *mut c_double,
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
    pub dsum: *mut c_double,
}

#[repr(C)]
pub struct gsl_wavelet_workspace {
    pub scratch: *mut c_double,
    pub n: size_t,
}

#[repr(C)]
pub struct gsl_wavelet {
    pub type_: *const gsl_wavelet_type,
    pub h1: *const c_double,
    pub g1: *const c_double,
    pub h2: *const c_double,
    pub g2: *const c_double,
    pub nc: size_t,
    pub offset: size_t,
}

#[repr(C)]
pub struct gsl_wavelet_type {
    pub name: *const c_char,
    pub init: Option<extern "C" fn(h1: *const *const c_double,
                                   g1: *const *const c_double,
                                   h2: *const *const c_double,
                                   g2: *const *const c_double,
                                   nc: *mut size_t,
                                   offset: *mut size_t,
                                   member: size_t)
                                   -> c_int>,
}
