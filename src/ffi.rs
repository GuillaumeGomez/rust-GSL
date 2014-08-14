//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use libc::{c_double, c_int, c_uint, c_float, c_void, size_t};
use enums;
use cblas;

pub type CBLAS_INDEX = c_uint;
pub type CBLAS_INDEX_t = CBLAS_INDEX;
pub type CBLAS_TRANSPOSE_t = cblas::Transpose;
pub type CBLAS_UPLO_t = cblas::Uplo;
pub type CBLAS_DIAG_t = cblas::Diag;
pub type CBLAS_SIDE_t = cblas::Side;

pub trait FFI<T> {
    fn wrap(r: *mut T) -> Self;
    fn unwrap(&Self) -> *mut T;
}

extern "C" {
    // Airy functions
    pub fn gsl_sf_airy_Ai(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Ai_scaled(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_scaled_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi_scaled(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_scaled_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
    // Derivatives of Airy Functions
    pub fn gsl_sf_airy_Ai_deriv(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi_deriv(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Ai_deriv_scaled(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Ai_deriv_scaled_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_airy_Bi_deriv_scaled(x: c_double, mode: enums::gsl_mode_t) -> c_double;
    pub fn gsl_sf_airy_Bi_deriv_scaled_e(x: c_double, mode: enums::gsl_mode_t, result: *mut gsl_sf_result) -> enums::Value;
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
    pub fn gsl_sf_bessel_sequence_Jnu_e(nu: c_double, mode: enums::gsl_mode_t, size: i64, v: *mut c_double) -> enums::Value;
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
    pub fn gsl_sf_bessel_InKnu(nu: c_double, x: c_double) -> c_double;
    pub fn gsl_sf_bessel_InKnu_e(nu: c_double, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
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
    pub fn gsl_sf_legendre_Plm_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_legendre_Plm_deriv_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double,
        result_deriv_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_legendre_sphPlm(l: c_int, m: c_int, x: c_double) -> c_double;
    pub fn gsl_sf_legendre_sphPlm_e(l: c_int, m: c_int, x: c_double, result: *mut gsl_sf_result) -> enums::Value;
    pub fn gsl_sf_legendre_sphPlm_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double) -> enums::Value;
    pub fn gsl_sf_legendre_sphPlm_deriv_array(lmax: c_int, m: c_int, x: c_double, result_array: *mut c_double,
        result_deriv_array: *mut c_double) -> enums::Value;
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
}

pub struct gsl_sf_result {
    pub val: c_double,
    pub err: c_double
}

pub struct gsl_sf_result_e10 {
    pub val: c_double,
    pub err: c_double,
    pub e10: c_int
}

pub struct gsl_vector_float {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_float,
    pub owner: c_int
}

pub struct gsl_block_float {
    pub size: size_t,
    pub data: *mut c_float
}

pub struct gsl_vector {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block,
    pub owner: c_int
}

pub struct gsl_block {
    pub size: size_t,
    pub data: *mut c_double
}

pub struct gsl_vector_complex_float {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_complex_float,
    pub owner: c_int
}

pub struct gsl_block_complex_float {
    pub size: size_t,
    pub data: *mut c_float
}

pub struct gsl_vector_complex {
    pub size: size_t,
    pub stride: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block_complex,
    pub owner: c_int
}

pub struct gsl_block_complex {
    pub size: size_t,
    pub data: *mut c_double
}

pub struct gsl_complex {
    pub data: [c_double, ..2]
}

pub struct gsl_complex_float {
    pub data: [c_float, ..2]
}

pub struct gsl_matrix {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block,
    pub owner: c_int
}

pub struct gsl_matrix_float {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_float,
    pub owner: c_int
}

pub struct gsl_matrix_complex {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_double,
    pub block: *mut gsl_block_complex,
    pub owner: c_int
}

pub struct gsl_matrix_complex_float {
    pub size1: size_t,
    pub size2: size_t,
    pub tda: size_t,
    pub data: *mut c_float,
    pub block: *mut gsl_block_complex_float,
    pub owner: c_int
}

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

pub struct gsl_eigen_symmv_workspace {
    pub size: size_t,
    pub d: *mut c_double,
    pub sd: *mut c_double,
    pub gc: *mut c_double,
    pub gs: *mut c_double
}