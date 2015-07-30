//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use enums;
use std::mem::zeroed;

/// This routine computes the regular modified cylindrical Bessel function of zeroth order, I_0(x)
pub fn I0(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_I0(x) }
}

/// This routine computes the regular modified cylindrical Bessel function of zeroth order, I_0(x)
pub fn I0_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_I0_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular modified cylindrical Bessel function of first order, I_1(x).
pub fn I1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_I1(x) }
}

/// This routine computes the regular modified cylindrical Bessel function of first order, I_1(x).
pub fn I1_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_I1_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular modified cylindrical Bessel function of order n, I_n(x).
pub fn In(n: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_In(n, x) }
}

/// This routine computes the regular modified cylindrical Bessel function of order n, I_n(x).
pub fn In_e(n: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_In_e(n, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the regular modified cylindrical Bessel functions I_n(x) for n from nmin to nmax inclusive, storing the results in the array result_array.
/// The start of the range nmin must be positive or zero.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn In_array(nmin: i32, nmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_In_array(nmin, nmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the scaled regular modified cylindrical Bessel function of zeroth order \exp(-|x|) I_0(x).
pub fn I0_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_I0_scaled(x) }
}

/// This routine computes the scaled regular modified cylindrical Bessel function of zeroth order \exp(-|x|) I_0(x).
pub fn I0_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_I0_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled regular modified cylindrical Bessel function of first order \exp(-|x|) I_1(x).
pub fn I1_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_I1_scaled(x) }
}

/// This routine computes the scaled regular modified cylindrical Bessel function of first order \exp(-|x|) I_1(x).
pub fn I1_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_I1_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled regular modified cylindrical Bessel function of order n, \exp(-|x|) I_n(x)
pub fn In_scaled(n: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_In_scaled(n, x) }
}

/// This routine computes the scaled regular modified cylindrical Bessel function of order n, \exp(-|x|) I_n(x)
pub fn In_scaled_e(n: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_In_scaled_e(n, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the scaled regular cylindrical Bessel functions \exp(-|x|) I_n(x) for n from nmin to nmax inclusive, storing the results in the array result_array.
/// The start of the range nmin must be positive or zero.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn In_scaled_array(nmin: i32, nmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_In_scaled_array(nmin, nmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the scaled regular modified spherical Bessel function of zeroth order, \exp(-|x|) i_0(x).
pub fn i0_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_i0_scaled(x) }
}

/// This routine computes the scaled regular modified spherical Bessel function of zeroth order, \exp(-|x|) i_0(x).
pub fn i0_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_i0_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled regular modified spherical Bessel function of first order, \exp(-|x|) i_1(x).
pub fn i1_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_i1_scaled(x) }
}

/// This routine computes the scaled regular modified spherical Bessel function of first order, \exp(-|x|) i_1(x).
pub fn i1_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_i1_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled regular modified spherical Bessel function of second order, \exp(-|x|) i_2(x)
pub fn i2_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_i2_scaled(x) }
}

/// This routine computes the scaled regular modified spherical Bessel function of second order, \exp(-|x|) i_2(x)
pub fn i2_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_i2_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled regular modified spherical Bessel function of order l, \exp(-|x|) i_l(x)
pub fn il_scaled(l: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_il_scaled(l, x) }
}

/// This routine computes the scaled regular modified spherical Bessel function of order l, \exp(-|x|) i_l(x)
pub fn il_scaled_e(l: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_il_scaled_e(l, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the scaled regular modified cylindrical Bessel functions \exp(-|x|) i_l(x) for l from 0 to lmax inclusive for lmax >= 0, storing the results in the array result_array. The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn il_scaled_array(lmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_il_scaled_array(lmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the regular modified Bessel function of fractional order \nu, I_\nu(x) for x>0, \nu>0.
pub fn Inu(nu: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Inu(nu, x) }
}

/// This routine computes the regular modified Bessel function of fractional order \nu, I_\nu(x) for x>0, \nu>0.
pub fn Inu_e(nu: f64, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Inu_e(nu, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled regular modified Bessel function of fractional order \nu, \exp(-|x|)I_\nu(x) for x>0, \nu>0.
pub fn Inu_scaled(nu: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Inu_scaled(nu, x) }
}

/// This routine computes the scaled regular modified Bessel function of fractional order \nu, \exp(-|x|)I_\nu(x) for x>0, \nu>0.
pub fn Inu_scaled_e(nu: f64, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Inu_scaled_e(nu, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular cylindrical Bessel function of zeroth order, J_0(x).
pub fn J0(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_J0(x) }
}

/// This routine computes the regular cylindrical Bessel function of zeroth order, J_0(x).
pub fn J0_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_J0_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular cylindrical Bessel function of first order, J_1(x).
pub fn J1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_J1(x) }
}

/// This routine computes the regular cylindrical Bessel function of first order, J_1(x).
pub fn J1_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_J1_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular cylindrical Bessel function of order n, J_n(x).
pub fn Jn(n: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Jn(n, x) }
}

/// This routine computes the regular cylindrical Bessel function of order n, J_n(x).
pub fn Jn_e(n: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Jn_e(n, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the regular cylindrical Bessel functions J_n(x) for n from nmin to nmax inclusive, storing the results in the array result_array.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn Jn_array(nmin: i32, nmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_Jn_array(nmin, nmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the regular spherical Bessel function of zeroth order, j_0(x) = \sin(x)/x.
pub fn j0(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_j0(x) }
}

/// This routine computes the regular spherical Bessel function of zeroth order, j_0(x) = \sin(x)/x.
pub fn j0_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_j0_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular spherical Bessel function of first order, j_1(x) = (\sin(x)/x - \cos(x))/x.
pub fn j1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_j1(x) }
}

/// This routine computes the regular spherical Bessel function of first order, j_1(x) = (\sin(x)/x - \cos(x))/x.
pub fn j1_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_j1_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular spherical Bessel function of second order, j_2(x) = ((3/x^2 - 1)\sin(x) - 3\cos(x)/x)/x.
pub fn j2(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_j2(x) }
}

/// This routine computes the regular spherical Bessel function of second order, j_2(x) = ((3/x^2 - 1)\sin(x) - 3\cos(x)/x)/x.
pub fn j2_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_j2_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the regular spherical Bessel function of order l, j_l(x), for l >= 0 and x >= 0.
pub fn jl(l: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_jl(l, x) }
}

/// This routine computes the regular spherical Bessel function of order l, j_l(x), for l >= 0 and x >= 0.
pub fn jl_e(l: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_jl_e(l, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the regular spherical Bessel functions j_l(x) for l from 0 to lmax inclusive for lmax >= 0 and x >= 0, storing the results in the array result_array.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn jl_array(lmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_jl_array(lmax, x, result_array.as_mut_ptr()) }
}

/// This routine uses Steed’s method to compute the values of the regular spherical Bessel functions j_l(x) for l from 0 to lmax inclusive for lmax >= 0 and x >= 0, storing the results in the array result_array.
/// The Steed/Barnett algorithm is described in Comp. Phys. Comm. 21, 297 (1981). Steed’s method is more stable than the recurrence used in the other functions but is also slower.
pub fn jl_steed_array(lmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_jl_steed_array(lmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the regular cylindrical Bessel function of fractional order \nu, J_\nu(x).
pub fn Jnu(nu: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Jnu(nu, x) }
}

/// This routine computes the regular cylindrical Bessel function of fractional order \nu, J_\nu(x).
pub fn Jnu_e(nu: f64, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Jnu_e(nu, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This function computes the regular cylindrical Bessel function of fractional order \nu, J_\nu(x), evaluated at a series of x values. The array v of length size contains the x values.
/// They are assumed to be strictly ordered and positive. The array is over-written with the values of J_\nu(x_i).
pub fn sequence_Jnu(nu: f64, mode: ::Mode, v: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_sequence_Jnu_e(nu, mode, v.len() as i64, v.as_mut_ptr()) }
}

/// This routine computes the irregular modified cylindrical Bessel function of zeroth order, K_0(x), for x > 0.
pub fn K0(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_K0(x) }
}

/// This routine computes the irregular modified cylindrical Bessel function of zeroth order, K_0(x), for x > 0.
pub fn K0_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_K0_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular modified cylindrical Bessel function of first order, K_1(x), for x > 0.
pub fn K1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_K1(x) }
}

/// This routine computes the irregular modified cylindrical Bessel function of first order, K_1(x), for x > 0.
pub fn K1_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_K1_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular modified cylindrical Bessel function of order n, K_n(x), for x > 0.
pub fn Kn(n: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Kn(n, x) }
}

/// This routine computes the irregular modified cylindrical Bessel function of order n, K_n(x), for x > 0.
pub fn Kn_e(n: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Kn_e(n, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the irregular modified cylindrical Bessel functions K_n(x) for n from nmin to nmax inclusive, storing the results in the array result_array.
/// The start of the range nmin must be positive or zero. The domain of the function is x>0.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn Kn_array(nmin: i32, nmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_Kn_array(nmin, nmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the scaled irregular modified cylindrical Bessel function of zeroth order \exp(x) K_0(x) for x>0.
pub fn K0_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_K0_scaled(x) }
}

/// This routine computes the scaled irregular modified cylindrical Bessel function of zeroth order \exp(x) K_0(x) for x>0.
pub fn K0_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_K0_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled irregular modified cylindrical Bessel function of first order \exp(x) K_1(x) for x>0.
pub fn K1_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_K1_scaled(x) }
}

/// This routine computes the scaled irregular modified cylindrical Bessel function of first order \exp(x) K_1(x) for x>0.
pub fn K1_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_K1_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled irregular modified cylindrical Bessel function of order n, \exp(x) K_n(x), for x>0.
pub fn Kn_scaled(n: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Kn_scaled(n, x) }
}

/// This routine computes the scaled irregular modified cylindrical Bessel function of order n, \exp(x) K_n(x), for x>0.
pub fn Kn_scaled_e(n: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Kn_scaled_e(n, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the scaled irregular cylindrical Bessel functions \exp(x) K_n(x) for n from nmin to nmax inclusive, storing the results in the array result_array.
/// The start of the range nmin must be positive or zero. The domain of the function is x>0.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn Kn_scaled_array(nmin: i32, nmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_Kn_scaled_array(nmin, nmax, x, result_array.as_mut_ptr()) }
}

/// The irregular modified spherical Bessel functions k_l(x) are related to the irregular modified Bessel functions of fractional order, k_l(x) = \sqrt{\pi/(2x)} K_{l+1/2}(x).
/// This routine computes the scaled irregular modified spherical Bessel function of zeroth order, \exp(x) k_0(x), for x>0.
pub fn k0_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_k0_scaled(x) }
}

/// The irregular modified spherical Bessel functions k_l(x) are related to the irregular modified Bessel functions of fractional order, k_l(x) = \sqrt{\pi/(2x)} K_{l+1/2}(x).
/// This routine computes the scaled irregular modified spherical Bessel function of zeroth order, \exp(x) k_0(x), for x>0.
pub fn k0_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_k0_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled irregular modified spherical Bessel function of first order, \exp(x) k_1(x), for x>0.
pub fn k1_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_k1_scaled(x) }
}

/// This routine computes the scaled irregular modified spherical Bessel function of first order, \exp(x) k_1(x), for x>0.
pub fn k1_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_k1_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled irregular modified spherical Bessel function of second order, \exp(x) k_2(x), for x>0.
pub fn k2_scaled(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_k2_scaled(x) }
}

/// This routine computes the scaled irregular modified spherical Bessel function of second order, \exp(x) k_2(x), for x>0.
pub fn k2_scaled_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_k2_scaled_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the scaled irregular modified spherical Bessel function of order l, \exp(x) k_l(x), for x>0.
pub fn kl_scaled(l: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_kl_scaled(l, x) }
}

/// This routine computes the scaled irregular modified spherical Bessel function of order l, \exp(x) k_l(x), for x>0.
pub fn kl_scaled_e(l: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_kl_scaled_e(l, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the scaled irregular modified spherical Bessel functions \exp(x) k_l(x) for l from 0 to lmax inclusive for lmax >= 0 and x>0, storing the results in the array result_array.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn kl_scaled_array(lmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_kl_scaled_array(lmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the irregular modified Bessel function of fractional order \nu, K_\nu(x) for x>0, \nu>0.
pub fn Knu(nu: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Knu(nu, x) }
}

/// This routine computes the irregular modified Bessel function of fractional order \nu, K_\nu(x) for x>0, \nu>0.
pub fn Knu_e(nu: f64, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Knu_e(nu, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

// /// This routine computes the logarithm of the irregular modified Bessel function of fractional order \nu, \ln(K_\nu(x)) for x>0, \nu>0.
// pub fn InKnu(nu: f64, x: f64) -> f64 {
//     unsafe { ffi::gsl_sf_bessel_InKnu(nu, x) }
// }

// /// This routine computes the logarithm of the irregular modified Bessel function of fractional order \nu, \ln(K_\nu(x)) for x>0, \nu>0.
// pub fn InKnu_e(nu: f64, x: f64) -> (enums::value::Value, ::types::Result) {
//     let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
//     let ret = unsafe { ffi::gsl_sf_bessel_InKnu_e(nu, x, &mut result) };

//     (ret, ::types::Result{val: result.val, err: result.err})
// }

/// This routine computes the scaled irregular modified Bessel function of fractional order \nu, \exp(+|x|) K_\nu(x) for x>0, \nu>0.
pub fn Knu_scaled(nu: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Knu_scaled(nu, x) }
}

/// This routine computes the scaled irregular modified Bessel function of fractional order \nu, \exp(+|x|) K_\nu(x) for x>0, \nu>0.
pub fn Knu_scaled_e(nu: f64, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Knu_scaled_e(nu, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular cylindrical Bessel function of zeroth order, Y_0(x), for x>0.
pub fn Y0(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Y0(x) }
}

/// This routine computes the irregular cylindrical Bessel function of zeroth order, Y_0(x), for x>0.
pub fn Y0_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Y0_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular cylindrical Bessel function of first order, Y_1(x), for x>0.
pub fn Y1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Y1(x) }
}

/// This routine computes the irregular cylindrical Bessel function of first order, Y_1(x), for x>0.
pub fn Y1_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Y1_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular cylindrical Bessel function of order n, Y_n(x), for x>0.
pub fn Yn(n: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Yn(n, x) }
}

/// This routine computes the irregular cylindrical Bessel function of order n, Y_n(x), for x>0.
pub fn Yn_e(n: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Yn_e(n, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the irregular cylindrical Bessel functions Y_n(x) for n from nmin to nmax inclusive, storing the results in the array result_array.
/// The domain of the function is x>0.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn Yn_array(nmin: i32, nmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_Yn_array(nmin, nmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the irregular spherical Bessel function of zeroth order, y_0(x) = -\cos(x)/x.
pub fn y0(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_y0(x) }
}

/// This routine computes the irregular spherical Bessel function of zeroth order, y_0(x) = -\cos(x)/x.
pub fn y0_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_y0_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular spherical Bessel function of first order, y_1(x) = -(\cos(x)/x + \sin(x))/x.
pub fn y1(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_y1(x) }
}

/// This routine computes the irregular spherical Bessel function of first order, y_1(x) = -(\cos(x)/x + \sin(x))/x.
pub fn y1_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_y1_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular spherical Bessel function of second order, y_2(x) = (-3/x^3 + 1/x)\cos(x) - (3/x^2)\sin(x).
pub fn y2(x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_y2(x) }
}

/// This routine computes the irregular spherical Bessel function of second order, y_2(x) = (-3/x^3 + 1/x)\cos(x) - (3/x^2)\sin(x).
pub fn y2_e(x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_y2_e(x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the irregular spherical Bessel function of order l, y_l(x), for l >= 0.
pub fn yl(l: i32, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_yl(l, x) }
}

/// This routine computes the irregular spherical Bessel function of order l, y_l(x), for l >= 0.
pub fn yl_e(l: i32, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_yl_e(l, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the values of the irregular spherical Bessel functions y_l(x) for l from 0 to lmax inclusive for lmax >= 0, storing the results in the array result_array.
/// The values are computed using recurrence relations for efficiency, and therefore may differ slightly from the exact values.
pub fn yl_array(lmax: i32, x: f64, result_array: &mut [f64]) -> enums::value::Value {
    unsafe { ffi::gsl_sf_bessel_yl_array(lmax, x, result_array.as_mut_ptr()) }
}

/// This routine computes the irregular cylindrical Bessel function of fractional order \nu, Y_\nu(x).
pub fn Ynu(nu: f64, x: f64) -> f64 {
    unsafe { ffi::gsl_sf_bessel_Ynu(nu, x) }
}

/// This routine computes the irregular cylindrical Bessel function of fractional order \nu, Y_\nu(x).
pub fn Ynu_e(nu: f64, x: f64) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_Ynu_e(nu, x, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the location of the s-th positive zero of the Bessel function J_0(x).
pub fn zero_J0(s: u32) -> f64 {
    unsafe { ffi::gsl_sf_bessel_zero_J0(s) }
}

/// This routine computes the location of the s-th positive zero of the Bessel function J_0(x).
pub fn zero_J0_e(s: u32) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_zero_J0_e(s, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the location of the s-th positive zero of the Bessel function J_1(x).
pub fn zero_J1(s: u32) -> f64 {
    unsafe { ffi::gsl_sf_bessel_zero_J1(s) }
}

/// This routine computes the location of the s-th positive zero of the Bessel function J_1(x).
pub fn zero_J1_e(s: u32) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_zero_J1_e(s, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}

/// This routine computes the location of the s-th positive zero of the Bessel function J_\nu(x).
/// The current implementation does not support negative values of nu.
pub fn zero_Jnu(nu: f64, s: u32) -> f64 {
    unsafe { ffi::gsl_sf_bessel_zero_Jnu(nu, s) }
}

/// This routine computes the location of the s-th positive zero of the Bessel function J_\nu(x).
/// The current implementation does not support negative values of nu.
pub fn zero_Jnu_e(nu: f64, s: u32) -> (enums::value::Value, ::types::Result) {
    let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
    let ret = unsafe { ffi::gsl_sf_bessel_zero_Jnu_e(nu, s, &mut result) };

    (ret, ::types::Result{val: result.val, err: result.err})
}