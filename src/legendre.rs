//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! The Legendre Functions and Legendre Polynomials are described in Abramowitz & Stegun, Chapter 8.

pub mod polynomials {
    use ffi;
    use enums;
    use std::mem::zeroed;

    /// This function evaluates the Legendre polynomials P_l(x) using explicit representations for l=1, 2, 3.
    pub fn legendre_P1(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_P1(x) }
    }

    /// This function evaluates the Legendre polynomials P_l(x) using explicit representations for l=1, 2, 3.
    pub fn legendre_P2(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_P2(x) }
    }

    /// This function evaluates the Legendre polynomials P_l(x) using explicit representations for l=1, 2, 3.
    pub fn legendre_P3(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_P3(x) }
    }

    /// This function evaluates the Legendre polynomials P_l(x) using explicit representations for l=1, 2, 3.
    pub fn legendre_P1_e(x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_P1_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function evaluates the Legendre polynomials P_l(x) using explicit representations for l=1, 2, 3.
    pub fn legendre_P2_e(x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_P2_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function evaluates the Legendre polynomials P_l(x) using explicit representations for l=1, 2, 3.
    pub fn legendre_P3_e(x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_P3_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function evaluates the Legendre polynomial P_l(x) for a specific value of l, x subject to l >= 0, |x| <= 1
    pub fn legendre_Pl(l: i32, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_Pl(l, x) }
    }

    /// This function evaluates the Legendre polynomial P_l(x) for a specific value of l, x subject to l >= 0, |x| <= 1
    pub fn legendre_Pl_e(l: i32, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_Pl_e(l, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function computes arrays of Legendre polynomials P_l(x) and derivatives dP_l(x)/dx, for l = 0, \dots, lmax, |x| <= 1
    pub fn legendre_Pl_array(x: f64, result_array: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_Pl_array(result_array.len() as i32, x, result_array.as_mut_ptr()) }
    }

    /// This function computes arrays of Legendre polynomials P_l(x) and derivatives dP_l(x)/dx, for l = 0, \dots, lmax, |x| <= 1
    pub fn legendre_Pl_deriv_array(x: f64, result_array: &mut [f64], result_deriv_array: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_Pl_deriv_array(result_array.len() as i32, x, result_array.as_mut_ptr(),
            result_deriv_array.as_mut_ptr()) }
    }

    /// This function computes the Legendre function Q_0(x) for x > -1, x != 1
    pub fn legendre_Q0(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_Q0(x) }
    }

    /// This function computes the Legendre function Q_0(x) for x > -1, x != 1
    pub fn legendre_Q0_e(x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_Q0_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function computes the Legendre function Q_0(x) for x > -1, x != 1.
    pub fn legendre_Q1(x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_Q1(x) }
    }

    /// This function computes the Legendre function Q_0(x) for x > -1, x != 1.
    pub fn legendre_Q1_e(x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_Q1_e(x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function computes the Legendre function Q_l(x) for x > -1, x != 1 and l >= 0.
    pub fn legendre_Ql(l: i32, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_Ql(l, x) }
    }

    /// This function computes the Legendre function Q_l(x) for x > -1, x != 1 and l >= 0.
    pub fn legendre_Ql_e(l: i32, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_Ql_e(l, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }
}

/// The following functions compute the associated Legendre Polynomials P_l^m(x).
/// Note that this function grows combinatorially with l and can overflow for l larger than about 150.
/// There is no trouble for small m, but overflow occurs when m and l are both large.
/// Rather than allow overflows, these functions refuse to calculate P_l^m(x) and return [`OvrFlw`](enums/type.Value.html) when they can sense that l and m are too big.
/// 
/// If you want to calculate a spherical harmonic, then do not use these functions. Instead use [`legendre_sphPlm`](fn.legendre_sphPlm.html) below, which uses a similar recursion, but with the normalized functions.
pub mod associated_polynomials {
    use ffi;
    use enums;
    use std::mem::zeroed;

    /// This routine computes the associated Legendre polynomial P_l^m(x) for m >= 0, l >= m, |x| <= 1.
    pub fn legendre_Plm(l: i32, m: i32, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_Plm(l, m, x) }
    }

    /// This routine computes the associated Legendre polynomial P_l^m(x) for m >= 0, l >= m, |x| <= 1.
    pub fn legendre_Plm_e(l: i32, m: i32, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_Plm_e(l, m, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function computes arrays of Legendre polynomials P_l^m(x) and derivatives dP_l^m(x)/dx, for m >= 0, l = |m|, ..., lmax, |x| <= 1.
    #[cfg(not(feature = "v2"))]
    pub fn legendre_Plm_array(m: i32, x: f64, result_array: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_Plm_array(result_array.len() as i32, m, x, result_array.as_mut_ptr()) }
    }

    /// This function computes arrays of Legendre polynomials P_l^m(x) and derivatives dP_l^m(x)/dx, for m >= 0, l = |m|, ..., lmax, |x| <= 1.
    #[cfg(not(feature = "v2"))]
    pub fn legendre_Plm_deriv_array(m: i32, x: f64, result_array: &mut [f64], result_deriv_array: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_Plm_deriv_array(result_array.len() as i32, m, x, result_array.as_mut_ptr(),
            result_deriv_array.as_mut_ptr()) }
    }

    /// This routine computes the normalized associated Legendre polynomial \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) suitable for use in spherical harmonics.
    /// The parameters must satisfy m >= 0, l >= m, |x| <= 1.
    /// This routine avoids the overflows that occur for the standard normalization of P_l^m(x).
    pub fn legendre_sphPlm(l: i32, m: i32, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_sphPlm(l, m, x) }
    }

    /// This routine computes the normalized associated Legendre polynomial \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x) suitable for use in spherical harmonics.
    /// The parameters must satisfy m >= 0, l >= m, |x| <= 1.
    /// This routine avoids the overflows that occur for the standard normalization of P_l^m(x).
    pub fn legendre_sphPlm_e(l: i32, m: i32, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_sphPlm_e(l, m, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function computes arrays of normalized associated Legendre functions \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x), and derivatives, for m >= 0, l = |m|, ..., lmax, |x| <= 1.0
    #[cfg(not(feature = "v2"))]
    pub fn legendre_sphPlm_array(m: i32, x: f64, result_array: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_sphPlm_array(result_array.len() as i32, m, x, result_array.as_mut_ptr()) }
    }

    /// This function computes arrays of normalized associated Legendre functions \sqrt{(2l+1)/(4\pi)} \sqrt{(l-m)!/(l+m)!} P_l^m(x), and derivatives, for m >= 0, l = |m|, ..., lmax, |x| <= 1.0
    #[cfg(not(feature = "v2"))]
    pub fn legendre_sphPlm_deriv_array(m: i32, x: f64, result_array: &mut [f64], result_deriv_array: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_sphPlm_deriv_array(result_array.len() as i32, m, x, result_array.as_mut_ptr(),
            result_deriv_array.as_mut_ptr()) }
    }

    /// This function returns the size of result_array[] needed for the array versions of P_l^m(x), lmax - m + 1.
    #[cfg(not(feature = "v2"))]
    pub fn legendre_array_size(lmax: i32, m: i32) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_array_size(lmax, m) }
    }
}

/// The Conical Functions P^\mu_{-(1/2)+i\lambda}(x) and Q^\mu_{-(1/2)+i\lambda} are described in Abramowitz & Stegun, Section 8.12.
pub mod conical {
    use ffi;
    use enums;
    use std::mem::zeroed;

    /// This routine computes the irregular Spherical Conical Function P^{1/2}_{-1/2 + i \lambda}(x) for x > -1.
    pub fn half(lambda: f64, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_conicalP_half(lambda, x) }
    }

    /// This routine computes the irregular Spherical Conical Function P^{1/2}_{-1/2 + i \lambda}(x) for x > -1.
    pub fn half_e(lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_conicalP_half_e(lambda, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the regular Spherical Conical Function P^{-1/2}_{-1/2 + i \lambda}(x) for x > -1.
    pub fn mhalf(lambda: f64, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_conicalP_mhalf(lambda, x) }
    }

    /// This routine computes the regular Spherical Conical Function P^{-1/2}_{-1/2 + i \lambda}(x) for x > -1.
    pub fn mhalf_e(lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_conicalP_mhalf_e(lambda, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the conical function P^0_{-1/2 + i \lambda}(x) for x > -1.
    pub fn _0(lambda: f64, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_conicalP_0(lambda, x) }
    }

    /// This routine computes the conical function P^0_{-1/2 + i \lambda}(x) for x > -1.
    pub fn _0_e(lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_conicalP_0_e(lambda, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the conical function P^1_{-1/2 + i \lambda}(x) for x > -1.
    pub fn _1(lambda: f64, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_conicalP_1(lambda, x) }
    }

    /// This routine computes the conical function P^1_{-1/2 + i \lambda}(x) for x > -1.
    pub fn _1_e(lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_conicalP_1_e(lambda, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the Regular Spherical Conical Function P^{-1/2-l}_{-1/2 + i \lambda}(x) for x > -1, l >= -1.
    pub fn sph_reg(l: i32, lambda: f64, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_conicalP_sph_reg(l, lambda, x) }
    }

    /// This routine computes the Regular Spherical Conical Function P^{-1/2-l}_{-1/2 + i \lambda}(x) for x > -1, l >= -1.
    pub fn sph_reg_e(l: i32, lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_conicalP_sph_reg_e(l, lambda, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the Regular Cylindrical Conical Function P^{-m}_{-1/2 + i \lambda}(x) for x > -1, m >= -1.
    pub fn cyl_reg(m: i32, lambda: f64, x: f64) -> f64 {
        unsafe { ffi::gsl_sf_conicalP_cyl_reg(m, lambda, x) }
    }

    /// This routine computes the Regular Cylindrical Conical Function P^{-m}_{-1/2 + i \lambda}(x) for x > -1, m >= -1.
    pub fn cyl_reg_e(m: i32, lambda: f64, x: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_conicalP_cyl_reg_e(m, lambda, x, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }
}

/// The following spherical functions are specializations of Legendre functions which give the regular eigenfunctions of the Laplacian on a 3-dimensional hyperbolic space H3d.
/// Of particular interest is the flat limit, \lambda \to \infty, \eta \to 0, \lambda\eta fixed.
pub mod radial {
    use ffi;
    use enums;
    use std::mem::zeroed;

    /// This routine computes the zeroth radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space, L^{H3d}_0(\lambda,\eta) := \sin(\lambda\eta)/(\lambda\sinh(\eta)) for \eta >= 0.
    /// In the flat limit this takes the form L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta).
    pub fn legendre_H3d_0(lambda: f64, eta: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_H3d_0(lambda, eta) }
    }

    /// This routine computes the zeroth radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space, L^{H3d}_0(\lambda,\eta) := \sin(\lambda\eta)/(\lambda\sinh(\eta)) for \eta >= 0.
    /// In the flat limit this takes the form L^{H3d}_0(\lambda,\eta) = j_0(\lambda\eta).
    pub fn legendre_H3d_0_e(lambda: f64, eta: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_H3d_0_e(lambda, eta, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the first radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space, L^{H3d}_1(\lambda,\eta) := 1/\sqrt{\lambda^2 + 1} \sin(\lambda \eta)/(\lambda \sinh(\eta))
    /// (\coth(\eta) - \lambda \cot(\lambda\eta)) for \eta >= 0.
    /// In the flat limit this takes the form L^{H3d}_1(\lambda,\eta) = j_1(\lambda\eta).
    pub fn legendre_H3d_1(lambda: f64, eta: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_H3d_1(lambda, eta) }
    }

    /// This routine computes the first radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space, L^{H3d}_1(\lambda,\eta) := 1/\sqrt{\lambda^2 + 1} \sin(\lambda \eta)/(\lambda \sinh(\eta))
    /// (\coth(\eta) - \lambda \cot(\lambda\eta)) for \eta >= 0.
    /// In the flat limit this takes the form L^{H3d}_1(\lambda,\eta) = j_1(\lambda\eta).
    pub fn legendre_H3d_1_e(lambda: f64, eta: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_H3d_1_e(lambda, eta, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This routine computes the l-th radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space \eta >= 0, l >= 0. In the flat limit this takes the form L^{H3d}_l(\lambda,\eta) = j_l(\lambda\eta).
    pub fn legendre_H3d(l: i32, lambda: f64, eta: f64) -> f64 {
        unsafe { ffi::gsl_sf_legendre_H3d(l, lambda, eta) }
    }

    /// This routine computes the l-th radial eigenfunction of the Laplacian on the 3-dimensional hyperbolic space \eta >= 0, l >= 0. In the flat limit this takes the form L^{H3d}_l(\lambda,\eta) = j_l(\lambda\eta).
    pub fn legendre_H3d_e(l: i32, lambda: f64, eta: f64) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_legendre_H3d_e(l, lambda, eta, &mut result) };

        (ret, ::types::Result{val: result.val, err: result.err})
    }

    /// This function computes an array of radial eigenfunctions L^{H3d}_l(\lambda, \eta) for 0 <= l <= lmax.
    pub fn legendre_H3d_array(lambda: f64, eta: f64, result_array: &mut [f64]) -> enums::Value {
        unsafe { ffi::gsl_sf_legendre_H3d_array(result_array.len() as i32, lambda, eta, result_array.as_mut_ptr()) }
    }
}