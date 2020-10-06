//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use enums;
use std::mem::MaybeUninit;

/// This routine computes the lowest-order normalized hydrogenic bound state radial wavefunction R_1 := 2Z \sqrt{Z} \exp(-Z r).
pub fn hydrogenicR_1(Z: f64, r: f64) -> f64 {
    unsafe { sys::gsl_sf_hydrogenicR_1(Z, r) }
}

/// This routine computes the lowest-order normalized hydrogenic bound state radial wavefunction R_1 := 2Z \sqrt{Z} \exp(-Z r).
pub fn hydrogenicR_1_e(Z: f64, r: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_hydrogenicR_1_e(Z, r, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This routine computes the n-th normalized hydrogenic bound state radial wavefunction,
///
/// R_n := 2 (Z^{3/2}/n^2) \sqrt{(n-l-1)!/(n+l)!} \exp(-Z r/n) (2Zr/n)^l
///           L^{2l+1}_{n-l-1}(2Zr/n).  
///
/// where L^a_b(x) is the generalized Laguerre polynomial (see [`Laguerre Functions`](http://www.gnu.org/software/gsl/manual/html_node/Laguerre-Functions.html#Laguerre-Functions)).
/// The normalization is chosen such that the wavefunction \psi is given by \psi(n,l,r) = R_n Y_{lm}.
pub fn hydrogenicR(n: i32, l: i32, Z: f64, r: f64) -> f64 {
    unsafe { sys::gsl_sf_hydrogenicR(n, l, Z, r) }
}

/// This routine computes the n-th normalized hydrogenic bound state radial wavefunction,
///
/// R_n := 2 (Z^{3/2}/n^2) \sqrt{(n-l-1)!/(n+l)!} \exp(-Z r/n) (2Zr/n)^l
///           L^{2l+1}_{n-l-1}(2Zr/n).  
///
/// where L^a_b(x) is the generalized Laguerre polynomial (see [`Laguerre Functions`](http://www.gnu.org/software/gsl/manual/html_node/Laguerre-Functions.html#Laguerre-Functions)).
/// The normalization is chosen such that the wavefunction \psi is given by \psi(n,l,r) = R_n Y_{lm}.
pub fn hydrogenicR_e(n: i32, l: i32, Z: f64, r: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_hydrogenicR_e(n, l, Z, r, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This function computes the Coulomb wave functions F_L(\eta,x), G_{L-k}(\eta,x) and their derivatives F'_L(\eta,x), G'_{L-k}(\eta,x) with respect to x. The parameters are restricted to L, L-k > -1/2, x > 0 and integer k. Note that L itself is not restricted to being an integer. The results are stored in the parameters F, G for the function values and Fp, Gp for the derivative values.
/// If an overflow occurs, GSL_EOVRFLW is returned and scaling exponents are stored in the modifiable parameters exp_F, exp_G.
pub fn wave_FG_e(
    eta: f64,
    x: f64,
    L_F: f64,
    k: i32,
    exp_F: &mut f64,
    exp_G: &mut f64,
) -> (
    enums::Value,
    ::types::Result,
    ::types::Result,
    ::types::Result,
    ::types::Result,
) {
    let mut F = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let mut Fp = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let mut G = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let mut Gp = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe {
        sys::gsl_sf_coulomb_wave_FG_e(
            eta,
            x,
            L_F,
            k,
            F.as_mut_ptr(),
            Fp.as_mut_ptr(),
            G.as_mut_ptr(),
            Gp.as_mut_ptr(),
            exp_F,
            exp_G,
        )
    };

    (
        enums::Value::from(ret),
        unsafe { F.assume_init() }.into(),
        unsafe { Fp.assume_init() }.into(),
        unsafe { G.assume_init() }.into(),
        unsafe { Gp.assume_init() }.into(),
    )
}

/// This function computes the Coulomb wave function F_L(\eta,x) for L = Lmin \dots Lmin + kmax, storing the results in fc_array.
/// In the case of overflow the exponent is stored in F_exponent.
pub fn wave_F_array(
    L_min: f64,
    eta: f64,
    x: f64,
    fc_array: &mut [f64],
    F_exponent: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        sys::gsl_sf_coulomb_wave_F_array(
            L_min,
            fc_array.len() as i32,
            eta,
            x,
            fc_array.as_mut_ptr(),
            F_exponent,
        )
    })
}

/// This function computes the functions F_L(\eta,x), G_L(\eta,x) for L = Lmin \dots Lmin + kmax storing the results in fc_array and gc_array.
/// In the case of overflow the exponents are stored in F_exponent and G_exponent.
pub fn wave_FG_array(
    L_min: f64,
    eta: f64,
    x: f64,
    fc_array: &mut [f64],
    gc_array: &mut [f64],
    F_exponent: &mut f64,
    G_exponent: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        sys::gsl_sf_coulomb_wave_FG_array(
            L_min,
            fc_array.len() as i32,
            eta,
            x,
            fc_array.as_mut_ptr(),
            gc_array.as_mut_ptr(),
            F_exponent,
            G_exponent,
        )
    })
}

/// This function computes the functions F_L(\eta,x), G_L(\eta,x) and their derivatives F'_L(\eta,x), G'_L(\eta,x) for L = Lmin \dots Lmin + kmax storing the results in fc_array, gc_array, fcp_array and gcp_array.
/// In the case of overflow the exponents are stored in F_exponent and G_exponent.
pub fn wave_FGp_array(
    L_min: f64,
    eta: f64,
    x: f64,
    fc_array: &mut [f64],
    fcp_array: &mut [f64],
    gc_array: &mut [f64],
    gcp_array: &mut [f64],
    F_exponent: &mut f64,
    G_exponent: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        sys::gsl_sf_coulomb_wave_FGp_array(
            L_min,
            fc_array.len() as i32,
            eta,
            x,
            fc_array.as_mut_ptr(),
            fcp_array.as_mut_ptr(),
            gc_array.as_mut_ptr(),
            gcp_array.as_mut_ptr(),
            F_exponent,
            G_exponent,
        )
    })
}

/// This function computes the Coulomb wave function divided by the argument F_L(\eta, x)/x for L = Lmin \dots Lmin + kmax, storing the results in fc_array.
/// In the case of overflow the exponent is stored in F_exponent. This function reduces to spherical Bessel functions in the limit \eta \to 0.
pub fn wave_sphF_array(
    L_min: f64,
    eta: f64,
    x: f64,
    fc_array: &mut [f64],
    F_exponent: &mut f64,
) -> enums::Value {
    enums::Value::from(unsafe {
        sys::gsl_sf_coulomb_wave_sphF_array(
            L_min,
            fc_array.len() as i32,
            eta,
            x,
            fc_array.as_mut_ptr(),
            F_exponent,
        )
    })
}

/// This function computes the Coulomb wave function normalization constant C_L(\eta) for L > -1.
pub fn CL_e(L: f64, eta: f64) -> (enums::Value, ::types::Result) {
    let mut result = unsafe { MaybeUninit::<sys::gsl_sf_result>::uninit() };
    let ret = unsafe { sys::gsl_sf_coulomb_CL_e(L, eta, result.as_mut_ptr()) };

    (::Value::from(ret), unsafe { result.assume_init() }.into())
}

/// This function computes the Coulomb wave function normalization constant C_L(\eta) for L = Lmin \dots Lmin + kmax, Lmin > -1.
pub fn CL_array(Lmin: f64, eta: f64, cl: &mut [f64]) -> enums::Value {
    enums::Value::from(unsafe {
        sys::gsl_sf_coulomb_CL_array(Lmin, cl.len() as i32, eta, cl.as_mut_ptr())
    })
}
