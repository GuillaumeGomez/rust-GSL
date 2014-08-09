//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod Coulomb {
    use ffi;
    use Gsl;
    use std::mem::zeroed;
    use enums;

    /// This routine computes the lowest-order normalized hydrogenic bound state radial wavefunction R_1 := 2Z \sqrt{Z} \exp(-Z r).
    pub fn hydrogenicR_1(Z: f64, r: f64) -> f64 {
        unsafe { ffi::gsl_sf_hydrogenicR_1(Z, r) }
    }

    /// This routine computes the lowest-order normalized hydrogenic bound state radial wavefunction R_1 := 2Z \sqrt{Z} \exp(-Z r).
    pub fn hydrogenicR_1_e(Z: f64, r: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_hydrogenicR_1_e(Z, r, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This routine computes the n-th normalized hydrogenic bound state radial wavefunction,
    /// 
    /// R_n := 2 (Z^{3/2}/n^2) \sqrt{(n-l-1)!/(n+l)!} \exp(-Z r/n) (2Zr/n)^l
    ///           L^{2l+1}_{n-l-1}(2Zr/n).  
    /// 
    /// where L^a_b(x) is the generalized Laguerre polynomial (see [`Laguerre Functions`](http://www.gnu.org/software/gsl/manual/html_node/Laguerre-Functions.html#Laguerre-Functions)).
    /// The normalization is chosen such that the wavefunction \psi is given by \psi(n,l,r) = R_n Y_{lm}.
    pub fn hydrogenicR(n: i32, l: i32, Z: f64, r: f64) -> f64 {
        unsafe { ffi::gsl_sf_hydrogenicR(n, l, Z, r) }
    }

    /// This routine computes the n-th normalized hydrogenic bound state radial wavefunction,
    /// 
    /// R_n := 2 (Z^{3/2}/n^2) \sqrt{(n-l-1)!/(n+l)!} \exp(-Z r/n) (2Zr/n)^l
    ///           L^{2l+1}_{n-l-1}(2Zr/n).  
    /// 
    /// where L^a_b(x) is the generalized Laguerre polynomial (see [`Laguerre Functions`](http://www.gnu.org/software/gsl/manual/html_node/Laguerre-Functions.html#Laguerre-Functions)).
    /// The normalization is chosen such that the wavefunction \psi is given by \psi(n,l,r) = R_n Y_{lm}.
    pub fn hydrogenicR_e(n: i32, l: i32, Z: f64, r: f64) -> (enums::GslValue, Gsl::Result) {
        let mut result = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_hydrogenicR_e(n, l, Z, r, &mut result) };

        (ret, Gsl::Result{val: result.val, err: result.err})
    }

    /// This function computes the Coulomb wave functions F_L(\eta,x), G_{L-k}(\eta,x) and their derivatives F'_L(\eta,x), G'_{L-k}(\eta,x) with respect to x. The parameters are restricted to L, L-k > -1/2, x > 0 and integer k. Note that L itself is not restricted to being an integer. The results are stored in the parameters F, G for the function values and Fp, Gp for the derivative values.
    /// If an overflow occurs, GSL_EOVRFLW is returned and scaling exponents are stored in the modifiable parameters exp_F, exp_G.
    pub fn wave_FG_e(eta: f64, x: f64, L_F: f64, k: i32, exp_F: &mut f64, exp_G: &mut f64) -> (enums::GslValue, Gsl::Result, Gsl::Result, Gsl::Result, Gsl::Result) {
        let mut F = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let mut Fp = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let mut G = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let mut Gp = unsafe { zeroed::<ffi::gsl_sf_result>() };
        let ret = unsafe { ffi::gsl_sf_coulomb_wave_FG_e(eta, x, L_F, k, &mut F, &mut Fp, &mut G, &mut Gp, exp_f, exp_G); }

        (ret, F, Fp, G, Gp)
    }
}