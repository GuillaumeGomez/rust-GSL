//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// The complete Fermi-Dirac integral F_j(x) is given by,
///
/// F_j(x)   := (1/\Gamma(j+1)) \int_0^\infty dt (t^j / (\exp(t-x) + 1))
///
/// Note that the Fermi-Dirac integral is sometimes defined without the normalisation factor in other texts.
pub mod complete_integrals {
    use crate::Value;
    use std::mem::MaybeUninit;

    /// This routine computes the complete Fermi-Dirac integral with an index of -1.
    /// This integral is given by F_{-1}(x) = e^x / (1 + e^x).
    pub fn fermi_dirac_m1(x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_m1(x) }
    }

    /// This routine computes the complete Fermi-Dirac integral with an index of -1.
    /// This integral is given by F_{-1}(x) = e^x / (1 + e^x).
    pub fn fermi_dirac_m1_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_m1_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complete Fermi-Dirac integral with an index of 0.
    /// This integral is given by F_0(x) = \ln(1 + e^x).
    pub fn fermi_dirac_0(x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_0(x) }
    }

    /// This routine computes the complete Fermi-Dirac integral with an index of 0.
    /// This integral is given by F_0(x) = \ln(1 + e^x).
    pub fn fermi_dirac_0_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_0_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complete Fermi-Dirac integral with an index of 1, F_1(x) = \int_0^\infty dt (t /(\exp(t-x)+1)).
    pub fn fermi_dirac_1(x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_1(x) }
    }

    /// This routine computes the complete Fermi-Dirac integral with an index of 1, F_1(x) = \int_0^\infty dt (t /(\exp(t-x)+1)).
    pub fn fermi_dirac_1_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_1_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complete Fermi-Dirac integral with an index of 2, F_2(x) = (1/2) \int_0^\infty dt (t^2 /(\exp(t-x)+1)).
    pub fn fermi_dirac_2(x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_2(x) }
    }

    /// This routine computes the complete Fermi-Dirac integral with an index of 2, F_2(x) = (1/2) \int_0^\infty dt (t^2 /(\exp(t-x)+1)).
    pub fn fermi_dirac_2_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_2_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complete Fermi-Dirac integral with an integer index of j, F_j(x) = (1/\Gamma(j+1)) \int_0^\infty dt (t^j /(\exp(t-x)+1)).
    pub fn fermi_dirac_int(j: i32, x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_int(j, x) }
    }

    /// This routine computes the complete Fermi-Dirac integral with an integer index of j, F_j(x) = (1/\Gamma(j+1)) \int_0^\infty dt (t^j /(\exp(t-x)+1)).
    pub fn fermi_dirac_int_e(j: i32, x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_int_e(j, x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complete Fermi-Dirac integral F_{-1/2}(x).
    pub fn fermi_dirac_mhalf(x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_mhalf(x) }
    }

    /// This routine computes the complete Fermi-Dirac integral F_{-1/2}(x).
    pub fn fermi_dirac_mhalf_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_mhalf_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complete Fermi-Dirac integral F_{1/2}(x).
    pub fn fermi_dirac_half(x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_half(x) }
    }

    /// This routine computes the complete Fermi-Dirac integral F_{1/2}(x).
    pub fn fermi_dirac_half_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_half_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complete Fermi-Dirac integral F_{3/2}(x).
    pub fn fermi_dirac_3half(x: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_3half(x) }
    }

    /// This routine computes the complete Fermi-Dirac integral F_{3/2}(x).
    pub fn fermi_dirac_3half_e(x: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_3half_e(x, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }
}

/// The incomplete Fermi-Dirac integral F_j(x,b) is given by,
///
/// F_j(x,b)   := (1/\Gamma(j+1)) \int_b^\infty dt (t^j / (\Exp(t-x) + 1))
pub mod incomplete_integrals {
    use crate::Value;
    use std::mem::MaybeUninit;

    /// This routine computes the incomplete Fermi-Dirac integral with an index of zero, F_0(x,b) = \ln(1 + e^{b-x}) - (b-x).
    pub fn fermi_dirac_inc_0(x: f64, b: f64) -> f64 {
        unsafe { sys::gsl_sf_fermi_dirac_inc_0(x, b) }
    }

    /// This routine computes the incomplete Fermi-Dirac integral with an index of zero, F_0(x,b) = \ln(1 + e^{b-x}) - (b-x).
    pub fn fermi_dirac_inc_0_e(x: f64, b: f64) -> (Value, ::types::Result) {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { ::sys::gsl_sf_fermi_dirac_inc_0_e(x, b, result.as_mut_ptr()) };

        (::Value::from(ret), unsafe { result.assume_init() }.into())
    }
}
