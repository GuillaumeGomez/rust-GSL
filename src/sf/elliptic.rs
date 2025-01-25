//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Further information about the elliptic integrals can be found in Abramowitz & Stegun, Chapter 17.

/// The Legendre forms of elliptic integrals F(\phi,k), E(\phi,k) and \Pi(\phi,k,n) are defined by,
///
/// F(\phi,k) = \int_0^\phi dt 1/\sqrt((1 - k^2 \sin^2(t)))
///
/// E(\phi,k) = \int_0^\phi dt   \sqrt((1 - k^2 \sin^2(t)))
///
/// Pi(\phi,k,n) = \int_0^\phi dt 1/((1 + n \sin^2(t))\sqrt(1 - k^2 \sin^2(t)))
///
/// The complete Legendre forms are denoted by K(k) = F(\pi/2, k) and E(k) = E(\pi/2, k).
///
/// The notation used here is based on Carlson, Numerische Mathematik 33 (1979) 1 and differs slightly from that used by Abramowitz & Stegun, where the functions are given in terms of the parameter m = k^2 and n is replaced by -n.
pub mod legendre {
    pub mod complete {
        use crate::{types, Error};
        use std::mem::MaybeUninit;

        /// This routine computes the complete elliptic integral K(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_Kcomp")]
        pub fn ellint_Kcomp(k: f64, mode: crate::Mode) -> f64 {
            unsafe { sys::gsl_sf_ellint_Kcomp(k, mode.into()) }
        }

        /// This routine computes the complete elliptic integral K(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_Kcomp_e")]
        pub fn ellint_Kcomp_e(k: f64, mode: crate::Mode) -> Result<types::Result, Error> {
            let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
            let ret = unsafe { sys::gsl_sf_ellint_Kcomp_e(k, mode.into(), result.as_mut_ptr()) };

            Error::handle(ret, unsafe { result.assume_init() }.into())
        }

        /// This routine computes the complete elliptic integral E(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_Ecomp")]
        pub fn ellint_Ecomp(k: f64, mode: crate::Mode) -> f64 {
            unsafe { sys::gsl_sf_ellint_Ecomp(k, mode.into()) }
        }

        /// This routine computes the complete elliptic integral E(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_Ecomp_e")]
        pub fn ellint_Ecomp_e(k: f64, mode: crate::Mode) -> Result<types::Result, Error> {
            let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
            let ret = unsafe { sys::gsl_sf_ellint_Ecomp_e(k, mode.into(), result.as_mut_ptr()) };

            Error::handle(ret, unsafe { result.assume_init() }.into())
        }

        /// This routine computes the complete elliptic integral \Pi(k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        #[doc(alias = "gsl_sf_ellint_Pcomp")]
        pub fn ellint_Pcomp(k: f64, n: f64, mode: crate::Mode) -> f64 {
            unsafe { sys::gsl_sf_ellint_Pcomp(k, n, mode.into()) }
        }

        /// This routine computes the complete elliptic integral \Pi(k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        #[doc(alias = "gsl_sf_ellint_Pcomp_e")]
        pub fn ellint_Pcomp_e(k: f64, n: f64, mode: crate::Mode) -> Result<types::Result, Error> {
            let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
            let ret = unsafe { sys::gsl_sf_ellint_Pcomp_e(k, n, mode.into(), result.as_mut_ptr()) };

            Error::handle(ret, unsafe { result.assume_init() }.into())
        }
    }

    pub mod incomplete {
        use crate::{types, Error};
        use std::mem::MaybeUninit;

        /// This routine computes the incomplete elliptic integral F(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_F")]
        pub fn ellint_F(phi: f64, k: f64, mode: crate::Mode) -> f64 {
            unsafe { sys::gsl_sf_ellint_F(phi, k, mode.into()) }
        }

        /// This routine computes the incomplete elliptic integral F(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_F_e")]
        pub fn ellint_F_e(phi: f64, k: f64, mode: crate::Mode) -> Result<types::Result, Error> {
            let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
            let ret = unsafe { sys::gsl_sf_ellint_F_e(phi, k, mode.into(), result.as_mut_ptr()) };

            Error::handle(ret, unsafe { result.assume_init() }.into())
        }

        /// This routine computes the incomplete elliptic integral E(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_E")]
        pub fn ellint_E(phi: f64, k: f64, mode: crate::Mode) -> f64 {
            unsafe { sys::gsl_sf_ellint_E(phi, k, mode.into()) }
        }

        /// This routine computes the incomplete elliptic integral E(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        #[doc(alias = "gsl_sf_ellint_E_e")]
        pub fn ellint_E_e(phi: f64, k: f64, mode: crate::Mode) -> Result<types::Result, Error> {
            let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
            let ret = unsafe { sys::gsl_sf_ellint_E_e(phi, k, mode.into(), result.as_mut_ptr()) };

            Error::handle(ret, unsafe { result.assume_init() }.into())
        }

        /// This routine computes the incomplete elliptic integral \Pi(\phi,k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        #[doc(alias = "gsl_sf_ellint_P")]
        pub fn ellint_P(phi: f64, k: f64, n: f64, mode: crate::Mode) -> f64 {
            unsafe { sys::gsl_sf_ellint_P(phi, k, n, mode.into()) }
        }

        /// This routine computes the incomplete elliptic integral \Pi(\phi,k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        #[doc(alias = "gsl_sf_ellint_P_e")]
        pub fn ellint_P_e(
            phi: f64,
            k: f64,
            n: f64,
            mode: crate::Mode,
        ) -> Result<types::Result, Error> {
            let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
            let ret =
                unsafe { sys::gsl_sf_ellint_P_e(phi, k, n, mode.into(), result.as_mut_ptr()) };

            Error::handle(ret, unsafe { result.assume_init() }.into())
        }

        /// This routine computes the incomplete elliptic integral D(\phi,k) which is defined through the Carlson form RD(x,y,z) by the following relation,
        ///
        /// D(\phi,k,n) = (1/3)(\sin(\phi))^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1).
        #[doc(alias = "gsl_sf_ellint_D")]
        pub fn ellint_D(phi: f64, k: f64, mode: crate::Mode) -> f64 {
            unsafe { sys::gsl_sf_ellint_D(phi, k, mode.into()) }
        }

        /// This routine computes the incomplete elliptic integral D(\phi,k) which is defined through the Carlson form RD(x,y,z) by the following relation,
        ///
        /// D(\phi,k,n) = (1/3)(\sin(\phi))^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1).
        ///
        /// The argument n is not used and will be removed in a future release.
        #[doc(alias = "gsl_sf_ellint_D_e")]
        pub fn ellint_D_e(phi: f64, k: f64, mode: crate::Mode) -> Result<types::Result, Error> {
            let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
            let ret = unsafe { sys::gsl_sf_ellint_D_e(phi, k, mode.into(), result.as_mut_ptr()) };

            Error::handle(ret, unsafe { result.assume_init() }.into())
        }
    }
}

/// The Carlson symmetric forms of elliptical integrals RC(x,y), RD(x,y,z), RF(x,y,z) and RJ(x,y,z,p) are defined by,
///
/// RC(x,y) = 1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1)
///
/// RD(x,y,z) = 3/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2)
///
/// RF(x,y,z) = 1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2)
///
/// RJ(x,y,z,p) = 3/2 \int_0^\infty dt
///                (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1)
pub mod carlson {
    use crate::{types, Error};
    use std::mem::MaybeUninit;

    /// This routine computes the incomplete elliptic integral RC(x,y) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RC")]
    pub fn ellint_RC(x: f64, y: f64, mode: crate::Mode) -> f64 {
        unsafe { sys::gsl_sf_ellint_RC(x, y, mode.into()) }
    }

    /// This routine computes the incomplete elliptic integral RC(x,y) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RC_e")]
    pub fn ellint_RC_e(x: f64, y: f64, mode: crate::Mode) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_ellint_RC_e(x, y, mode.into(), result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the incomplete elliptic integral RD(x,y,z) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RD")]
    pub fn ellint_RD(x: f64, y: f64, z: f64, mode: crate::Mode) -> f64 {
        unsafe { sys::gsl_sf_ellint_RD(x, y, z, mode.into()) }
    }

    /// This routine computes the incomplete elliptic integral RD(x,y,z) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RD_e")]
    pub fn ellint_RD_e(x: f64, y: f64, z: f64, mode: crate::Mode) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_ellint_RD_e(x, y, z, mode.into(), result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the incomplete elliptic integral RF(x,y,z) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RF")]
    pub fn ellint_RF(x: f64, y: f64, z: f64, mode: crate::Mode) -> f64 {
        unsafe { sys::gsl_sf_ellint_RF(x, y, z, mode.into()) }
    }

    /// This routine computes the incomplete elliptic integral RF(x,y,z) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RF_e")]
    pub fn ellint_RF_e(x: f64, y: f64, z: f64, mode: crate::Mode) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_ellint_RF_e(x, y, z, mode.into(), result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the incomplete elliptic integral RJ(x,y,z,p) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RJ")]
    pub fn ellint_RJ(x: f64, y: f64, z: f64, p: f64, mode: crate::Mode) -> f64 {
        unsafe { sys::gsl_sf_ellint_RJ(x, y, z, p, mode.into()) }
    }

    /// This routine computes the incomplete elliptic integral RJ(x,y,z,p) to the accuracy specified by the mode variable mode.
    #[doc(alias = "gsl_sf_ellint_RJ_e")]
    pub fn ellint_RJ_e(
        x: f64,
        y: f64,
        z: f64,
        p: f64,
        mode: crate::Mode,
    ) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_ellint_RJ_e(x, y, z, p, mode.into(), result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }
}
