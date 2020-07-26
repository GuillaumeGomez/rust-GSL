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
        use enums;
        use ffi;
        use std::mem::zeroed;

        /// This routine computes the complete elliptic integral K(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_Kcomp(k: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_Kcomp(k, mode) }
        }

        /// This routine computes the complete elliptic integral K(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_Kcomp_e(k: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
            let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
            let ret = unsafe { ::ffi::gsl_sf_ellint_Kcomp_e(k, mode, &mut result) };

            (
                enums::Value::from(ret),
                ::types::Result {
                    val: result.val,
                    err: result.err,
                },
            )
        }

        /// This routine computes the complete elliptic integral E(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_Ecomp(k: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_Ecomp(k, mode) }
        }

        /// This routine computes the complete elliptic integral E(k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_Ecomp_e(k: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
            let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
            let ret = unsafe { ::ffi::gsl_sf_ellint_Ecomp_e(k, mode, &mut result) };

            (
                enums::Value::from(ret),
                ::types::Result {
                    val: result.val,
                    err: result.err,
                },
            )
        }

        /// This routine computes the complete elliptic integral \Pi(k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        pub fn ellint_Pcomp(k: f64, n: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_Pcomp(k, n, mode) }
        }

        /// This routine computes the complete elliptic integral \Pi(k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        pub fn ellint_Pcomp_e(k: f64, n: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
            let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
            let ret = unsafe { ::ffi::gsl_sf_ellint_Pcomp_e(k, n, mode, &mut result) };

            (
                enums::Value::from(ret),
                ::types::Result {
                    val: result.val,
                    err: result.err,
                },
            )
        }
    }

    pub mod incomplete {
        use enums;
        use ffi;
        use std::mem::zeroed;

        /// This routine computes the incomplete elliptic integral F(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_F(phi: f64, k: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_F(phi, k, mode) }
        }

        /// This routine computes the incomplete elliptic integral F(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_F_e(phi: f64, k: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
            let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
            let ret = unsafe { ::ffi::gsl_sf_ellint_F_e(phi, k, mode, &mut result) };

            (
                enums::Value::from(ret),
                ::types::Result {
                    val: result.val,
                    err: result.err,
                },
            )
        }

        /// This routine computes the incomplete elliptic integral E(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_E(phi: f64, k: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_E(phi, k, mode) }
        }

        /// This routine computes the incomplete elliptic integral E(\phi,k) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameter m = k^2.
        pub fn ellint_E_e(phi: f64, k: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
            let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
            let ret = unsafe { ::ffi::gsl_sf_ellint_E_e(phi, k, mode, &mut result) };

            (
                enums::Value::from(ret),
                ::types::Result {
                    val: result.val,
                    err: result.err,
                },
            )
        }

        /// This routine computes the incomplete elliptic integral \Pi(\phi,k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        pub fn ellint_P(phi: f64, k: f64, n: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_P(phi, k, n, mode) }
        }

        /// This routine computes the incomplete elliptic integral \Pi(\phi,k,n) to the accuracy specified by the mode variable mode.
        /// Note that Abramowitz & Stegun define this function in terms of the parameters m = k^2 and \sin^2(\alpha) = k^2, with the change of sign n \to -n.
        pub fn ellint_P_e(
            phi: f64,
            k: f64,
            n: f64,
            mode: ::Mode,
        ) -> (enums::Value, ::types::Result) {
            let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
            let ret = unsafe { ::ffi::gsl_sf_ellint_P_e(phi, k, n, mode, &mut result) };

            (
                enums::Value::from(ret),
                ::types::Result {
                    val: result.val,
                    err: result.err,
                },
            )
        }

        /// This routine computes the incomplete elliptic integral D(\phi,k) which is defined through the Carlson form RD(x,y,z) by the following relation,
        ///
        /// D(\phi,k,n) = (1/3)(\sin(\phi))^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1).
        #[cfg(feature = "v2")]
        pub fn ellint_D(phi: f64, k: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_D(phi, k, mode) }
        }

        /// This routine computes the incomplete elliptic integral D(\phi,k) which is defined through the Carlson form RD(x,y,z) by the following relation,
        ///
        /// D(\phi,k,n) = (1/3)(\sin(\phi))^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1).
        ///
        /// The argument n is not used and will be removed in a future release.
        #[cfg(not(feature = "v2"))]
        pub fn ellint_D(phi: f64, k: f64, n: f64, mode: ::Mode) -> f64 {
            unsafe { ffi::gsl_sf_ellint_D(phi, k, n, mode) }
        }

        /// This routine computes the incomplete elliptic integral D(\phi,k) which is defined through the Carlson form RD(x,y,z) by the following relation,
        ///
        /// D(\phi,k,n) = (1/3)(\sin(\phi))^3 RD (1-\sin^2(\phi), 1-k^2 \sin^2(\phi), 1).
        ///
        /// The argument n is not used and will be removed in a future release.
        pub fn ellint_D_e(
            phi: f64,
            k: f64,
            n: f64,
            mode: ::Mode,
        ) -> (enums::Value, ::types::Result) {
            let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
            let ret = unsafe { ::ffi::gsl_sf_ellint_D_e(phi, k, n, mode, &mut result) };

            (
                enums::Value::from(ret),
                ::types::Result {
                    val: result.val,
                    err: result.err,
                },
            )
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
    use enums;
    use ffi;
    use std::mem::zeroed;

    /// This routine computes the incomplete elliptic integral RC(x,y) to the accuracy specified by the mode variable mode.
    pub fn ellint_RC(x: f64, y: f64, mode: ::Mode) -> f64 {
        unsafe { ffi::gsl_sf_ellint_RC(x, y, mode) }
    }

    /// This routine computes the incomplete elliptic integral RC(x,y) to the accuracy specified by the mode variable mode.
    pub fn ellint_RC_e(x: f64, y: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_ellint_RC_e(x, y, mode, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }

    /// This routine computes the incomplete elliptic integral RD(x,y,z) to the accuracy specified by the mode variable mode.
    pub fn ellint_RD(x: f64, y: f64, z: f64, mode: ::Mode) -> f64 {
        unsafe { ffi::gsl_sf_ellint_RD(x, y, z, mode) }
    }

    /// This routine computes the incomplete elliptic integral RD(x,y,z) to the accuracy specified by the mode variable mode.
    pub fn ellint_RD_e(x: f64, y: f64, z: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_ellint_RD_e(x, y, z, mode, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }

    /// This routine computes the incomplete elliptic integral RF(x,y,z) to the accuracy specified by the mode variable mode.
    pub fn ellint_RF(x: f64, y: f64, z: f64, mode: ::Mode) -> f64 {
        unsafe { ffi::gsl_sf_ellint_RF(x, y, z, mode) }
    }

    /// This routine computes the incomplete elliptic integral RF(x,y,z) to the accuracy specified by the mode variable mode.
    pub fn ellint_RF_e(x: f64, y: f64, z: f64, mode: ::Mode) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_ellint_RF_e(x, y, z, mode, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }

    /// This routine computes the incomplete elliptic integral RJ(x,y,z,p) to the accuracy specified by the mode variable mode.
    pub fn ellint_RJ(x: f64, y: f64, z: f64, p: f64, mode: ::Mode) -> f64 {
        unsafe { ffi::gsl_sf_ellint_RJ(x, y, z, p, mode) }
    }

    /// This routine computes the incomplete elliptic integral RJ(x,y,z,p) to the accuracy specified by the mode variable mode.
    pub fn ellint_RJ_e(
        x: f64,
        y: f64,
        z: f64,
        p: f64,
        mode: ::Mode,
    ) -> (enums::Value, ::types::Result) {
        let mut result = unsafe { zeroed::<::ffi::gsl_sf_result>() };
        let ret = unsafe { ::ffi::gsl_sf_ellint_RJ_e(x, y, z, p, mode, &mut result) };

        (
            enums::Value::from(ret),
            ::types::Result {
                val: result.val,
                err: result.err,
            },
        )
    }
}
