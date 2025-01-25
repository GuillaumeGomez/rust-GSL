//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! This following routines compute the gamma and beta functions in their full and incomplete forms, as well as various kinds of factorials.

/// The Gamma function is defined by the following integral,
///
/// \Gamma(x) = \int_0^\infty dt  t^{x-1} \exp(-t)
///
/// It is related to the factorial function by \Gamma(n)=(n-1)! for positive integer n.
/// Further information on the Gamma function can be found in Abramowitz & Stegun, Chapter 6.
pub mod gamma {
    use crate::{types, Error};
    use std::mem::MaybeUninit;

    /// These routines compute the Gamma function \Gamma(x), subject to x not being a negative integer or zero. The function is computed using the real Lanczos method.
    /// The maximum value of x such that \Gamma(x) is not considered an overflow is given by the macro GSL_SF_GAMMA_XMAX and is 171.0.
    #[doc(alias = "gsl_sf_gamma")]
    pub fn gamma(x: f64) -> f64 {
        unsafe { sys::gsl_sf_gamma(x) }
    }

    /// This routine provides an exponential function \exp(x) using GSL semantics and error checking.
    #[doc(alias = "gsl_sf_gamma_e")]
    pub fn gamma_e(x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_gamma_e(x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the Gamma function \Gamma(x), subject to x not being a negative integer or zero.
    /// The function is computed using the real Lanczos method. The maximum value of x such that \Gamma(x) is not considered an overflow is given by the macro GSL_SF_GAMMA_XMAX and is 171.0.
    #[doc(alias = "gsl_sf_lngamma")]
    pub fn lngamma(x: f64) -> f64 {
        unsafe { sys::gsl_sf_lngamma(x) }
    }

    /// This routine computes the Gamma function \Gamma(x), subject to x not being a negative integer or zero.
    /// The function is computed using the real Lanczos method. The maximum value of x such that \Gamma(x) is not considered an overflow is given by the macro GSL_SF_GAMMA_XMAX and is 171.0.
    #[doc(alias = "gsl_sf_lngamma_e")]
    pub fn lngamma_e(x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lngamma_e(x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the sign of the gamma function and the logarithm of its magnitude, subject to x not being a negative integer or zero.
    /// The function is computed using the real Lanczos method.
    /// The value of the gamma function and its error can be reconstructed using the relation \Gamma(x) = sgn * \exp(result\_lg), taking into account the two components of result_lg.
    #[doc(alias = "gsl_sf_lngamma_sgn_e")]
    pub fn lngamma_sgn_e(x: f64, sgn: &mut f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lngamma_sgn_e(x, result.as_mut_ptr(), sgn) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the regulated Gamma Function \Gamma^*(x) for x > 0. The regulated gamma function is given by,
    ///
    /// ```latex
    /// \Gamma^*(x) = \Gamma(x)/(\sqrt{2\pi} x^{(x-1/2)} \exp(-x))
    ///
    ///             = (1 + (1/12x) + ...)  for x to infty
    /// ```
    ///
    /// and is a useful suggestion of Temme.
    #[doc(alias = "gsl_sf_gammastar")]
    pub fn gammastar(x: f64) -> f64 {
        unsafe { sys::gsl_sf_gammastar(x) }
    }

    /// This routine computes the regulated Gamma Function \Gamma^*(x) for x > 0. The regulated gamma function is given by,
    ///
    /// ```latex
    /// \Gamma^*(x) = \Gamma(x)/(\sqrt{2\pi} x^{(x-1/2)} \exp(-x))
    ///
    ///             = (1 + (1/12x) + ...)  for x to infty
    /// ```
    ///
    /// and is a useful suggestion of Temme.
    #[doc(alias = "gsl_sf_gammastar_e")]
    pub fn gammastar_e(x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_gammastar_e(x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the reciprocal of the gamma function, 1/\Gamma(x) using the real Lanczos method.
    #[doc(alias = "gsl_sf_gammainv")]
    pub fn gammainv(x: f64) -> f64 {
        unsafe { sys::gsl_sf_gammainv(x) }
    }

    /// This routine computes the reciprocal of the gamma function, 1/\Gamma(x) using the real Lanczos method.
    #[doc(alias = "gsl_sf_gammainv_e")]
    pub fn gammainv_e(x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_gammainv_e(x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes \log(\Gamma(z)) for complex z=z_r+i z_i and z not a negative integer or zero, using the complex Lanczos method.
    /// The returned parameters are lnr = \log|\Gamma(z)| and arg = \arg(\Gamma(z)) in (-\pi,\pi]. Note that the phase part (arg) is not well-determined when |z| is very large, due to inevitable roundoff in restricting to (-\pi,\pi].
    /// This will result in a GSL_ELOSS error when it occurs. The absolute value part (lnr), however, never suffers from loss of precision.
    #[doc(alias = "gsl_sf_lngamma_complex_e")]
    pub fn lngamma_complex_e(zr: f64, zi: f64) -> Result<(types::Result, types::Result), Error> {
        let mut lnr = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let mut arg = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret =
            unsafe { sys::gsl_sf_lngamma_complex_e(zr, zi, lnr.as_mut_ptr(), arg.as_mut_ptr()) };

        Error::handle(
            ret,
            (
                unsafe { lnr.assume_init() }.into(),
                unsafe { arg.assume_init() }.into(),
            ),
        )
    }
}

/// Although factorials can be computed from the Gamma function, using the relation n! = \Gamma(n+1) for non-negative integer n, it is usually more
/// efficient to call the functions in this section, particularly for small values of n, whose factorial values are maintained in hardcoded tables.
pub mod factorials {
    use crate::{types, Error};
    use std::mem::MaybeUninit;

    /// This routine computes the factorial n!. The factorial is related to the Gamma function by n! = \Gamma(n+1).
    /// The maximum value of n such that n! is not considered an overflow is given by the macro SF_FACT_NMAX and is 170.
    #[doc(alias = "gsl_sf_fact")]
    pub fn fact(n: u32) -> f64 {
        unsafe { sys::gsl_sf_fact(n) }
    }

    /// This routine computes the factorial n!. The factorial is related to the Gamma function by n! = \Gamma(n+1).
    /// The maximum value of n such that n! is not considered an overflow is given by the macro SF_FACT_NMAX and is 170.
    #[doc(alias = "gsl_sf_fact_e")]
    pub fn fact_e(n: u32) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_fact_e(n, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the double factorial n!! = n(n-2)(n-4) \dots.
    /// The maximum value of n such that n!! is not considered an overflow is given by the macro SF_DOUBLEFACT_NMAX and is 297.
    #[doc(alias = "gsl_sf_doublefact")]
    pub fn doublefact(n: u32) -> f64 {
        unsafe { sys::gsl_sf_doublefact(n) }
    }

    /// This routine computes the double factorial n!! = n(n-2)(n-4) \dots.
    /// The maximum value of n such that n!! is not considered an overflow is given by the macro SF_DOUBLEFACT_NMAX and is 297.
    #[doc(alias = "gsl_sf_doublefact_e")]
    pub fn doublefact_e(n: u32) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_doublefact_e(n, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the logarithm of the factorial of n, \log(n!).
    /// The algorithm is faster than computing \ln(\Gamma(n+1)) via gsl_sf_lngamma for n < 170, but defers for larger n.
    #[doc(alias = "gsl_sf_lnfact")]
    pub fn lnfact(n: u32) -> f64 {
        unsafe { sys::gsl_sf_lnfact(n) }
    }

    /// This routine computes the logarithm of the factorial of n, \log(n!).
    /// The algorithm is faster than computing \ln(\Gamma(n+1)) via gsl_sf_lngamma for n < 170, but defers for larger n.
    #[doc(alias = "gsl_sf_lnfact_e")]
    pub fn lnfact_e(n: u32) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lnfact_e(n, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the logarithm of the double factorial of n, \log(n!!).
    #[doc(alias = "gsl_sf_lndoublefact")]
    pub fn lndoublefact(n: u32) -> f64 {
        unsafe { sys::gsl_sf_lndoublefact(n) }
    }

    /// This routine computes the logarithm of the double factorial of n, \log(n!!).
    #[doc(alias = "gsl_sf_lndoublefact_e")]
    pub fn lndoublefact_e(n: u32) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lndoublefact_e(n, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the combinatorial factor n choose m = n!/(m!(n-m)!)
    #[doc(alias = "gsl_sf_choose")]
    pub fn choose(n: u32, m: u32) -> f64 {
        unsafe { sys::gsl_sf_choose(n, m) }
    }

    /// This routine computes the combinatorial factor n choose m = n!/(m!(n-m)!)
    #[doc(alias = "gsl_sf_choose_e")]
    pub fn choose_e(n: u32, m: u32) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_choose_e(n, m, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the logarithm of n choose m. This is equivalent to the sum \log(n!) - \log(m!) - \log((n-m)!).
    #[doc(alias = "gsl_sf_lnchoose")]
    pub fn lnchoose(n: u32, m: u32) -> f64 {
        unsafe { sys::gsl_sf_lnchoose(n, m) }
    }

    /// This routine computes the logarithm of n choose m. This is equivalent to the sum \log(n!) - \log(m!) - \log((n-m)!).
    #[doc(alias = "gsl_sf_lnchoose_e")]
    pub fn lnchoose_e(n: u32, m: u32) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lnchoose_e(n, m, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the Taylor coefficient x^n / n! for x >= 0, n >= 0.
    #[doc(alias = "gsl_sf_taylorcoeff")]
    pub fn taylorcoeff(n: i32, x: f64) -> f64 {
        unsafe { sys::gsl_sf_taylorcoeff(n, x) }
    }

    /// This routine computes the Taylor coefficient x^n / n! for x >= 0, n >= 0.
    #[doc(alias = "gsl_sf_taylorcoeff_e")]
    pub fn taylorcoeff_e(n: i32, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_taylorcoeff_e(n, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }
}

pub mod pochhammer_symbol {
    use crate::{types, Error};
    use std::mem::MaybeUninit;

    /// This routine computes the Pochhammer symbol (a)_x = \Gamma(a + x)/\Gamma(a).
    /// The Pochhammer symbol is also known as the Apell symbol and sometimes written as (a,x).
    /// When a and a+x are negative integers or zero, the limiting value of the ratio is returned.
    #[doc(alias = "gsl_sf_poch")]
    pub fn poch(a: f64, x: f64) -> f64 {
        unsafe { sys::gsl_sf_poch(a, x) }
    }

    /// This routine computes the Pochhammer symbol (a)_x = \Gamma(a + x)/\Gamma(a).
    /// The Pochhammer symbol is also known as the Apell symbol and sometimes written as (a,x).
    /// When a and a+x are negative integers or zero, the limiting value of the ratio is returned.
    #[doc(alias = "gsl_sf_poch_e")]
    pub fn poch_e(a: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_poch_e(a, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the logarithm of the Pochhammer symbol, \log((a)_x) = \log(\Gamma(a + x)/\Gamma(a)).
    #[doc(alias = "gsl_sf_lnpoch")]
    pub fn lnpoch(a: f64, x: f64) -> f64 {
        unsafe { sys::gsl_sf_lnpoch(a, x) }
    }

    /// This routine computes the logarithm of the Pochhammer symbol, \log((a)_x) = \log(\Gamma(a + x)/\Gamma(a)).
    #[doc(alias = "gsl_sf_lnpoch_e")]
    pub fn lnpoch_e(a: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lnpoch_e(a, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// These routines compute the sign of the Pochhammer symbol and the logarithm of its magnitude.
    /// The computed parameters are result = \log(|(a)_x|) with a corresponding error term, and sgn = \sgn((a)_x) where (a)_x = \Gamma(a + x)/\Gamma(a).
    #[doc(alias = "gsl_sf_lnpoch_sgn_e")]
    pub fn lnpoch_sgn_e(a: f64, x: f64, sgn: &mut f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lnpoch_sgn_e(a, x, result.as_mut_ptr(), sgn) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the relative Pochhammer symbol ((a)_x - 1)/x where (a)_x = \Gamma(a + x)/\Gamma(a).
    #[doc(alias = "gsl_sf_pochrel")]
    pub fn pochrel(a: f64, x: f64) -> f64 {
        unsafe { sys::gsl_sf_pochrel(a, x) }
    }

    /// This routine computes the relative Pochhammer symbol ((a)_x - 1)/x where (a)_x = \Gamma(a + x)/\Gamma(a).
    #[doc(alias = "gsl_sf_pochrel_e")]
    pub fn pochrel_e(a: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_pochrel_e(a, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }
}

pub mod beta {
    use crate::{types, Error};
    use std::mem::MaybeUninit;

    /// This routine computes the Beta Function, B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b) subject to a and b not being negative integers.
    #[doc(alias = "gsl_sf_beta")]
    pub fn beta(a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_sf_beta(a, b) }
    }

    /// This routine computes the Beta Function, B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b) subject to a and b not being negative integers.
    #[doc(alias = "gsl_sf_beta_e")]
    pub fn beta_e(a: f64, b: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_beta_e(a, b, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the logarithm of the Beta Function, \log(B(a,b)) subject to a and b not being negative integers.
    #[doc(alias = "gsl_sf_lnbeta")]
    pub fn lnbeta(a: f64, b: f64) -> f64 {
        unsafe { sys::gsl_sf_lnbeta(a, b) }
    }

    /// This routine computes the logarithm of the Beta Function, \log(B(a,b)) subject to a and b not being negative integers.
    #[doc(alias = "gsl_sf_lnbeta_e")]
    pub fn lnbeta_e(a: f64, b: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_lnbeta_e(a, b, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }
}

pub mod incomplete_gamma {
    use crate::{types, Error};
    use std::mem::MaybeUninit;

    /// This routine computes the unnormalized incomplete Gamma Function \Gamma(a,x) = \int_x^\infty dt t^{a-1} \exp(-t) for a real and x >= 0.
    #[doc(alias = "gsl_sf_gamma_inc")]
    pub fn gamma_inc(a: f64, x: f64) -> f64 {
        unsafe { sys::gsl_sf_gamma_inc(a, x) }
    }

    /// This routine computes the unnormalized incomplete Gamma Function \Gamma(a,x) = \int_x^\infty dt t^{a-1} \exp(-t) for a real and x >= 0.
    #[doc(alias = "gsl_sf_gamma_inc_e")]
    pub fn gamma_inc_e(a: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_gamma_inc_e(a, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the normalized incomplete Gamma Function Q(a,x) = 1/\Gamma(a) \int_x^\infty dt t^{a-1} \exp(-t) for a > 0, x >= 0.
    #[doc(alias = "gsl_sf_gamma_inc_Q")]
    pub fn gamma_inc_Q(a: f64, x: f64) -> f64 {
        unsafe { sys::gsl_sf_gamma_inc_Q(a, x) }
    }

    /// This routine computes the normalized incomplete Gamma Function Q(a,x) = 1/\Gamma(a) \int_x^\infty dt t^{a-1} \exp(-t) for a > 0, x >= 0.
    #[doc(alias = "gsl_sf_gamma_inc_Q_e")]
    pub fn gamma_inc_Q_e(a: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_gamma_inc_Q_e(a, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the complementary normalized incomplete Gamma Function P(a,x) = 1 - Q(a,x) = 1/\Gamma(a) \int_0^x dt t^{a-1} \exp(-t) for a > 0, x >= 0.
    ///
    /// Note that Abramowitz & Stegun call P(a,x) the incomplete gamma function (section 6.5).
    #[doc(alias = "gsl_sf_gamma_inc_P")]
    pub fn gamma_inc_P(a: f64, x: f64) -> f64 {
        unsafe { sys::gsl_sf_gamma_inc_P(a, x) }
    }

    /// This routine computes the complementary normalized incomplete Gamma Function P(a,x) = 1 - Q(a,x) = 1/\Gamma(a) \int_0^x dt t^{a-1} \exp(-t) for a > 0, x >= 0.
    ///
    /// Note that Abramowitz & Stegun call P(a,x) the incomplete gamma function (section 6.5).
    #[doc(alias = "gsl_sf_gamma_inc_P_e")]
    pub fn gamma_inc_P_e(a: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_gamma_inc_P_e(a, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }
}

pub mod incomplete_beta {
    use crate::{types, Error};
    use std::mem::MaybeUninit;

    /// This routine computes the normalized incomplete Beta function I_x(a,b)=B_x(a,b)/B(a,b) where B_x(a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt for 0 <= x <= 1.
    /// For a > 0, b > 0 the value is computed using a continued fraction expansion.
    /// For all other values it is computed using the relation I_x(a,b,x) = (1/a) x^a 2F1(a,1-b,a+1,x)/B(a,b).
    #[doc(alias = "gsl_sf_beta_inc")]
    pub fn beta_inc(a: f64, b: f64, x: f64) -> f64 {
        unsafe { sys::gsl_sf_beta_inc(a, b, x) }
    }

    /// This routine computes the normalized incomplete Beta function I_x(a,b)=B_x(a,b)/B(a,b) where B_x(a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt for 0 <= x <= 1.
    /// For a > 0, b > 0 the value is computed using a continued fraction expansion.
    /// For all other values it is computed using the relation I_x(a,b,x) = (1/a) x^a 2F1(a,1-b,a+1,x)/B(a,b).
    #[doc(alias = "gsl_sf_beta_inc_e")]
    pub fn beta_inc_e(a: f64, b: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_beta_inc_e(a, b, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }
}
