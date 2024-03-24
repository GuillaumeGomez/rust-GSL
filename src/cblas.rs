//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use crate::vector::Vector;

/// Return the length of `x` as a `i32` value (to use in CBLAS calls).
#[inline]
fn len<F, T: Vector<F>>(x: &T) -> i32 {
    x.len().try_into().expect("Length must fit in `i32`")
}

#[inline]
fn as_ptr<F, T: Vector<F>>(x: &T) -> *const F {
    x.as_slice().as_ptr()
}

#[inline]
fn as_mut_ptr<F, T: Vector<F>>(x: &mut T) -> *mut F {
    x.as_mut_slice().as_mut_ptr()
}

/// Return the stride of `x` as a `i32` value (to use in CBLAS calls).
#[inline]
fn stride<F, T: Vector<F>>(x: &T) -> i32 {
    x.stride().try_into().expect("Stride must fit in `i32`")
}

pub mod level1 {
    use super::{as_mut_ptr, as_ptr, len, stride};
    use crate::vector::{check_equal_len, Vector};
    #[cfg(feature = "complex")]
    use num_complex::Complex;

    /// Return the sum of `alpha` and the dot product of `x` and `y`.
    #[doc(alias = "cblas_sdsdot")]
    pub fn sdsdot<T: Vector<f32>>(alpha: f32, x: &T, y: &T) -> f32 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_sdsdot(len(x), alpha, as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }

    /// Return the dot product of `x` and `y`.
    #[doc(alias = "cblas_dsdot")]
    pub fn dsdot<T: Vector<f32>>(x: &T, y: &T) -> f64 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_dsdot(len(x), as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }

    /// Return the dot product of `x` and `y`.
    #[doc(alias = "cblas_sdot")]
    pub fn sdot<T: Vector<f32>>(x: &T, y: &T) -> f32 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_sdot(len(x), as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }

    /// Return the dot product of `x` and `y`.
    #[doc(alias = "cblas_ddot")]
    pub fn ddot<T: Vector<f64>>(x: &T, y: &T) -> f64 {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        unsafe { sys::cblas_ddot(len(x), as_ptr(x), stride(x), as_ptr(y), stride(y)) }
    }

    #[cfg(feature = "complex")]
    /// Return the unconjugated dot product between `x` and `y`, that
    /// is ∑ xᵢ yᵢ.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::level1::cdotu;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(cdotu(&x, &x), Complex::new(3., 6.))
    /// ```
    #[doc(alias = "cblas_cdotu_sub")]
    pub fn cdotu<T: Vector<Complex<f32>>>(x: &T, y: &T) -> Complex<f32> {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotu: Complex<f32> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_cdotu_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotu as *mut Complex<f32> as *mut _,
            )
        }
        dotu
    }

    #[cfg(feature = "complex")]
    /// Return the (conjugated) dot product between `x` and `y`, that
    /// is ∑ x̅ᵢ yᵢ.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::level1::cdotc;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(cdotc(&x, &x), Complex::new(7., 0.))
    /// ```
    #[doc(alias = "cblas_cdotc_sub")]
    pub fn cdotc<T: Vector<Complex<f32>>>(x: &T, y: &T) -> Complex<f32> {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotc: Complex<f32> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_cdotc_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotc as *mut Complex<f32> as *mut _,
            )
        }
        dotc
    }

    #[cfg(feature = "complex")]
    /// Return the unconjugated dot product between `x` and `y`, that
    /// is ∑ xᵢ yᵢ.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::level1::zdotu;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(zdotu(&x, &x), Complex::new(3., 6.))
    /// ```
    #[doc(alias = "cblas_zdotu_sub")]
    pub fn zdotu<T: Vector<Complex<f64>>>(x: &T, y: &T) -> Complex<f64> {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotu: Complex<f64> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_zdotu_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotu as *mut Complex<f64> as *mut _,
            )
        }
        dotu
    }

    #[cfg(feature = "complex")]
    /// Return the (conjugated) dot product between `x` and `y`, that
    /// is ∑ x̅ᵢ yᵢ.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::level1::zdotc;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(zdotc(&x, &x), Complex::new(7., 0.))
    /// ```
    #[doc(alias = "cblas_zdotc_sub")]
    pub fn zdotc<T: Vector<Complex<f64>>>(x: &T, y: &T) -> Complex<f64> {
        check_equal_len(x, y).expect("The length of `x` and `y` must be equal");
        let mut dotc: Complex<f64> = Complex::new(0., 0.);
        unsafe {
            sys::cblas_zdotc_sub(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_ptr(y) as *const _,
                stride(y),
                &mut dotc as *mut Complex<f64> as *mut _,
            )
        }
        dotc
    }

    /// Return the Euclidean norm of `x`.
    #[doc(alias = "cblas_snrm2")]
    pub fn snrm2<T: Vector<f32>>(x: &T) -> f32 {
        unsafe { sys::cblas_snrm2(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the sum of the absolute values of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_sasum")]
    pub fn sasum<T: Vector<f32>>(x: &T) -> f32 {
        unsafe { sys::cblas_sasum(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the Euclidean norm of `x`.
    #[doc(alias = "cblas_dnrm2")]
    pub fn dnrm2<T: Vector<f64>>(x: &T) -> f64 {
        unsafe { sys::cblas_dnrm2(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the sum of the absolute values of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_dasum")]
    pub fn dasum<T: Vector<f64>>(x: &T) -> f64 {
        unsafe { sys::cblas_dasum(len(x), as_ptr(x), stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Return the Euclidean norm of `x`.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::level1::scnrm2;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(scnrm2(&x), 7f32.sqrt())
    /// ```
    #[doc(alias = "cblas_scnrm2")]
    pub fn scnrm2<T: Vector<Complex<f32>>>(x: &T) -> f32 {
        unsafe { sys::cblas_scnrm2(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Return the sum of the modulus of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_scasum")]
    pub fn scasum<T: Vector<Complex<f32>>>(x: &T) -> f32 {
        unsafe { sys::cblas_scasum(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Return the Euclidean norm of `x`.
    ///
    /// # Example
    ///
    /// ```
    /// use num_complex::Complex;
    /// use rgsl::cblas::level1::dznrm2;
    /// let x = [Complex::new(1., 1.), Complex::new(2., 1.)];
    /// assert_eq!(dznrm2(&x), 7f64.sqrt())
    /// ```
    #[doc(alias = "cblas_dznrm2")]
    pub fn dznrm2<T: Vector<Complex<f64>>>(x: &T) -> f64 {
        unsafe { sys::cblas_dznrm2(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Return the sum of the modulus of the elements of `x`
    /// (i.e., its L¹-norm).
    #[doc(alias = "cblas_dzasum")]
    pub fn dzasum<T: Vector<Complex<f64>>>(x: &T) -> f64 {
        unsafe { sys::cblas_dzasum(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Return the index of the element with maximum absolute value.
    #[doc(alias = "cblas_isamax")]
    pub fn isamax<T: Vector<f32>>(x: &T) -> usize {
        unsafe { sys::cblas_isamax(len(x), as_ptr(x), stride(x)) }
    }

    /// Return the index of the element with maximum absolute value.
    #[doc(alias = "cblas_idamax")]
    pub fn idamax<T: Vector<f64>>(x: &T) -> usize {
        unsafe { sys::cblas_idamax(len(x), as_ptr(x), stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Return the index of the element with maximum modulus.
    #[doc(alias = "cblas_icamax")]
    pub fn icamax<T: Vector<Complex<f32>>>(x: &T) -> usize {
        unsafe { sys::cblas_icamax(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Return the index of the element with maximum modulus.
    #[doc(alias = "cblas_izamax")]
    pub fn izamax<T: Vector<Complex<f64>>>(x: &T) -> usize {
        unsafe { sys::cblas_izamax(len(x), as_ptr(x) as *const _, stride(x)) }
    }

    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_sswap")]
    pub fn sswap<T: Vector<f32>>(x: &mut T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_sswap(len(x), as_mut_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_scopy")]
    pub fn scopy<T: Vector<f32>>(x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_scopy(len(x), as_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_saxpy")]
    pub fn saxpy<T: Vector<f32>>(alpha: f32, x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_saxpy(
                len(x),
                alpha,
                as_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
            )
        }
    }

    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_dswap")]
    pub fn dswap<T: Vector<f64>>(x: &mut T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_dswap(len(x), as_mut_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_dcopy")]
    pub fn dcopy<T: Vector<f64>>(x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe { sys::cblas_dcopy(len(x), as_ptr(x), stride(x), as_mut_ptr(y), stride(y)) }
    }

    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_daxpy")]
    pub fn daxpy<T: Vector<f64>>(alpha: f64, x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_daxpy(
                len(x),
                alpha,
                as_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_cswap")]
    pub fn cswap<T: Vector<Complex<f32>>>(x: &mut T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_cswap(
                len(x),
                as_mut_ptr(x) as *mut _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_ccopy")]
    pub fn ccopy<T: Vector<Complex<f32>>>(x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_ccopy(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_caxpy")]
    pub fn caxpy<T: Vector<Complex<f32>>>(alpha: &Complex<f32>, x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_caxpy(
                len(x),
                alpha as *const Complex<f32> as *const _,
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// Swap vectors `x` and `y`.
    #[doc(alias = "cblas_zswap")]
    pub fn zswap<T: Vector<Complex<f64>>>(x: &mut T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_zswap(
                len(x),
                as_mut_ptr(x) as *mut _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// Copy the content of `x` into `y`.
    #[doc(alias = "cblas_zcopy")]
    pub fn zcopy<T: Vector<Complex<f64>>>(x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_zcopy(
                len(x),
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// `y` := `alpha` * `x` + `y`.
    #[doc(alias = "cblas_zaxpy")]
    pub fn zaxpy<T: Vector<Complex<f64>>>(alpha: &Complex<f64>, x: &T, y: &mut T) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_zaxpy(
                len(x),
                alpha as *const Complex<f64> as *const _,
                as_ptr(x) as *const _,
                stride(x),
                as_mut_ptr(y) as *mut _,
                stride(y),
            )
        }
    }

    /// Given the Cartesian coordinates (`a`, `b`), returns
    /// (c, s, r, z) such that
    ///
    /// ⎧c  s⎫ ⎧a⎫ = ⎧r⎫
    /// ⎩s  c⎭ ⎩b⎭   ⎩0⎭
    ///
    /// The value of z is defined such that if |`a`| > |`b`|, z is s;
    /// otherwise if c ≠ 0, z is 1/c; otherwise z is 1.
    #[doc(alias = "cblas_srotg")]
    pub fn srotg(a: f32, b: f32) -> (f32, f32, f32, f32) {
        let mut c = f32::NAN;
        let mut s = f32::NAN;
        let mut r = a;
        let mut z = b;
        unsafe {
            sys::cblas_srotg(
                &mut r as *mut _,
                &mut z as *mut _,
                &mut c as *mut _,
                &mut s as *mut _,
            )
        }
        (c, s, r, z)
    }

    /// Modified matrix transformation (for the mathematical field `F`).
    #[derive(Clone, Copy)]
    pub enum H<F> {
        /// Specify that H is the matrix
        ///
        /// ⎧`h11`  `h12`⎫
        /// ⎩`h21`  `h22`⎭
        Full {
            h11: F,
            h21: F,
            h12: F,
            h22: F,
        },
        /// Specify that H is the matrix
        ///
        /// ⎧1.0  `h12`⎫
        /// ⎩`h21`  1.0⎭
        OffDiag {
            h21: F,
            h12: F,
        },
        /// Specify that H is the matrix
        ///
        /// ⎧`h11`   1.0⎫
        /// ⎩-1.0  `h22`⎭
        Diag {
            h11: F,
            h22: F,
        },
        Id,
    }

    /// Given Cartesian coordinates (`x1`, `x2`), return the
    /// transformation matrix H that zeros the second component or the
    /// vector (`x1` √`d1`, `x2` √`d2`):
    ///
    /// H ⎧`x1` √`d1`⎫ = ⎧y1⎫
    ///   ⎩`x2` √`d2`⎭   ⎩0.⎭
    ///
    /// The second component of the return value is `y1`.
    #[doc(alias = "cblas_srotmg")]
    pub fn srotmg(mut d1: f32, mut d2: f32, mut x1: f32, x2: f32) -> (H<f32>, f32) {
        let mut h: [f32; 5] = [0.; 5];
        unsafe {
            sys::cblas_srotmg(
                &mut d1 as *mut _,
                &mut d2 as *mut _,
                &mut x1 as *mut _,
                x2,
                &mut h as *mut _,
            )
        }
        let h = match h[0] {
            -1.0 => H::Full {
                h11: h[1],
                h21: h[2],
                h12: h[3],
                h22: h[4],
            },
            0.0 => H::OffDiag {
                h21: h[2],
                h12: h[3],
            },
            1.0 => H::Diag {
                h11: h[1],
                h22: h[4],
            },
            -2.0 => H::Id,
            _ => unreachable!("srotmg: incorrect flag value"),
        };
        (h, x1)
    }

    /// Apply plane rotation.  More specifically, perform the
    /// following transformation in place :
    ///
    /// ⎧`x`ᵢ⎫ = ⎧`c`  `s`⎫ ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭   ⎩-`s` `c`⎭ ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_srot")]
    pub fn srot<T: Vector<f32>>(x: &mut T, y: &mut T, c: f32, s: f32) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_srot(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                c,
                s,
            )
        }
    }

    /// Apply the matrix rotation `h` to `x`, `y`.
    ///
    /// ⎧`x`ᵢ⎫ = `h` ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭       ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_srotm")]
    pub fn srotm<T: Vector<f32>>(x: &mut T, y: &mut T, h: H<f32>) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        let p = match h {
            H::Full { h11, h21, h12, h22 } => [-1.0, h11, h21, h12, h22],
            H::OffDiag { h21, h12 } => [0.0, 1., h21, h12, 1.],
            H::Diag { h11, h22 } => [1.0, h11, -1., 1., h22],
            H::Id => [-2.0, 1., 0., 0., 1.],
        };
        unsafe {
            sys::cblas_srotm(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                &p as *const _,
            )
        }
    }

    /// Given the Cartesian coordinates (`a`, `b`), returns
    /// (c, s, r, z) such that
    ///
    /// ⎧c  s⎫ ⎧a⎫ = ⎧r⎫
    /// ⎩s  c⎭ ⎩b⎭   ⎩0⎭
    ///
    /// The value of z is defined such that if |`a`| > |`b`|, z is s;
    /// otherwise if c ≠ 0, z is 1/c; otherwise z is 1.
    #[doc(alias = "cblas_drotg")]
    pub fn drotg(a: f64, b: f64) -> (f64, f64, f64, f64) {
        let mut c = f64::NAN;
        let mut s = f64::NAN;
        let mut r = a;
        let mut z = b;
        unsafe {
            sys::cblas_drotg(
                &mut r as *mut _,
                &mut z as *mut _,
                &mut c as *mut _,
                &mut s as *mut _,
            )
        }
        (c, s, r, z)
    }

    /// Given Cartesian coordinates (`x1`, `x2`), return the
    /// transformation matrix H that zeros the second component or the
    /// vector (`x1` √`d1`, `x2` √`d2`):
    ///
    /// H ⎧`x1` √`d1`⎫ = ⎧y1⎫
    ///   ⎩`x2` √`d2`⎭   ⎩0.⎭
    ///
    /// The second component of the return value is `y1`.
    #[doc(alias = "cblas_drotmg")]
    pub fn drotmg(mut d1: f64, mut d2: f64, mut x1: f64, x2: f64) -> (H<f64>, f64) {
        let mut h: [f64; 5] = [0.; 5];
        unsafe {
            sys::cblas_drotmg(
                &mut d1 as *mut _,
                &mut d2 as *mut _,
                &mut x1 as *mut _,
                x2,
                &mut h as *mut _,
            )
        }
        let h = match h[0] {
            -1.0 => H::Full {
                h11: h[1],
                h21: h[2],
                h12: h[3],
                h22: h[4],
            },
            0.0 => H::OffDiag {
                h21: h[2],
                h12: h[3],
            },
            1.0 => H::Diag {
                h11: h[1],
                h22: h[4],
            },
            -2.0 => H::Id,
            _ => unreachable!("srotmg: incorrect flag value"),
        };
        (h, x1)
    }

    /// Apply plane rotation.  More specifically, perform the
    /// following transformation in place :
    ///
    /// ⎧`x`ᵢ⎫ = ⎧`c`  `s`⎫ ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭   ⎩-`s` `c`⎭ ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_drot")]
    pub fn drot<T: Vector<f64>>(x: &mut T, y: &mut T, c: f64, s: f64) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        unsafe {
            sys::cblas_drot(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                c,
                s,
            )
        }
    }

    /// Apply the matrix rotation `h` to `x`, `y`.
    ///
    /// ⎧`x`ᵢ⎫ = `h` ⎧`x`ᵢ⎫
    /// ⎩`y`ᵢ⎭       ⎩`y`ᵢ⎭
    ///
    /// for all indices i.
    #[doc(alias = "cblas_drotm")]
    pub fn drotm<T: Vector<f64>>(x: &mut T, y: &mut T, h: H<f64>) {
        check_equal_len(x, y).expect("Vectors `x` and `y` must have the same length");
        let p = match h {
            H::Full { h11, h21, h12, h22 } => [-1.0, h11, h21, h12, h22],
            H::OffDiag { h21, h12 } => [0.0, 1., h21, h12, 1.],
            H::Diag { h11, h22 } => [1.0, h11, -1., 1., h22],
            H::Id => [-2.0, 1., 0., 0., 1.],
        };
        unsafe {
            sys::cblas_drotm(
                len(x),
                as_mut_ptr(x),
                stride(x),
                as_mut_ptr(y),
                stride(y),
                &p as *const _,
            )
        }
    }

    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_sscal")]
    pub fn sscal<T: Vector<f32>>(alpha: f32, x: &mut T) {
        unsafe { sys::cblas_sscal(len(x), alpha, as_mut_ptr(x), stride(x)) }
    }

    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_dscal")]
    pub fn dscal<T: Vector<f64>>(alpha: f64, x: &mut T) {
        unsafe { sys::cblas_dscal(len(x), alpha, as_mut_ptr(x), stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_cscal")]
    pub fn cscal<T: Vector<Complex<f32>>>(alpha: &Complex<f32>, x: &mut T) {
        unsafe {
            sys::cblas_cscal(
                len(x),
                alpha as *const Complex<f32> as *const _,
                as_mut_ptr(x) as *mut _,
                stride(x),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_zscal")]
    pub fn zscal<T: Vector<Complex<f64>>>(alpha: &Complex<f64>, x: &mut T) {
        unsafe {
            sys::cblas_zscal(
                len(x),
                alpha as *const Complex<f64> as *const _,
                as_mut_ptr(x) as *mut _,
                stride(x),
            )
        }
    }

    #[cfg(feature = "complex")]
    /// Multiply each element of `x` by `alpha`.
    #[doc(alias = "cblas_csscal")]
    pub fn csscal<T: Vector<Complex<f32>>>(alpha: f32, x: &mut T) {
        unsafe { sys::cblas_csscal(len(x), alpha, as_mut_ptr(x) as *mut _, stride(x)) }
    }

    #[cfg(feature = "complex")]
    /// Multiple each element of a matrix/vector by a constant.
    #[doc(alias = "cblas_zdscal")]
    pub fn zdscal<T: Vector<Complex<f64>>>(alpha: f64, x: &mut T) {
        unsafe { sys::cblas_zdscal(len(x), alpha, as_mut_ptr(x) as *mut _, stride(x)) }
    }
}

pub mod level2 {
    use crate::enums;

    /// Multiplies a matrix and a vector.
    ///
    /// * order : Whether matrices are row major order (C-Style) for column major order (Fortran-style). One of enum CblasRowMajor or CblasColMajor
    /// * transA :  Whether to transpose matrix A. One of enum CblasNoTrans, CBlasTrans.
    /// * M : Rows in matrix A
    /// * N : Columns in matrix A
    /// * alpha : scalar factor for (sigma * op(A) * x)
    /// * A : matrix A
    /// * lda : The size of the first dimension of matrix A
    /// * X : vector X
    /// * incx : use every incX'th element of X
    /// * beta : scalar factor y
    /// * Y : vector Y
    /// * incy : use every incY'th element of Y
    ///
    /// For parameter lda, if you are passing a matrix `A[m][n]`, the value of parameter lda should be m.
    #[doc(alias = "cblas_sgemv")]
    pub fn sgemv(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        X: &[f32],
        incx: i32,
        beta: f32,
        Y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_sgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_sgbmv")]
    pub fn sgbmv(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        X: &[f32],
        incx: i32,
        beta: f32,
        Y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_sgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_strmv")]
    pub fn strmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_strmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stbmv")]
    pub fn stbmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stpmv")]
    pub fn stpmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[f32],
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_strsv")]
    pub fn strsv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_strsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stbsv")]
    pub fn stbsv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[f32],
        lda: i32,
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_stpsv")]
    pub fn stpsv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[f32],
        X: &mut [f32],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_stpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dgemv")]
    pub fn dgemv(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        X: &[f64],
        incx: i32,
        beta: f64,
        Y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dgbmv")]
    pub fn dgbmv(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        X: &[f64],
        incx: i32,
        beta: f64,
        Y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha,
                A.as_ptr(),
                lda,
                X.as_ptr(),
                incx,
                beta,
                Y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dtrmv")]
    pub fn dtrmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtrmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtbmv")]
    pub fn dtbmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtpmv")]
    pub fn dtpmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[f64],
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtrsv")]
    pub fn dtrsv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtrsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtbsv")]
    pub fn dtbsv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[f64],
        lda: i32,
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr(),
                lda,
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_dtpsv")]
    pub fn dtpsv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[f64],
        X: &mut [f64],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_dtpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr(),
                X.as_mut_ptr(),
                incx,
            )
        }
    }

    #[doc(alias = "cblas_cgemv")]
    pub fn cgemv<T>(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_cgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_cgbmv")]
    pub fn cgbmv<T>(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_cgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_ctrmv")]
    pub fn ctrmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctrmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctbmv")]
    pub fn ctbmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctpmv")]
    pub fn ctpmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctrsv")]
    pub fn ctrsv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctrsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctbsv")]
    pub fn ctbsv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ctpsv")]
    pub fn ctpsv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ctpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_zgemv")]
    pub fn zgemv<T>(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zgemv(
                order.into(),
                transA.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zgbmv")]
    pub fn zgbmv<T>(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        M: i32,
        N: i32,
        KL: i32,
        KU: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        X: &[T],
        incx: i32,
        beta: &[T],
        Y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zgbmv(
                order.into(),
                transA.into(),
                M,
                N,
                KL,
                KU,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                X.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                Y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_ztrmv")]
    pub fn ztrmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztrmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztbmv")]
    pub fn ztbmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztbmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztpmv")]
    pub fn ztpmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztpmv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztrsv")]
    pub fn ztrsv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztrsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztbsv")]
    pub fn ztbsv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        K: i32,
        A: &[T],
        lda: i32,
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztbsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                K,
                A.as_ptr() as *const _,
                lda,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ztpsv")]
    pub fn ztpsv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        N: i32,
        Ap: &[T],
        X: &mut [T],
        incx: i32,
    ) {
        unsafe {
            sys::cblas_ztpsv(
                order.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                N,
                Ap.as_ptr() as *const _,
                X.as_mut_ptr() as *mut _,
                incx,
            )
        }
    }

    #[doc(alias = "cblas_ssymv")]
    pub fn ssymv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        x: &[f32],
        incx: i32,
        beta: f32,
        y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_ssymv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_ssbmv")]
    pub fn ssbmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        x: &[f32],
        incx: i32,
        beta: f32,
        y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_ssbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_sspmv")]
    pub fn sspmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        Ap: &[f32],
        x: &[f32],
        incx: i32,
        beta: f32,
        y: &mut [f32],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_sspmv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                Ap.as_ptr(),
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_sger")]
    pub fn sger(
        order: enums::CblasOrder,
        M: i32,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        y: &[f32],
        incy: i32,
        A: &mut [f32],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_sger(
                order.into(),
                M,
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_ssyr")]
    pub fn ssyr(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        A: &mut [f32],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_ssyr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_sspr")]
    pub fn sspr(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        Ap: &mut [f32],
    ) {
        unsafe {
            sys::cblas_sspr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                Ap.as_mut_ptr(),
            )
        }
    }

    #[doc(alias = "cblas_ssyr2")]
    pub fn ssyr2(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        y: &[f32],
        incy: i32,
        A: &mut [f32],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_ssyr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_sspr2")]
    pub fn sspr2(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        x: &[f32],
        incx: i32,
        y: &[f32],
        incy: i32,
        A: &mut [f32],
    ) {
        unsafe {
            sys::cblas_sspr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
            )
        }
    }

    #[doc(alias = "cblas_dsymv")]
    pub fn dsymv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        x: &[f64],
        incx: i32,
        beta: f64,
        y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dsymv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dsbmv")]
    pub fn dsbmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        x: &[f64],
        incx: i32,
        beta: f64,
        y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dsbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dspmv")]
    pub fn dspmv(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        Ap: &[f64],
        x: &[f64],
        incx: i32,
        beta: f64,
        y: &mut [f64],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_dspmv(
                order.into(),
                uplo.into(),
                N,
                alpha,
                Ap.as_ptr(),
                x.as_ptr(),
                incx,
                beta,
                y.as_mut_ptr(),
                incy,
            )
        }
    }

    #[doc(alias = "cblas_dger")]
    pub fn dger(
        order: enums::CblasOrder,
        M: i32,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        y: &[f64],
        incy: i32,
        A: &mut [f64],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_dger(
                order.into(),
                M,
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_dsyr")]
    pub fn dsyr(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        A: &mut [f64],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_dsyr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_dspr")]
    pub fn dspr(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        Ap: &mut [f64],
    ) {
        unsafe {
            sys::cblas_dspr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                Ap.as_mut_ptr(),
            )
        }
    }

    #[doc(alias = "cblas_dsyr2")]
    pub fn dsyr2(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        y: &[f64],
        incy: i32,
        A: &mut [f64],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_dsyr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
                lda,
            )
        }
    }

    #[doc(alias = "cblas_dspr2")]
    pub fn dspr2(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        x: &[f64],
        incx: i32,
        y: &[f64],
        incy: i32,
        A: &mut [f64],
    ) {
        unsafe {
            sys::cblas_dspr2(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr(),
                incx,
                y.as_ptr(),
                incy,
                A.as_mut_ptr(),
            )
        }
    }

    #[doc(alias = "cblas_chemv")]
    pub fn chemv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_chemv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_chbmv")]
    pub fn chbmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_chbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_chpmv")]
    pub fn chpmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        Ap: &[T],
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_chpmv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                Ap.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_cgeru")]
    pub fn cgeru<T>(
        order: enums::CblasOrder,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cgeru(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_cgerc")]
    pub fn cgerc<T>(
        order: enums::CblasOrder,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cgerc(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_cher")]
    pub fn cher<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        x: &[T],
        incx: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cher(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_chpr")]
    pub fn chpr<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f32,
        x: &[T],
        incx: i32,
        Ap: &mut [T],
    ) {
        unsafe {
            sys::cblas_chpr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }

    #[doc(alias = "cblas_cher2")]
    pub fn cher2<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_cher2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_chpr2")]
    pub fn chpr2<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[f64],
        incy: i32,
        Ap: &mut [f64],
    ) {
        unsafe {
            sys::cblas_chpr2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }

    #[doc(alias = "cblas_zhemv")]
    pub fn zhemv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zhemv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zhbmv")]
    pub fn zhbmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zhbmv(
                order.into(),
                uplo.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zhpmv")]
    pub fn zhpmv<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        Ap: &[T],
        x: &[T],
        incx: i32,
        beta: &[T],
        y: &mut [T],
        incy: i32,
    ) {
        unsafe {
            sys::cblas_zhpmv(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                Ap.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                beta.as_ptr() as *const _,
                y.as_mut_ptr() as *mut _,
                incy,
            )
        }
    }

    #[doc(alias = "cblas_zgeru")]
    pub fn zgeru<T>(
        order: enums::CblasOrder,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zgeru(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zgerc")]
    pub fn zgerc<T>(
        order: enums::CblasOrder,
        M: i32,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zgerc(
                order.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zher")]
    pub fn zher<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        x: &[T],
        incx: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zher(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zhpr")]
    pub fn zhpr<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: f64,
        x: &[T],
        incx: i32,
        Ap: &mut [T],
    ) {
        unsafe {
            sys::cblas_zhpr(
                order.into(),
                uplo.into(),
                N,
                alpha,
                x.as_ptr() as *const _,
                incx,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }

    #[doc(alias = "cblas_zher2")]
    pub fn zher2<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[T],
        incy: i32,
        A: &mut [T],
        lda: i32,
    ) {
        unsafe {
            sys::cblas_zher2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                A.as_mut_ptr() as *mut _,
                lda,
            )
        }
    }

    #[doc(alias = "cblas_zhpr2")]
    pub fn zhpr2<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        N: i32,
        alpha: &[T],
        x: &[T],
        incx: i32,
        y: &[f64],
        incy: i32,
        Ap: &mut [f64],
    ) {
        unsafe {
            sys::cblas_zhpr2(
                order.into(),
                uplo.into(),
                N,
                alpha.as_ptr() as *const _,
                x.as_ptr() as *const _,
                incx,
                y.as_ptr() as *const _,
                incy,
                Ap.as_mut_ptr() as *mut _,
            )
        }
    }
}

pub mod level3 {
    use crate::enums;

    /// General crate::types::Matrix-MatrixF64 multiplication for single precision float.
    ///
    /// __Parameters:__
    ///
    /// * order : Whether matrices are row major order (C-Style) for column major order (Fortran-style). One of enum CblasRowMajor or CblasColMajor.
    /// * transA : Whether to transpose matrix A. One of enum CblasNoTrans, CBlasTrans, CBlasConjTrans.
    /// * transB : Whether to transpose matrix B. One of enum CblasNoTrans, CBlasTrans, CBlasConjTrans.
    /// * M : Rows in matrices A and C
    /// * N : Columns in Matrices B and C
    /// * K : Columns in matrix A and Rows in matrix B
    /// * alpha : scalar factor for op(A)op(B)
    /// * A : matrix A
    /// * lda : The size of the first dimension of matrix A
    /// * B : matrix B
    /// * ldb : The size of the first dimension of matrix B
    /// * beta : scalar factor for C
    /// * C : matrix C
    /// * ldc : The size of the first dimension of matrix C
    ///
    /// For parameters lda, ldb, and ldc, if you are passing a matrix `D[m][n]`, the value of parameter lda, ldb, or ldc should be m.
    #[doc(alias = "cblas_sgemm")]
    pub fn sgemm(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        M: i32,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &[f32],
        ldb: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_sgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    /// Symmetric crate::types::Matrix-MatrixF64 multiplication for single precision float.
    ///
    /// __Parameters:__
    ///
    /// * order : Whether matrices are row major order (C-Style) for column major order (Fortran-style). One of enum CblasRowMajor or CblasColMajor.
    /// * side : If CBlasSideLeft, perform (sigma(A)(B) + beta C). If CBlasSideRight, perform (sigma (B)(A) + beta C)
    /// * uplo : Indicates whether to use the upper (CBlasUpper) or lower (CBlasLower) triangle of matrix A
    /// * M : Rows in matrices A and C
    /// * N : Columns in Matrices B and C
    /// * alpha : scalar factor for op(A)op(B)
    /// * A : matrix A
    /// * lda : The size of the first dimension of matrix A
    /// * B : matrix B
    /// * ldb : The size of the first dimension of matrix B
    /// * beta : scalar factor for C
    /// * C : matrix C
    /// * ldc : The size of the first dimension of matrix C
    #[doc(alias = "cblas_ssymm")]
    pub fn ssymm(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &[f32],
        ldb: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_ssymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ssyrk")]
    pub fn ssyrk(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_ssyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ssyr2k")]
    pub fn ssyr2k(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &[f32],
        ldb: i32,
        beta: f32,
        C: &mut [f32],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_ssyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_strmm")]
    pub fn strmm(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &mut [f32],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_strmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_strsm")]
    pub fn strsm(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: f32,
        A: &[f32],
        lda: i32,
        B: &mut [f32],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_strsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_dgemm")]
    pub fn dgemm(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        M: i32,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &[f64],
        ldb: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dsymm")]
    pub fn dsymm(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &[f64],
        ldb: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dsymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dsyrk")]
    pub fn dsyrk(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dsyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dsyr2k")]
    pub fn dsyr2k(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &[f64],
        ldb: i32,
        beta: f64,
        C: &mut [f64],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_dsyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr(),
                lda,
                B.as_ptr(),
                ldb,
                beta,
                C.as_mut_ptr(),
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_dtrmm")]
    pub fn dtrmm(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &mut [f64],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_dtrmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_dtrsm")]
    pub fn dtrsm(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: f64,
        A: &[f64],
        lda: i32,
        B: &mut [f64],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_dtrsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha,
                A.as_ptr(),
                lda,
                B.as_mut_ptr(),
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_cgemm")]
    pub fn cgemm<T>(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        M: i32,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_cgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_csymm")]
    pub fn csymm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_csymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_csyrk")]
    pub fn csyrk<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_csyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_csyr2k")]
    pub fn csyr2k<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_csyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ctrmm")]
    pub fn ctrmm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ctrmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_ctrsm")]
    pub fn ctrsm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ctrsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_zgemm")]
    pub fn zgemm<T>(
        order: enums::CblasOrder,
        transA: enums::CblasTranspose,
        transB: enums::CblasTranspose,
        M: i32,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zgemm(
                order.into(),
                transA.into(),
                transB.into(),
                M,
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zsymm")]
    pub fn zsymm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zsymm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zsyrk")]
    pub fn zsyrk<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zsyrk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zsyr2k")]
    pub fn zsyr2k<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zsyr2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_ztrmm")]
    pub fn ztrmm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ztrmm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_ztrsm")]
    pub fn ztrsm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        transA: enums::CblasTranspose,
        diag: enums::CblasDiag,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &mut [T],
        ldb: i32,
    ) {
        unsafe {
            sys::cblas_ztrsm(
                order.into(),
                side.into(),
                uplo.into(),
                transA.into(),
                diag.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_mut_ptr() as *mut _,
                ldb,
            )
        }
    }

    #[doc(alias = "cblas_chemm")]
    pub fn chemm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_chemm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_cherk")]
    pub fn cherk<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: f32,
        A: &[T],
        lda: i32,
        beta: f32,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_cherk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr() as *const _,
                lda,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_cher2k")]
    pub fn cher2k<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: f32,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_cher2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zhemm")]
    pub fn zhemm<T>(
        order: enums::CblasOrder,
        side: enums::CblasSide,
        uplo: enums::CblasUplo,
        M: i32,
        N: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: &[T],
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zhemm(
                order.into(),
                side.into(),
                uplo.into(),
                M,
                N,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta.as_ptr() as *const _,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zherk")]
    pub fn zherk<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: f64,
        A: &[T],
        lda: i32,
        beta: f64,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zherk(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha,
                A.as_ptr() as *const _,
                lda,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }

    #[doc(alias = "cblas_zher2k")]
    pub fn zher2k<T>(
        order: enums::CblasOrder,
        uplo: enums::CblasUplo,
        trans: enums::CblasTranspose,
        N: i32,
        K: i32,
        alpha: &[T],
        A: &[T],
        lda: i32,
        B: &[T],
        ldb: i32,
        beta: f64,
        C: &mut [T],
        ldc: i32,
    ) {
        unsafe {
            sys::cblas_zher2k(
                order.into(),
                uplo.into(),
                trans.into(),
                N,
                K,
                alpha.as_ptr() as *const _,
                A.as_ptr() as *const _,
                lda,
                B.as_ptr() as *const _,
                ldb,
                beta,
                C.as_mut_ptr() as *mut _,
                ldc,
            )
        }
    }
}
