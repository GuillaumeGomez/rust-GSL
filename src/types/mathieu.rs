//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The routines described in this section compute the angular and radial Mathieu functions, and their characteristic values. Mathieu functions are the solutions of the following two differential equations:

d^2y/dv^2 + (a - 2q\cos 2v)y = 0

d^2f/du^2 - (a - 2q\cosh 2u)f = 0

The angular Mathieu functions ce_r(x,q), se_r(x,q) are the even and odd periodic solutions of the first equation, which is known as Mathieu’s equation. These exist only for the discrete sequence of characteristic values a=a_r(q) (even-periodic) and a=b_r(q) (odd-periodic).

The radial Mathieu functions Mc^{(j)}_{r}(z,q), Ms^{(j)}_{r}(z,q) are the solutions of the second equation, which is referred to as Mathieu’s modified equation. The radial Mathieu functions of the first, second, third and fourth kind are denoted by the parameter j, which takes the value 1, 2, 3 or 4.

For more information on the Mathieu functions, see Abramowitz and Stegun, Chapter 20.
!*/

use crate::ffi::FFI;
use crate::{types, Error};
use std::mem::MaybeUninit;

ffi_wrapper!(MathieuWorkspace, *mut sys::gsl_sf_mathieu_workspace, gsl_sf_mathieu_free,
"The Mathieu functions can be computed for a single order or for multiple orders, using array-based
routines. The array-based routines require a preallocated workspace.");

impl MathieuWorkspace {
    /// This function returns a workspace for the array versions of the Mathieu routines.
    /// The arguments n and qmax specify the maximum order and q-value of Mathieu functions which can be computed with this workspace.
    #[doc(alias = "gsl_sf_mathieu_alloc")]
    pub fn new(n: usize, qmax: f64) -> Option<MathieuWorkspace> {
        let tmp = unsafe { sys::gsl_sf_mathieu_alloc(n, qmax) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This routine computes the characteristic values a_n(q), b_n(q) of the Mathieu functions ce_n(q,x) and se_n(q,x), respectively.
    #[doc(alias = "gsl_sf_mathieu_a_e")]
    pub fn mathieu_a(n: i32, q: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_a_e(n, q, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the characteristic values a_n(q), b_n(q) of the Mathieu functions ce_n(q,x) and se_n(q,x), respectively.
    #[doc(alias = "gsl_sf_mathieu_b_e")]
    pub fn mathieu_b(n: i32, q: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_b_e(n, q, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes a series of Mathieu characteristic values a_n(q), b_n(q) for n from order_min to order_max inclusive, storing the results in the array result_array.
    #[doc(alias = "gsl_sf_mathieu_a_array")]
    pub fn mathieu_a_array(
        &mut self,
        order_min: i32,
        order_max: i32,
        q: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_sf_mathieu_a_array(
                order_min,
                order_max,
                q,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// This routine computes a series of Mathieu characteristic values a_n(q), b_n(q) for n from order_min to order_max inclusive, storing the results in the array result_array.
    #[doc(alias = "gsl_sf_mathieu_b_array")]
    pub fn mathieu_b_array(
        &mut self,
        order_min: i32,
        order_max: i32,
        q: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_sf_mathieu_b_array(
                order_min,
                order_max,
                q,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// This routine computes the angular Mathieu functions ce_n(q,x) and se_n(q,x), respectively.
    #[doc(alias = "gsl_sf_mathieu_ce_e")]
    pub fn mathieu_ce(n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_ce_e(n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the angular Mathieu functions ce_n(q,x) and se_n(q,x), respectively.
    #[doc(alias = "gsl_sf_mathieu_se_e")]
    pub fn mathieu_se(n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_se_e(n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes a series of the angular Mathieu functions ce_n(q,x) and se_n(q,x) of order n from nmin to nmax inclusive, storing the results in the array result_array.
    #[doc(alias = "gsl_sf_mathieu_ce_array")]
    pub fn mathieu_ce_array(
        &mut self,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_sf_mathieu_ce_array(
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// This routine computes a series of the angular Mathieu functions ce_n(q,x) and se_n(q,x) of order n from nmin to nmax inclusive, storing the results in the array result_array.
    #[doc(alias = "gsl_sf_mathieu_se_array")]
    pub fn mathieu_se_array(
        &mut self,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_sf_mathieu_se_array(
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// This routine computes the radial j-th kind Mathieu functions Mc_n^{(j)}(q,x) and Ms_n^{(j)}(q,x) of order n.
    ///
    /// The allowed values of j are 1 and 2. The functions for j = 3,4 can be computed as M_n^{(3)} = M_n^{(1)} + iM_n^{(2)} and M_n^{(4)} = M_n^{(1)} - iM_n^{(2)}, where M_n^{(j)} = Mc_n^{(j)} or Ms_n^{(j)}.
    #[doc(alias = "gsl_sf_mathieu_Mc_e")]
    pub fn mathieu_Mc(j: i32, n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_Mc_e(j, n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes the radial j-th kind Mathieu functions Mc_n^{(j)}(q,x) and Ms_n^{(j)}(q,x) of order n.
    ///
    /// The allowed values of j are 1 and 2. The functions for j = 3,4 can be computed as M_n^{(3)} = M_n^{(1)} + iM_n^{(2)} and M_n^{(4)} = M_n^{(1)} - iM_n^{(2)}, where M_n^{(j)} = Mc_n^{(j)} or Ms_n^{(j)}.
    #[doc(alias = "gsl_sf_mathieu_Ms_e")]
    pub fn mathieu_Ms(j: i32, n: i32, q: f64, x: f64) -> Result<types::Result, Error> {
        let mut result = MaybeUninit::<sys::gsl_sf_result>::uninit();
        let ret = unsafe { sys::gsl_sf_mathieu_Ms_e(j, n, q, x, result.as_mut_ptr()) };

        Error::handle(ret, unsafe { result.assume_init() }.into())
    }

    /// This routine computes a series of the radial Mathieu functions of kind j, with order from nmin to nmax inclusive, storing the results in the array result_array.
    #[doc(alias = "gsl_sf_mathieu_Mc_array")]
    pub fn mathieu_Mc_array(
        &mut self,
        j: i32,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_sf_mathieu_Mc_array(
                j,
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }

    /// This routine computes a series of the radial Mathieu functions of kind j, with order from nmin to nmax inclusive, storing the results in the array result_array.
    #[doc(alias = "gsl_sf_mathieu_Ms_array")]
    pub fn mathieu_Ms_array(
        &mut self,
        j: i32,
        nmin: i32,
        nmax: i32,
        q: f64,
        x: f64,
        result_array: &mut [f64],
    ) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_sf_mathieu_Ms_array(
                j,
                nmin,
                nmax,
                q,
                x,
                self.unwrap_unique(),
                result_array.as_mut_ptr(),
            )
        };
        Error::handle(ret, ())
    }
}
