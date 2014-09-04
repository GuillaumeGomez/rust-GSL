//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use enums;
use std::intrinsics::{fabsf64, logf64, powf64};
use std::c_vec::CVec;

/// The QAG algorithm is a simple adaptive integration procedure. The integration region is divided into subintervals, and on each iteration
/// the subinterval with the largest estimated error is bisected. This reduces the overall error rapidly, as the subintervals become concentrated
/// around local difficulties in the integrand. These subintervals are managed by a gsl_integration_workspace struct, which handles the memory
/// for the subinterval ranges, results and error estimates.
pub struct IntegrationWorkspace {
    w: *mut ffi::gsl_integration_workspace
}

struct InternParam<'r, T:'r> {
    func: ::function<T>,
    param: &'r mut T,
    p2: f64
}

impl IntegrationWorkspace {
    /// This function allocates a workspace sufficient to hold n double precision intervals, their integration results and error estimates. One
    /// workspace may be used multiple times as all necessary reinitialization is performed automatically by the integration routines.
    pub fn new(n: u64) -> Option<IntegrationWorkspace> {
        let tmp = unsafe { ffi::gsl_integration_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(IntegrationWorkspace {
                w: tmp
            })
        }
    }

    /// This function applies an integration rule adaptively until an estimate of the integral of f over (a,b) is achieved within the desired
    /// absolute and relative error limits, epsabs and epsrel. The function returns the final approximation, result, and an estimate of the
    /// absolute error, abserr. The integration rule is determined by the value of key, which should be chosen from the following symbolic names,
    /// 
    /// GSL_INTEG_GAUSS15  (key = 1)
    /// 
    /// GSL_INTEG_GAUSS21  (key = 2)
    /// 
    /// GSL_INTEG_GAUSS31  (key = 3)
    /// 
    /// GSL_INTEG_GAUSS41  (key = 4)
    /// 
    /// GSL_INTEG_GAUSS51  (key = 5)
    /// 
    /// GSL_INTEG_GAUSS61  (key = 6)
    /// 
    /// corresponding to the 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod rules. The higher-order rules give better accuracy for smooth functions,
    /// while lower-order rules save time when the function contains local difficulties, such as discontinuities.
    /// 
    /// On each iteration the adaptive integration strategy bisects the interval with the largest error estimate. The subintervals and their
    /// results are stored in the memory provided by workspace. The maximum number of subintervals is given by limit, which may not exceed the
    /// allocated size of the workspace.
    pub fn qag<T>(&self, f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, key: ::GaussKonrodRule,
        result: &mut f64, abserr: &mut f64) -> enums::Value {
        let integration_rule = match key {
            enums::Gauss15 => ::integration::qk15,
            enums::Gauss21 => ::integration::qk21,
            enums::Gauss31 => ::integration::qk31,
            enums::Gauss41 => ::integration::qk41,
            enums::Gauss51 => ::integration::qk51,
            enums::Gauss61 => ::integration::qk61,
            /*_ => {
                let file = file!();
                "value of key does specify a known integration rule".with_c_str(|c_str|{
                    file.with_c_str(|c_file|{
                        unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::Inval as i32) }
                    });
                });
                // this line is not used but just for compilation...
                ::integration::qk15
            }*/
        };
        intern_qag(f, arg, a, b, epsabs, epsrel, limit, self, result, abserr, integration_rule)
    }

    /// This function applies the Gauss-Kronrod 21-point integration rule adaptively until an estimate of the integral of f over (a,b) is achieved
    /// within the desired absolute and relative error limits, epsabs and epsrel. The results are extrapolated using the epsilon-algorithm, which
    /// accelerates the convergence of the integral in the presence of discontinuities and integrable singularities. The function returns the
    /// final approximation from the extrapolation, result, and an estimate of the absolute error, abserr. The subintervals and their results are
    /// stored in the memory provided by workspace. The maximum number of subintervals is given by limit, which may not exceed the allocated size
    /// of the workspace.
    pub fn qags<T>(&self, f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, result: &mut f64,
        abserr: &mut f64) -> enums::Value {
        unsafe { intern_qags(f, arg, a, b, epsabs, epsrel, limit, self, result, abserr, ::integration::qk21) }
    }

    /// This function applies the adaptive integration algorithm QAGS taking account of the user-supplied locations of singular points. The array
    /// pts of length npts should contain the endpoints of the integration ranges defined by the integration region and locations of the singularities.
    /// For example, to integrate over the region (a,b) with break-points at x_1, x_2, x_3 (where a < x_1 < x_2 < x_3 < b) the following pts array
    /// should be used
    /// 
    /// pts[0] = a
    /// pts[1] = x_1
    /// pts[2] = x_2
    /// pts[3] = x_3
    /// pts[4] = b
    /// with npts = 5.
    /// 
    /// If you know the locations of the singular points in the integration region then this routine will be faster than QAGS.
    pub fn qagp<T>(&self, f: ::function<T>, arg: &mut T, pts: &mut [f64], epsabs: f64, epsrel: f64, limit: u64, result: &mut f64,
        abserr: &mut f64) -> enums::Value {
        unsafe { intern_qagp(f, arg, pts, epsabs, epsrel, limit, self, result, abserr, ::integration::qk21) }
    }

    /// This function computes the integral of the function f over the infinite interval (-\infty,+\infty). The integral is mapped onto the
    /// semi-open interval (0,1] using the transformation x = (1-t)/t,
    /// 
    /// \int_{-\infty}^{+\infty} dx f(x) = 
    ///      \int_0^1 dt (f((1-t)/t) + f((-1+t)/t))/t^2.
    /// 
    /// It is then integrated using the QAGS algorithm. The normal 21-point Gauss-Kronrod rule of QAGS is replaced by a 15-point rule, because
    /// the transformation can generate an integrable singularity at the origin. In this case a lower-order rule is more efficient.
    pub fn qagi<T>(&self, f: ::function<T>, arg: &mut T, epsabs: f64, epsrel: f64, limit: u64, result: &mut f64, abserr: &mut f64) -> enums::Value {
        let mut s = InternParam{func: f, param: arg, p2: 0f64};

        unsafe { intern_qags(i_transform, &mut s, 0f64, 1f64, epsabs, epsrel, limit, self, result, abserr, ::integration::qk15) }
    }

    /// This function computes the integral of the function f over the semi-infinite interval (a,+\infty). The integral is mapped onto the
    /// semi-open interval (0,1] using the transformation x = a + (1-t)/t,
    /// 
    /// \int_{a}^{+\infty} dx f(x) = 
    ///      \int_0^1 dt f(a + (1-t)/t)/t^2
    /// 
    /// and then integrated using the QAGS algorithm.
    pub fn qagiu<T>(&self, f: ::function<T>, arg: &mut T, a: f64, epsabs: f64, epsrel: f64, limit: u64, result: &mut f64, abserr: &mut f64) -> enums::Value {
        let mut s = InternParam{func: f, param: arg, p2: a};

        unsafe { intern_qags(iu_transform, &mut s, 0f64, 1f64, epsabs, epsrel, limit, self, result, abserr, ::integration::qk15) }
    }

    /// This function computes the integral of the function f over the semi-infinite interval (-\infty,b). The integral is mapped onto the semi-open
    /// interval (0,1] using the transformation x = b - (1-t)/t,
    /// 
    /// \int_{-\infty}^{b} dx f(x) = 
    ///      \int_0^1 dt f(b - (1-t)/t)/t^2
    /// 
    /// and then integrated using the QAGS algorithm.
    pub fn qagil<T>(&self, f: ::function<T>, arg: &mut T, b: f64, epsabs: f64, epsrel: f64, limit: u64, result: &mut f64, abserr: &mut f64) -> enums::Value {
        let mut s = InternParam{func: f, param: arg, p2: b};

        unsafe { intern_qags(il_transform, &mut s, 0f64, 1f64, epsabs, epsrel, limit, self, result, abserr, ::integration::qk15) }
    }

    /// This function computes the Cauchy principal value of the integral of f over (a,b), with a singularity at c,
    /// 
    /// I = \int_a^b dx f(x) / (x - c)
    /// 
    /// The adaptive bisection algorithm of QAG is used, with modifications to ensure that subdivisions do not occur at the singular point x = c.
    /// When a subinterval contains the point x = c or is close to it then a special 25-point modified Clenshaw-Curtis rule is used to control
    /// the singularity. Further away from the singularity the algorithm uses an ordinary 15-point Gauss-Kronrod integration rule.
    #[allow(dead_assignment)]
    pub fn qawc<T>(&self, f: ::function<T>, arg: &mut T, a: f64, b: f64, c: f64, epsabs: f64, epsrel: f64, limit: u64,
        result: &mut f64, abserr: &mut f64) -> enums::Value {
        let mut result0 = 0f64;
        let mut abserr0 = 0f64;
        let mut roundoff_type1 = 0i32;
        let mut roundoff_type2 = 0i32;
        let mut error_type = 0i32;
        let mut err_reliable = 0i32;
        let mut sign = 1f64;
        let mut lower = 0f64;
        let mut higher = 0f64;

        /* Initialize results */
        *result = 0f64;
        *abserr = 0f64;

        unsafe {
            if limit > (*self.w).limit {
                rgsl_error!("iteration limit exceeds available workspace", enums::Inval);
            }

            if b < a {
                lower = b;
                higher = a;
                sign = -1f64;
            } else {
                lower = a;
                higher = b;
            }

            self.initialise(lower, higher);

            if epsabs <= 0f64 && (epsrel < 50f64 * ::DBL_EPSILON || epsrel < 0.5e-28f64) {
                rgsl_error!("tolerance cannot be acheived with given epsabs and epsrel", enums::BadTol);
            }

            if c == a || c == b {
                rgsl_error!("cannot integrate with singularity on endpoint", enums::Inval);
            }      

            /* perform the first integration */
            qc25c(f, arg, lower, higher, c, &mut result0, &mut abserr0, &mut err_reliable);

            self.set_initial_result(result0, abserr0);

            /* Test on accuracy, use 0.01 relative error as an extra safety margin on the first iteration (ignored for subsequent iterations) */
            let mut tolerance = epsabs.max(epsrel * fabsf64(result0));

            if abserr0 < tolerance && abserr0 < 0.01f64 * fabsf64(result0) {
                *result = sign * result0;
                *abserr = abserr0;

                return enums::Success;
            } else if limit == 1 {
                *result = sign * result0;
                *abserr = abserr0;

                rgsl_error!("a maximum of one iteration was insufficient", enums::MaxIter);
            }

            let mut area = result0;
            let mut errsum = abserr0;

            let mut iteration = 1;

            loop {
                let mut a_i = 0f64;
                let mut b_i = 0f64;
                let mut r_i = 0f64;
                let mut e_i = 0f64;
                let mut area1 = 0f64;
                let mut area2 = 0f64;
                let mut error1 = 0f64;
                let mut error2 = 0f64;
                let mut err_reliable1 = 0i32;
                let mut err_reliable2 = 0i32;

                /* Bisect the subinterval with the largest error estimate */
                self.retrieve(&mut a_i, &mut b_i, &mut r_i, &mut e_i);

                let a1 = a_i; 
                let mut b1 = 0.5 * (a_i + b_i);
                let mut a2 = b1;
                let b2 = b_i;

                if c > a1 && c <= b1 {
                    b1 = 0.5f64 * (c + b2);
                    a2 = b1;
                } else if c > b1 && c < b2 {
                    b1 = 0.5f64 * (a1 + c);
                    a2 = b1;
                }

                qc25c(f, arg, a1, b1, c, &mut area1, &mut error1, &mut err_reliable1);
                qc25c(f, arg, a2, b2, c, &mut area2, &mut error2, &mut err_reliable2);

                let area12 = area1 + area2;
                let error12 = error1 + error2;

                errsum += error12 - e_i;
                area += area12 - r_i;

                if err_reliable1 != 0 && err_reliable2 != 0 {
                    let delta = r_i - area12;

                    if fabsf64(delta) <= 1.0e-5f64 * fabsf64(area12) && error12 >= 0.99f64 * e_i {
                        roundoff_type1 += 1;
                    }
                    if iteration >= 10 && error12 > e_i {
                        roundoff_type2 += 1;
                    }
                }

                tolerance = epsabs.max(epsrel * fabsf64(area));

                if errsum > tolerance {
                    if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
                        /* round off error */
                        error_type = 2;
                    }

                    /* set error flag in the case of bad integrand behaviour at a point of the integration range */
                    if ::util::subinterval_too_small(a1, a2, b2) {
                        error_type = 3;
                    }
                }

                self.update(a1, b1, area1, error1, a2, b2, area2, error2);

                self.retrieve(&mut a_i, &mut b_i, &mut r_i, &mut e_i);

                iteration += 1;

                if iteration < limit && error_type == 0 && errsum > tolerance {
                } else {
                    break;
                }
            }

            *result = sign * self.sum_results();
            *abserr = errsum;

            if errsum <= tolerance {
              enums::Success
            } else if error_type == 2 {
                rgsl_error!("roundoff error prevents tolerance from being achieved", enums::Round);
                enums::Round
            } else if error_type == 3 {
                rgsl_error!("bad integrand behavior found in the integration interval", enums::Sing);
                enums::Sing
            } else if iteration == limit {
                rgsl_error!("maximum number of subdivisions reached", enums::MaxIter);
                enums::MaxIter
            } else {
                rgsl_error!("could not integrate function", enums::Failed);
                enums::Failed
            }
        }
    }

    pub fn sort_results(&self) {
        unsafe {
            let nint = (*self.w).size as uint;
            let mut t_elist = CVec::new((*self.w).elist, nint);
            let mut t_order = CVec::new((*self.w).order, nint);
            let elist = t_elist.as_mut_slice();
            let order = t_order.as_mut_slice();

            for i in range(0u, nint) {
                let i1 = order[i] as uint;
                let mut e1 = elist[i1];
                let mut i_max = i1;

                for j in range(i + 1, nint) {
                    let i2 = order[j] as uint;
                    let e2 = elist[i2];

                    if e2 >= e1 {
                        i_max = i2;
                        e1 = e2;
                    }
                }

                if i_max != i1 {
                    order[i] = order[i_max];
                    order[i_max] = i1 as u64;
                }
            }
            (*self.w).i = order[0];
        }
    }

    pub fn qpsrt(&self) {
        let w = self.w;

        unsafe {
            let mut order = CVec::new((*w).order, (*w).size as uint + 1u);
            let last = (*w).size - 1;
            let limit = (*w).limit;

            let mut i_nrmax = (*w).nrmax;
            let mut i_maxerr = order.as_slice()[i_nrmax as uint];

            // Check whether the list contains more than two error estimates
            if last < 2 {
                order.as_mut_slice()[0u] = 0;
                order.as_mut_slice()[1u] = 1;
                (*w).i = i_maxerr;
                return ;
            }

            // Compute the number of elements in the list to be maintained in descending order. This number depends on the number of
            // subdivisions still allowed.
            let top =  if last < (limit / 2 + 2) {
                last
            } else {
                limit - last + 1
            };

            let elist = CVec::new((*w).elist, top as uint + 1u);
            let errmax = elist.as_slice()[i_maxerr as uint];

            // This part of the routine is only executed if, due to a difficult integrand, subdivision increased the error estimate. In the normal
            // case the insert procedure should start after the nrmax-th largest error estimate.
            while i_nrmax > 0 && errmax > elist.as_slice()[order.as_slice()[i_nrmax as uint - 1u] as uint] {
                order.as_mut_slice()[i_nrmax as uint] = order.as_slice()[i_nrmax as uint - 1u];
                i_nrmax -= 1;
            }

            // Insert errmax by traversing the list top-down, starting comparison from the element elist(order(i_nrmax+1)).
            let mut i = i_nrmax + 1;

            // The order of the tests in the following line is important to prevent a segmentation fault
            while i < top && errmax < elist.as_slice()[order.as_slice()[i as uint] as uint] {
                order.as_mut_slice()[i as uint - 1u] = order.as_slice()[i as uint];
                i += 1;
            }

            order.as_mut_slice()[i as uint - 1u] = i_maxerr;

            // Insert errmin by traversing the list bottom-up
            let errmin = elist.as_slice()[last as uint];

            let mut k = top - 1;
            while k > i - 2 && errmin >= elist.as_slice()[order.as_slice()[k as uint] as uint] {
                order.as_mut_slice()[k as uint + 1u] = order.as_slice()[k as uint];
                k -= 1;
            }
      
            order.as_mut_slice()[k as uint + 1u] = last;

            // Set i_max and e_max
            i_maxerr = order.as_slice()[i_nrmax as uint] ;

            (*w).i = i_maxerr ;
            (*w).nrmax = i_nrmax ;
        }
    }

    pub fn sum_results(&self) -> f64 {
        unsafe {
            let f_w = self.w;
            let mut result_sum = 0f64;
            let v_rlist = CVec::new((*f_w).rlist, (*f_w).size as uint);
            let rlist = v_rlist.as_slice();

            for k in range(0u, (*f_w).size as uint) {
                result_sum += rlist[k];
            }

            result_sum
        }
    }

    pub fn retrieve(&self, a: &mut f64, b: &mut f64, r: &mut f64, e: &mut f64) {
        unsafe {
            let w = self.w;
            let alist = CVec::new((*w).alist, (*w).i as uint + 1u);
            let blist = CVec::new((*w).blist, (*w).i as uint + 1u);
            let rlist = CVec::new((*w).rlist, (*w).i as uint + 1u);
            let elist = CVec::new((*w).elist, (*w).i as uint + 1u);

            let i = (*w).i as uint;

            *a = alist.as_slice()[i];
            *b = blist.as_slice()[i];
            *r = rlist.as_slice()[i];
            *e = elist.as_slice()[i];
        }
    }

    pub fn update(&self, a1: f64, b1: f64, area1: f64, error1: f64, a2: f64, b2: f64, area2: f64, error2: f64) {
        let w = self.w;

        unsafe {
            let i_max = (*w).i as uint;
            let i_new = (*w).size as uint;
            let tmp = if i_max > i_new {
                i_max + 1u
            } else {
                i_new + 1u
            };
            let mut alist = CVec::new((*w).alist, tmp);
            let mut blist = CVec::new((*w).blist, tmp);
            let mut rlist = CVec::new((*w).rlist, tmp);
            let mut elist = CVec::new((*w).elist, tmp);
            let mut level = CVec::new((*w).level, tmp);

            let new_level = level.as_slice()[i_max] + 1;

            /* append the newly-created intervals to the list */

            if error2 > error1 {
                // blist[maxerr] is already == b2
                alist.as_mut_slice()[i_max] = a2;
                rlist.as_mut_slice()[i_max] = area2;
                elist.as_mut_slice()[i_max] = error2;
                level.as_mut_slice()[i_max] = new_level;

                alist.as_mut_slice()[i_new] = a1;
                blist.as_mut_slice()[i_new] = b1;
                rlist.as_mut_slice()[i_new] = area1;
                elist.as_mut_slice()[i_new] = error1;
                level.as_mut_slice()[i_new] = new_level;
            } else {
                // alist[maxerr] is already == a1
                blist.as_mut_slice()[i_max] = b1;
                rlist.as_mut_slice()[i_max] = area1;
                elist.as_mut_slice()[i_max] = error1;
                level.as_mut_slice()[i_max] = new_level;

                alist.as_mut_slice()[i_new] = a2;
                blist.as_mut_slice()[i_new] = b2;
                rlist.as_mut_slice()[i_new] = area2;
                elist.as_mut_slice()[i_new] = error2;
                level.as_mut_slice()[i_new] = new_level;
            }

            (*w).size += 1;

            if new_level > (*w).maximum_level {
                (*w).maximum_level = new_level;
            }

            self.qpsrt();
        }
    }

    pub fn set_initial_result(&self, result: f64, error: f64) {
        unsafe {
            let mut rlist = CVec::new((*self.w).rlist, 1);
            let mut elist = CVec::new((*self.w).elist, 1);

            (*self.w).size = 1;
            rlist.as_mut_slice()[0] = result;
            elist.as_mut_slice()[0] = error;
        }
    }

    pub fn initialise(&self, a: f64, b: f64) {
        let w = self.w;

        unsafe {
            let mut alist = CVec::new((*w).alist, 1);
            let mut blist = CVec::new((*w).blist, 1);
            let mut rlist = CVec::new((*w).rlist, 1);
            let mut elist = CVec::new((*w).elist, 1);
            let mut order = CVec::new((*w).order, 1);
            let mut level = CVec::new((*w).level, 1);

            (*w).size = 0;
            (*w).nrmax = 0;
            (*w).i = 0;
            alist.as_mut_slice()[0] = a;
            blist.as_mut_slice()[0] = b;
            rlist.as_mut_slice()[0] = 0f64;
            elist.as_mut_slice()[0] = 0f64;
            order.as_mut_slice()[0] = 0u64;
            level.as_mut_slice()[0] = 0u64;

            (*w).maximum_level = 0;
        }
    }

    fn append_interval(&self, a1: f64, b1: f64, area1: f64, error1: f64) {
        unsafe {
            let w = self.w;
            let i_new = (*w).size as uint;
            let mut alist = CVec::new((*w).alist, i_new + 1u);
            let mut blist = CVec::new((*w).blist, i_new + 1u);
            let mut rlist = CVec::new((*w).rlist, i_new + 1u);
            let mut elist = CVec::new((*w).elist, i_new + 1u);
            let mut order = CVec::new((*w).order, i_new + 1u);
            let mut level = CVec::new((*w).level, i_new + 1u);

            alist.as_mut_slice()[i_new] = a1;
            blist.as_mut_slice()[i_new] = b1;
            rlist.as_mut_slice()[i_new] = area1;
            elist.as_mut_slice()[i_new] = error1;
            order.as_mut_slice()[i_new] = i_new as u64;
            level.as_mut_slice()[i_new] = 0;

            (*w).size += 1;
        }
    }

    pub fn limit(&self) -> u64 {
        unsafe { (*self.w).limit }
    }
}

impl Drop for IntegrationWorkspace {
    fn drop(&mut self) {
        unsafe { ffi::gsl_integration_workspace_free(self.w) };
        self.w = ::std::ptr::mut_null();
    }
}

impl ffi::FFI<ffi::gsl_integration_workspace> for IntegrationWorkspace {
    fn wrap(w: *mut ffi::gsl_integration_workspace) -> IntegrationWorkspace {
        IntegrationWorkspace {
            w: w
        }
    }

    fn unwrap(w: &IntegrationWorkspace) -> *mut ffi::gsl_integration_workspace {
        w.w
    }
}

/// The QAWS algorithm is designed for integrands with algebraic-logarithmic singularities at the end-points of an integration region. In order
/// to work efficiently the algorithm requires a precomputed table of Chebyshev moments.
pub struct IntegrationQawsTable {
    w: *mut ffi::gsl_integration_qaws_table
}

impl IntegrationQawsTable {
    /// This function allocates space for a gsl_integration_qaws_table struct describing a singular weight function W(x) with the parameters
    /// (\alpha, \beta, \mu, \nu),
    /// 
    /// W(x) = (x-a)^alpha (b-x)^beta log^mu (x-a) log^nu (b-x)
    /// 
    /// where \alpha > -1, \beta > -1, and \mu = 0, 1, \nu = 0, 1. The weight function can take four different forms depending on the values
    /// of \mu and \nu,
    /// 
    /// W(x) = (x-a)^alpha (b-x)^beta                   (mu = 0, nu = 0)
    /// W(x) = (x-a)^alpha (b-x)^beta log(x-a)          (mu = 1, nu = 0)
    /// W(x) = (x-a)^alpha (b-x)^beta log(b-x)          (mu = 0, nu = 1)
    /// W(x) = (x-a)^alpha (b-x)^beta log(x-a) log(b-x) (mu = 1, nu = 1)
    /// 
    /// The singular points (a,b) do not have to be specified until the integral is computed, where they are the endpoints of the integration
    /// range.
    /// 
    /// The function returns a pointer to the newly allocated table gsl_integration_qaws_table if no errors were detected, and 0 in the case
    /// of error.
    pub fn new(alpha: f64, beta: f64, mu: i32, nu: i32) -> Option<IntegrationQawsTable> {
        let tmp = unsafe { ffi::gsl_integration_qaws_table_alloc(alpha, beta, mu, nu) };

        if tmp.is_null() {
            None
        } else {
            Some(IntegrationQawsTable {
                w: tmp
            })
        }
    }

    /// This function modifies the parameters (\alpha, \beta, \mu, \nu)
    pub fn set(&self, alpha: f64, beta: f64, mu: i32, nu: i32) -> enums::Value {
        unsafe { ffi::gsl_integration_qaws_table_set(self.w, alpha, beta, mu, nu) }
    }

    /// This function computes the integral of the function f(x) over the interval (a,b) with the singular weight function (x-a)^\alpha
    /// (b-x)^\beta \log^\mu (x-a) \log^\nu (b-x). The parameters of the weight function (\alpha, \beta, \mu, \nu) are taken from the
    /// table self. The integral is,
    /// 
    /// I = \int_a^b dx f(x) (x-a)^alpha (b-x)^beta log^mu (x-a) log^nu (b-x).
    /// 
    /// The adaptive bisection algorithm of QAG is used. When a subinterval contains one of the endpoints then a special 25-point modified
    /// Clenshaw-Curtis rule is used to control the singularities. For subintervals which do not include the endpoints an ordinary 15-point
    /// Gauss-Kronrod integration rule is used.
    #[allow(dead_assignment)]
    pub fn qaws<T>(&self, f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, workspace: &IntegrationWorkspace,
        result: &mut f64, abserr: &mut f64) -> enums::Value {
        let mut result0 = 0f64;
        let mut abserr0 = 0f64;
        let mut roundoff_type1 = 0i32;
        let mut roundoff_type2 = 0i32;
        let mut error_type = 0i32;

        /* Initialize results */
        workspace.initialise(a, b);

        *result = 0f64;
        *abserr = 0f64;

        unsafe {
            if limit > (*workspace.w).limit {
                rgsl_error!("iteration limit exceeds available workspace", enums::Inval);
            }

            if b <= a {
                rgsl_error!("limits must form an ascending sequence, a < b", enums::Inval);
            }

            if epsabs <= 0f64 && (epsrel < 50f64 * ::DBL_EPSILON || epsrel < 0.5e-28f64) {
                rgsl_error!("tolerance cannot be acheived with given epsabs and epsrel", enums::BadTol);
            }

            /* perform the first integration */
            {
                let mut area1 = 0f64;
                let mut area2 = 0f64;
                let mut error1 = 0f64;
                let mut error2 = 0f64;
                let mut err_reliable1 = false;
                let mut err_reliable2 = false;
                let a1 = a;
                let b1 = 0.5f64 * (a + b);
                let a2 = b1;
                let b2 = b;

                qc25s(f, arg, a, b, a1, b1, self.w, &mut area1, &mut error1, &mut err_reliable1);
                qc25s(f, arg, a, b, a2, b2, self.w, &mut area2, &mut error2, &mut err_reliable2);

                if error1 > error2 {
                    workspace.append_interval(a1, b1, area1, error1);
                    workspace.append_interval(a2, b2, area2, error2);
                } else {
                    workspace.append_interval(a2, b2, area2, error2);
                    workspace.append_interval(a1, b1, area1, error1);
                }

                result0 = area1 + area2;
                abserr0 = error1 + error2;
            }

            /* Test on accuracy */
            let mut tolerance = epsabs.max(epsrel * fabsf64(result0));

            // Test on accuracy, use 0.01 relative error as an extra safety margin on the first iteration (ignored for subsequent iterations)
            if abserr0 < tolerance && abserr0 < 0.01 * fabsf64(result0) {
                *result = result0;
                *abserr = abserr0;

                return enums::Success;
            } else if limit == 1 {
                *result = result0;
                *abserr = abserr0;

                rgsl_error!("a maximum of one iteration was insufficient", enums::MaxIter);
            }

            let mut area = result0;
            let mut errsum = abserr0;

            let mut iteration = 2;

            loop {
                let mut a_i = 0f64;
                let mut b_i = 0f64;
                let mut r_i = 0f64;
                let mut e_i = 0f64;
                let mut area1 = 0f64;
                let mut area2 = 0f64;
                let mut error1 = 0f64;
                let mut error2 = 0f64;
                let mut err_reliable1 = false;
                let mut err_reliable2 = false;

                /* Bisect the subinterval with the largest error estimate */
                workspace.retrieve(&mut a_i, &mut b_i, &mut r_i, &mut e_i);

                let a1 = a_i; 
                let b1 = 0.5f64 * (a_i + b_i);
                let a2 = b1;
                let b2 = b_i;

                qc25s(f, arg, a, b, a1, b1, self.w, &mut area1, &mut error1, &mut err_reliable1);
                qc25s(f, arg, a, b, a2, b2, self.w, &mut area2, &mut error2, &mut err_reliable2);

                let area12 = area1 + area2;
                let error12 = error1 + error2;

                errsum += error12 - e_i;
                area += area12 - r_i;

                if err_reliable1 && err_reliable2 {
                    let delta = r_i - area12;

                    if fabsf64(delta) <= 1.0e-5f64 * fabsf64(area12) && error12 >= 0.99 * e_i {
                        roundoff_type1 += 1;
                    }
                    if iteration >= 10 && error12 > e_i {
                        roundoff_type2 += 1;
                    }
                }

                tolerance = epsabs.max(epsrel * fabsf64(area));

                if errsum > tolerance {
                    if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
                        /* round off error */
                        error_type = 2;
                    }

                    /* set error flag in the case of bad integrand behaviour at a point of the integration range */
                    if ::util::subinterval_too_small(a1, a2, b2) {
                        error_type = 3;
                    }
                }

                workspace.update(a1, b1, area1, error1, a2, b2, area2, error2);
                workspace.retrieve(&mut a_i, &mut b_i, &mut r_i, &mut e_i);

                iteration += 1;

                if iteration < limit && error_type == 0 && errsum > tolerance {
                } else {
                    break;
                }
            }

            *result = workspace.sum_results();
            *abserr = errsum;

            if errsum <= tolerance {
                enums::Success
            } else if error_type == 2 {
                rgsl_error!("roundoff error prevents tolerance from being achieved", enums::Round);
                enums::Round
            } else if error_type == 3 {
                rgsl_error!("bad integrand behavior found in the integration interval", enums::Sing);
                enums::Sing
            } else if iteration == limit {
                rgsl_error!("maximum number of subdivisions reached", enums::MaxIter);
                enums::MaxIter
            } else {
                rgsl_error!("could not integrate function", enums::Failed);
                enums::Failed
            }
        }
    }
}

impl Drop for IntegrationQawsTable {
    fn drop(&mut self) {
        unsafe { ffi::gsl_integration_qaws_table_free(self.w) };
        self.w = ::std::ptr::mut_null();
    }
}

impl ffi::FFI<ffi::gsl_integration_qaws_table> for IntegrationQawsTable {
    fn wrap(w: *mut ffi::gsl_integration_qaws_table) -> IntegrationQawsTable {
        IntegrationQawsTable {
            w: w
        }
    }

    fn unwrap(w: &IntegrationQawsTable) -> *mut ffi::gsl_integration_qaws_table {
        w.w
    }
}

fn i_transform<T>(t: f64, params: &mut InternParam<T>) -> f64 {
    let f = params.func;
    let x = (1f64 - t) / t;
    let y = f(x, params.param) + f(-x, params.param);
  
    (y / t) / t
}

fn iu_transform<T>(t: f64, p: &mut InternParam<T>) -> f64 {
    let a = p.p2;
    let f = p.func;
    let x = a + (1f64 - t) / t;
    let y = f(x, p.param);

    (y / t) / t
}

fn il_transform<T>(t: f64, p: &mut InternParam<T>) -> f64 {
    let b = p.p2;
    let f = p.func;
    let x = b - (1f64 - t) / t;
    let y = f(x, p.param);

    (y / t) / t
}

fn intern_qag<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, f_w: &IntegrationWorkspace,
    result: &mut f64, abserr: &mut f64, q: ::integration_function<T>) -> enums::Value {
    let w = f_w.w;
    let mut roundoff_type1 = 0i32;
    let mut roundoff_type2 = 0i32;
    let mut error_type = 0i32;

    *result = 0f64;
    *abserr = 0f64;

    // Initialize results
    f_w.initialise(a, b);

    if unsafe { limit > (*w).limit } {
        rgsl_error!("iteration limit exceeds available workspace", enums::Inval);
    }
    if epsabs <= 0f64 && (epsrel < 50f64 * ::DBL_EPSILON || epsrel < 0.5e-28f64) {
        rgsl_error!("tolerance cannot be achieved with given epsabs and epsrel", enums::BadTol);
    }

    // perform the first integration
    let mut result0 = 0f64;
    let mut abserr0 = 0f64;
    let mut resabs0 = 0f64;
    let mut resasc0 = 0f64;
    q(f, arg, a, b, &mut result0, &mut abserr0, &mut resabs0, &mut resasc0);

    f_w.set_initial_result(result0, abserr0);

    // Test on accuracy
    let mut tolerance = unsafe { epsabs.max(epsrel * fabsf64(result0)) };

    // need IEEE rounding here to match original quadpack behavior
    let round_off = 50f64 * ::DBL_EPSILON * resabs0;

    if abserr0 <= round_off && abserr0 > tolerance {
        *result = result0;
        *abserr = abserr0;

        rgsl_error!("cannot reach tolerance because of roundoff error on first attempt", enums::Round);
    } else if (abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0f64 {
        *result = result0;
        *abserr = abserr0;

        return enums::Success;
    } else if limit == 1 {
        *result = result0;
        *abserr = abserr0;

        rgsl_error!("a maximum of one iteration was insufficient", enums::MaxIter);
    }

    let mut area = result0;
    let mut errsum = abserr0;

    let iteration = 1u64;

    loop {
        let mut area1 = 0f64;
        let mut area2 = 0f64;
        let mut error1 = 0f64;
        let mut error2 = 0f64;
        let mut a_i = 0f64;
        let mut b_i = 0f64;
        let mut r_i = 0f64;
        let mut e_i = 0f64;
        let mut resabs1 = 0f64;
        let mut resasc1 = 0f64;
        let mut resabs2 = 0f64;
        let mut resasc2 = 0f64;

        // Bisect the subinterval with the largest error estimate
        f_w.retrieve(&mut a_i, &mut b_i, &mut r_i, &mut e_i);

        let a1 = a_i; 
        let b1 = 0.5 * (a_i + b_i);
        let a2 = b1;
        let b2 = b_i;

        q(f, arg, a1, b1, &mut area1, &mut error1, &mut resabs1, &mut resasc1);
        q(f, arg, a2, b2, &mut area2, &mut error2, &mut resabs2, &mut resasc2);

        let area12 = area1 + area2;
        let error12 = error1 + error2;

        errsum += error12 - e_i;
        area += area12 - r_i;

        if resasc1 != error1 && resasc2 != error2 {
            let delta = r_i - area12;

            if unsafe { fabsf64(delta) <= 1.0e-5f64 * fabsf64(area12) && error12 >= 0.99f64 * e_i } {
                roundoff_type1 += 1;
            }
            if iteration >= 10 && error12 > e_i {
                roundoff_type2 += 1;
            }
        }

        tolerance = unsafe { epsabs.max(epsrel * fabsf64(area)) };

        if errsum > tolerance {
            if roundoff_type1 >= 6 || roundoff_type2 >= 20 {
                error_type = 2;   /* round off error */
            }

            // set error flag in the case of bad integrand behaviour at
            // a point of the integration range

            if ::util::subinterval_too_small(a1, a2, b2) {
                error_type = 3;
            }
        }

        if iteration < limit && error_type == 0i32 && errsum > tolerance {
        } else {
            break;
        }
    }
    *result = f_w.sum_results();
    *abserr = errsum;

    if errsum <= tolerance {
        enums::Success
    } else if error_type == 2 {
        rgsl_error!("roundoff error prevents tolerance from being achieved", enums::Round);
        enums::Round
    } else if error_type == 3 {
        rgsl_error!("bad integrand behavior found in the integration interval", enums::Sing);
        enums::Sing
    } else if iteration == limit {
        rgsl_error!("maximum number of subdivisions reached", enums::MaxIter);
        enums::MaxIter
    } else {
        rgsl_error!("could not integrate function", enums::Failed);
        enums::Failed
    }
}

unsafe fn initialise_table(table: *mut ffi::extrapolation_table) {
    (*table).n = 0;
    (*table).nres = 0;
}

unsafe fn append_table(table: &mut ffi::extrapolation_table, y: f64) {
  let n = (*table).n as uint;

  (*table).rlist2[n] = y;
  (*table).n += 1;
}

unsafe fn intern_qelg(table: &mut ffi::extrapolation_table, result: &mut f64, abserr: &mut f64) {
    let mut epstab = (*table).rlist2;//CVec::new((*table).rlist2 as *mut f64, (*table).n as uint + 3);
    let mut res3la = (*table).res3la;//CVec::new((*table).res3la as *mut f64, 3u);
    let n = (*table).n as uint - 1u;

    let current = epstab[n];

    let mut absolute = ::DBL_MAX;
    let mut relative = 5f64 * ::DBL_EPSILON * fabsf64(current);

    let newelm = n / 2u;
    let n_orig = n;
    let mut n_final = n;

    let nres_orig = (*table).nres;

    *result = current;
    *abserr = ::DBL_MAX;

    if n < 2 {
        *result = current;
        *abserr = absolute.max(relative);
        return;
    }

    epstab[n + 2] = epstab[n];
    epstab[n] = ::DBL_MAX;

    for i in range(0, newelm) {
        let mut res = epstab[n - 2 * i + 2];
        let e0 = epstab[n - 2 * i - 2];
        let e1 = epstab[n - 2 * i - 1];
        let e2 = res;

        let e1abs = fabsf64(e1);
        let delta2 = e2 - e1;
        let err2 = fabsf64(delta2);
        let tol2 = fabsf64(e2).max(e1abs) * ::DBL_EPSILON;
        let delta3 = e1 - e0;
        let err3 = fabsf64(delta3);
        let tol3 = e1abs.max(fabsf64(e0)) * ::DBL_EPSILON;

        if err2 <= tol2 && err3 <= tol3 {
            /* If e0, e1 and e2 are equal to within machine accuracy, convergence is assumed.  */
            *result = res;
            absolute = err2 + err3;
            relative = 5f64 * ::DBL_EPSILON * fabsf64(res);
            *abserr = absolute.max(relative);
            return;
        }

        let e3 = epstab[n - 2 * i];
        epstab[n - 2 * i] = e1;
        let delta1 = e1 - e3;
        let err1 = fabsf64(delta1);
        let tol1 = e1abs.max(fabsf64(e3)) * ::DBL_EPSILON;

        /* If two elements are very close to each other, omit a part of the table by adjusting the value of n */
        if err1 <= tol1 || err2 <= tol2 || err3 <= tol3 {
            n_final = 2 * i;
            break;
        }

        let ss = (1f64 / delta1 + 1f64 / delta2) - 1f64 / delta3;

        /* Test to detect irregular behaviour in the table, and eventually omit a part of the table by adjusting the value of n. */

        if fabsf64(ss * e1) <= 0.0001f64 {
            n_final = 2 * i;
            break;
        }

        /* Compute a new element and eventually adjust the value of result. */

        res = e1 + 1f64 / ss;
        epstab[n - 2 * i] = res;

        {
            let error = err2 + fabsf64(res - e2) + err3;

            if error <= *abserr {
                *abserr = error;
                *result = res;
            }
        }
    }

    /* Shift the table */
    {
        let limexp = 49u64;

        if n_final == limexp as uint {
            n_final = 2u * (limexp as uint / 2u);
        }
    }

    if n_orig & 1 == 1 {
        for i in range(0, newelm + 1) {
          epstab[1 + i * 2] = epstab[i * 2 + 3];
        }
    } else {
        for i in range(0, newelm + 1) {
            epstab[i * 2] = epstab[i * 2 + 2];
        }
    }

    if n_orig != n_final {
        for i in range(0, n_final + 1) {
            epstab[i] = epstab[n_orig - n_final + i];
        }
    }

    (*table).n = n_final as u64 + 1;

    if nres_orig < 3 {
        res3la.as_mut_slice()[nres_orig as uint] = *result;
        *abserr = ::DBL_MAX;
    } else {
        /* Compute error estimate */
        *abserr = fabsf64(*result - res3la[2]) + fabsf64(*result - res3la[1]) + fabsf64(*result - res3la[0]);

        res3la[0] = res3la[1];
        res3la[1] = res3la[2];
        res3la[2] = *result;
    }

    /* In QUADPACK the variable table->nres is incremented at the top of qelg, so it increases on every call. This leads to the array
       res3la being accessed when its elements are still undefined, so I have moved the update to this point so that its value more
       useful. */

    (*table).nres = nres_orig + 1;  

    *abserr = (*abserr).max(5f64 * ::DBL_EPSILON * fabsf64(*result));
}

unsafe fn test_positivity(result: f64, resabs: f64) -> bool {
    (fabsf64(result) >= (1f64 - 50f64 * ::DBL_EPSILON) * resabs)
}

unsafe fn increase_nrmax(workspace: *mut ffi::gsl_integration_workspace) -> bool {
    let id = (*workspace).nrmax;

    let t_order = CVec::new((*workspace).order, (*workspace).nrmax as uint + 1u);
    let order = t_order.as_slice();
    let t_level = CVec::new((*workspace).level, order[(*workspace).nrmax as uint] as uint + 1u);
    let level = t_level.as_slice();

    let limit = (*workspace).limit;
    let last = (*workspace).size - 1;

    let jupbnd = if last > (1 + limit / 2) {
        limit + 1 - last
    } else {
        last
    };
  
    for k in range(id, jupbnd + 1) {
        let i_max = order[(*workspace).nrmax as uint];
      
        (*workspace).i = i_max ;
        if level[i_max as uint] < (*workspace).maximum_level {
            return true;
        }
        (*workspace).nrmax += 1;
    }
    false
}

unsafe fn large_interval(workspace: *mut ffi::gsl_integration_workspace) -> bool {
    let i = (*workspace).i ;
    let level = CVec::new((*workspace).level, i as uint + 1u);
  
    if level.as_slice()[i as uint] < (*workspace).maximum_level {
        true
    } else {
        false
    }
}

unsafe fn reset_nrmax(workspace: *mut ffi::gsl_integration_workspace) {
    (*workspace).nrmax = 0;
    (*workspace).i = *(*workspace).order;
}

unsafe fn compute_result(w: &IntegrationWorkspace, result: &mut f64, abserr: &mut f64, errsum: f64,
    error_type: i32) -> enums::Value {
    *result = w.sum_results();
    *abserr = errsum;
    return_error(error_type)
}

unsafe fn return_error(t_error_type: i32) -> enums::Value {
    let error_type = if t_error_type > 2 {
        t_error_type - 1
    } else {
        t_error_type
    };

    match error_type {
        0 => enums::Success,
        1 => {
            rgsl_error!("number of iterations was insufficient", enums::MaxIter);
            enums::MaxIter
        }
        2 => {
            rgsl_error!("cannot reach tolerance because of roundoff error", enums::Round);
            enums::Round
        }
        3 => {
            rgsl_error!("bad integrand behavior found in the integration interval", enums::Sing);
            enums::Sing
        }
        4 => {
            rgsl_error!("roundoff error detected in the extrapolation table", enums::Round);
            enums::Round
        }
        5 => {
            rgsl_error!("integral is divergent, or slowly convergent", enums::Round);
            enums::Round
        }
        _ => {
            rgsl_error!("could not integrate function", enums::Failed);
            enums::Failed
        }
    }
}

unsafe fn intern_qags<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, f_w: &IntegrationWorkspace,
    result: &mut f64, abserr: &mut f64, q: ::integration_function<T>) -> enums::Value {
    let w = f_w.w;
    let mut ertest = 0f64;
    let mut error_over_large_intervals = 0f64;
    let mut reseps = 0f64;
    let mut abseps = 0f64;
    let mut correc = 0f64;
    let mut ktmin = 0u64;
    let mut roundoff_type1 = 0i32;
    let mut roundoff_type2 = 0i32;
    let mut roundoff_type3 = 0i32;
    let mut error_type = 0i32;
    let mut error_type2 = 0i32;
    let mut result0 = 0f64;
    let mut abserr0 = 0f64;
    let mut resabs0 = 0f64;
    let mut resasc0 = 0f64;

    let mut extrapolate = 0i32;
    let mut disallow_extrapolation = 0i32;

    let mut table : ffi::extrapolation_table = ::std::mem::zeroed();

    /* Initialize results */
    f_w.initialise(a, b);

    *result = 0f64;
    *abserr = 0f64;

    if limit > (*w).limit {
        rgsl_error!("iteration limit exceeds available workspace", enums::Inval);
    }

    /* Test on accuracy */
    if epsabs <= 0f64 && (epsrel < 50f64 * ::DBL_EPSILON || epsrel < 0.5e-28f64) {
        rgsl_error!("tolerance cannot be acheived with given epsabs and epsrel", enums::BadTol);
    }

    /* Perform the first integration */
    q(f, arg, a, b, &mut result0, &mut abserr0, &mut resabs0, &mut resasc0);

    f_w.set_initial_result(result0, abserr0);

    let mut tolerance = epsabs.max(epsrel * fabsf64(result0));

    if abserr0 <= 100f64 * ::DBL_EPSILON * resabs0 && abserr0 > tolerance {
        *result = result0;
        *abserr = abserr0;

        rgsl_error!("cannot reach tolerance because of roundoff error on first attempt", enums::Round);
    } else if (abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0f64 {
        *result = result0;
        *abserr = abserr0;

        return enums::Success;
    } else if limit == 1 {
        *result = result0;
        *abserr = abserr0;

        rgsl_error!("a maximum of one iteration was insufficient", enums::MaxIter);
    }

    /* Initialization */
    initialise_table(&mut table);
    append_table(&mut table, result0);

    let mut area = result0;
    let mut errsum = abserr0;

    let mut res_ext = result0;
    let mut err_ext = ::DBL_MAX;

    let positive_integrand = test_positivity(result0, resabs0);

    let mut iteration = 1u64;

    loop {
        let mut a_i = 0f64;
        let mut b_i = 0f64;
        let mut r_i = 0f64;
        let mut e_i = 0f64;
        let mut area1 = 0f64;
        let mut area2 = 0f64;
        let mut error1 = 0f64;
        let mut error2 = 0f64;
        let mut resasc1 = 0f64;
        let mut resasc2 = 0f64;
        let mut resabs1 = 0f64;
        let mut resabs2 = 0f64;

        /* Bisect the subinterval with the largest error estimate */
        f_w.retrieve(&mut a_i, &mut b_i, &mut r_i, &mut e_i);

        let t_level = CVec::new((*w).level, (*w).i as uint + 1);
        let current_level = t_level.as_slice()[(*w).i as uint] + 1;

        let a1 = a_i;
        let b1 = 0.5 * (a_i + b_i);
        let a2 = b1;
        let b2 = b_i;

        iteration += 1;

        q(f, arg, a1, b1, &mut area1, &mut error1, &mut resabs1, &mut resasc1);
        q(f, arg, a2, b2, &mut area2, &mut error2, &mut resabs2, &mut resasc2);

        let area12 = area1 + area2;
        let error12 = error1 + error2;
        let last_e_i = e_i;

        /* Improve previous approximations to the integral and test for accuracy.

        We write these expressions in the same way as the original
        QUADPACK code so that the rounding errors are the same, which
        makes testing easier. */

        errsum = errsum + error12 - e_i;
        area = area + area12 - r_i;

        tolerance = epsabs.max(epsrel * fabsf64(area));

        if resasc1 != error1 && resasc2 != error2 {
            let delta = r_i - area12;

            if fabsf64(delta) <= 1.0e-5f64 * fabsf64(area12) && error12 >= 0.99f64 * e_i {
                if extrapolate == 0 {
                  roundoff_type1 += 1;
                } else {
                  roundoff_type2 += 1;
                }
            }
            if iteration > 10 && error12 > e_i {
                roundoff_type3 += 1;
            }
        }

        /* Test for roundoff and eventually set error flag */
        if roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20 {
            /* round off error */
            error_type = 2;
        }

        if roundoff_type2 >= 5 {
            error_type2 = 1;
        }

        /* set error flag in the case of bad integrand behaviour at a point of the integration range */
        if ::util::subinterval_too_small(a1, a2, b2) {
            error_type = 4;
        }

        /* append the newly-created intervals to the list */
        f_w.update(a1, b1, area1, error1, a2, b2, area2, error2);

        if errsum <= tolerance {
            return compute_result(f_w, result, abserr, errsum, error_type);
        }

        if error_type != 0 {
            break;
        }

        if iteration >= limit - 1 {
            error_type = 1;
            break;
        }

        /* set up variables on first iteration */
        if iteration == 2 {
            error_over_large_intervals = errsum;
            ertest = tolerance;
            append_table(&mut table, area);
            continue;
        }

        if disallow_extrapolation != 0 {
            continue;
        }

        error_over_large_intervals += -last_e_i;

        if current_level < (*w).maximum_level {
            error_over_large_intervals += error12;
        }

        if extrapolate == 0 {
            /* test whether the interval to be bisected next is the
             smallest interval. */

            if large_interval(w) {
                continue;
            }

            extrapolate = 1;
            (*w).nrmax = 1;
        }

        if error_type2 == 0 && error_over_large_intervals > ertest {
            if increase_nrmax(w) {
                continue;
            }
        }

        /* Perform extrapolation */
        append_table(&mut table, area);

        intern_qelg(&mut table, &mut reseps, &mut abseps);

        ktmin += 1;

        if ktmin > 5 && err_ext < 0.001f64 * errsum {
            error_type = 5;
        }

        if abseps < err_ext {
            ktmin = 0;
            err_ext = abseps;
            res_ext = reseps;
            correc = error_over_large_intervals;
            ertest = epsabs.max(epsrel * fabsf64(reseps));
            if err_ext <= ertest {
                break;
            }
        }

        /* Prepare bisection of the smallest interval. */
        if table.n == 1 {
            disallow_extrapolation = 1;
        }

        if error_type == 5 {
            break;
        }

        /* work on interval with largest error */
        reset_nrmax(w);
        extrapolate = 0;
        error_over_large_intervals = errsum;
        if iteration >= limit {
            break;
        }
    }

    *result = res_ext;
    *abserr = err_ext;

    if err_ext == ::DBL_MAX {
        return compute_result(f_w, result, abserr, errsum, error_type);
    }

    if error_type != 0 || error_type2 != 0 {
        if error_type2 != 0 {
            err_ext += correc;
        }

        if error_type == 0 {
            error_type = 3;
        }

        if res_ext != 0f64 && area != 0f64 {
            if err_ext / fabsf64(res_ext) > errsum / fabsf64(area) {
                return compute_result(f_w, result, abserr, errsum, error_type);
            }
        } else if err_ext > errsum {
            return compute_result(f_w, result, abserr, errsum, error_type);
        } else if area == 0f64 {
            return return_error(error_type);
        }
    }

    /*  Test on divergence. */
    let max_area = fabsf64(res_ext).max(fabsf64(area));

    if !positive_integrand && max_area < 0.01f64 * resabs0 {
        return return_error(error_type);
    }

    let ratio = res_ext / area;
    if ratio < 0.01f64 || ratio > 100f64 || errsum > fabsf64(area) {
        error_type = 6;
    }
    return_error(error_type)
}

unsafe fn intern_qagp<T>(f: ::function<T>, arg: &mut T, pts: &mut [f64], epsabs: f64, epsrel: f64, limit: u64, f_w: &IntegrationWorkspace,
    result: &mut f64, abserr: &mut f64, q: ::integration_function<T>) -> enums::Value {
    let w = f_w.w;
    let mut reseps = 0f64;
    let mut abseps = 0f64;
    let mut correc = 0f64;
    let mut ktmin = 0u64;
    let mut roundoff_type1 = 0i32;
    let mut roundoff_type2 = 0i32;
    let mut roundoff_type3 = 0i32;
    let mut error_type = 0i32;
    let mut error_type2 = 0i32;

    let mut extrapolate = 0i32;
    let mut disallow_extrapolation = 0i32;

    let mut table : ffi::extrapolation_table = ::std::mem::zeroed();

    /* number of intervals */
    let nint = pts.len() as u64 - 1u64;

    /* temporarily alias ndin to level */
    let mut t_ndin = CVec::new((*w).level, pts.len());
    let ndin = t_ndin.as_mut_slice();

    /* Initialize results */
    *result = 0f64;
    *abserr = 0f64;

    /* Test on validity of parameters */
    if limit > (*w).limit {
        rgsl_error!("iteration limit exceeds available workspace", enums::Inval);
    }

    if pts.len() as u64 > (*w).limit {
        rgsl_error!("pts length exceeds size of workspace", enums::Inval);
    }

    if epsabs <= 0f64 && (epsrel < 50f64 * ::DBL_EPSILON || epsrel < 0.5e-28f64) {
        rgsl_error!("tolerance cannot be acheived with given epsabs and epsrel", enums::BadTol);
    }

    /* Check that the integration range and break points are an ascending sequence */
    for i in range(0u, nint as uint) {
        if pts[i + 1] < pts[i] {
            rgsl_error!("points are not in an ascending sequence", enums::Inval);
        }
    }

    /* Perform the first integration */
    let mut result0 = 0f64;
    let mut abserr0 = 0f64;
    let mut resabs0 = 0f64;

    f_w.initialise(0f64, 0f64);

    for i in range(0u, nint as uint) {
        let mut area1 = 0f64;
        let mut error1 = 0f64;
        let mut resabs1 = 0f64;
        let mut resasc1 = 0f64;
        let a1 = pts[i];
        let b1 = pts[i + 1];

        q(f, arg, a1, b1, &mut area1, &mut error1, &mut resabs1, &mut resasc1);

        result0 = result0 + area1;
        abserr0 = abserr0 + error1;
        resabs0 = resabs0 + resabs1;

        f_w.append_interval(a1, b1, area1, error1);

        if error1 == resasc1 && error1 != 0f64 {
            ndin[i] = 1;
        } else {
            ndin[i] = 0;
        }
    }

    /* Compute the initial error estimate */
    let mut errsum = 0f64;
    let mut t_elist = CVec::new((*w).elist, nint as uint);
    let elist = t_elist.as_mut_slice();
    let mut t_level = CVec::new((*w).level, nint as uint);
    let level = t_level.as_mut_slice();

    for i in range(0u, nint as uint) {
        if ndin[i] != 0 {
            elist[i] = abserr0;
        }
        errsum = errsum + elist[i];
    }

    for i in range(0u, nint as uint) {
        level[i] = 0u64;
    }

    /* Sort results into order of decreasing error via the indirection array order[] */
    f_w.sort_results();

    /* Test on accuracy */
    let mut tolerance = epsabs.max(epsrel * fabsf64(result0));

    if abserr0 <= 100f64 * ::DBL_EPSILON * resabs0 && abserr0 > tolerance {
        *result = result0;
        *abserr = abserr0;

        rgsl_error!("cannot reach tolerance because of roundoff error on first attempt", enums::Round);
    } else if abserr0 <= tolerance {
        *result = result0;
        *abserr = abserr0;

        return enums::Success;
    } else if limit == 1 {
        *result = result0;
        *abserr = abserr0;

        rgsl_error!("a maximum of one iteration was insufficient", enums::MaxIter);
    }

    /* Initialization */
    initialise_table(&mut table);
    append_table(&mut table, result0);

    let mut area = result0;

    let mut res_ext = result0;
    let mut err_ext = ::DBL_MAX;

    let mut error_over_large_intervals = errsum;
    let mut ertest = tolerance;

    let positive_integrand = test_positivity(result0, resabs0);

    let mut iteration = nint - 1; 

    loop {
        let mut a_i = 0f64;
        let mut b_i = 0f64;
        let mut r_i = 0f64;
        let mut e_i = 0f64;
        let mut area1 = 0f64;
        let mut area2 = 0f64;
        let mut error1 = 0f64;
        let mut error2 = 0f64;
        let mut resasc1 = 0f64;
        let mut resasc2 = 0f64;
        let mut resabs1 = 0f64;
        let mut resabs2 = 0f64;

        /* Bisect the subinterval with the largest error estimate */
        f_w.retrieve(&mut a_i, &mut b_i, &mut r_i, &mut e_i);

        let current_level = level[(*w).i as uint] + 1u64;

        let a1 = a_i;
        let b1 = 0.5f64 * (a_i + b_i);
        let a2 = b1;
        let b2 = b_i;

        iteration += 1;

        q(f, arg, a1, b1, &mut area1, &mut error1, &mut resabs1, &mut resasc1);
        q(f, arg, a2, b2, &mut area2, &mut error2, &mut resabs2, &mut resasc2);

        let area12 = area1 + area2;
        let error12 = error1 + error2;
        let last_e_i = e_i;

        /* Improve previous approximations to the integral and test for accuracy.

         We write these expressions in the same way as the original QUADPACK code so that the rounding errors are the same, which
         makes testing easier. */

        errsum = errsum + error12 - e_i;
        area = area + area12 - r_i;

        tolerance = epsabs.max(epsrel * fabsf64(area));

        if resasc1 != error1 && resasc2 != error2 {
            let delta = r_i - area12;

            if fabsf64(delta) <= 1.0e-5f64 * fabsf64(area12) && error12 >= 0.99f64 * e_i {
                if extrapolate == 0 {
                    roundoff_type1 += 1;
                } else {
                  roundoff_type2 += 1;
                }
            }

            if nint + 1 > 10 && error12 > e_i {
                roundoff_type3 += 1;
            }
        }

        /* Test for roundoff and eventually set error flag */
        if roundoff_type1 + roundoff_type2 >= 10 || roundoff_type3 >= 20 {
            /* round off error */
            error_type = 2;
        }

        if roundoff_type2 >= 5 {
            error_type2 = 1;
        }

        /* set error flag in the case of bad integrand behaviour at a point of the integration range */
        if ::util::subinterval_too_small(a1, a2, b2) {
            error_type = 4;
        }

        /* append the newly-created intervals to the list */
        f_w.update(a1, b1, area1, error1, a2, b2, area2, error2);

        if errsum <= tolerance {
            return compute_result(f_w, result, abserr, errsum, error_type);
        }

        if error_type != 0 {
            break;
        }

        if iteration >= limit - 1 {
            error_type = 1;
            break;
        }

        if disallow_extrapolation != 0 {
            continue;
        }

        error_over_large_intervals += -last_e_i;

        if current_level < (*w).maximum_level {
            error_over_large_intervals += error12;
        }

        if extrapolate == 0 {
            /* test whether the interval to be bisected next is the smallest interval. */
            if large_interval(w) {
                continue;
            }

            extrapolate = 1;
            (*w).nrmax = 1;
        }

        /* The smallest interval has the largest error.  Before bisecting decrease the sum of the errors over the larger
         intervals (error_over_large_intervals) and perform extrapolation. */
        if error_type2 == 0 && error_over_large_intervals > ertest {
            if increase_nrmax(w) {
                continue;
            }
        }

        /* Perform extrapolation */
        append_table (&mut table, area);

        if table.n > 2 {
            intern_qelg(&mut table, &mut reseps, &mut abseps);

            ktmin += 1;

            if ktmin > 5 && err_ext < 0.001f64 * errsum {
                error_type = 5;
            }

            if abseps < err_ext {
                ktmin = 0;
                err_ext = abseps;
                res_ext = reseps;
                correc = error_over_large_intervals;
                ertest = epsabs.max(epsrel * fabsf64(reseps));
                if err_ext <= ertest {
                    break;
                }
            }

            /* Prepare bisection of the smallest interval. */
            if table.n == 1 {
                disallow_extrapolation = 1;
            }

            if error_type == 5 {
                break;
            }
        }

        reset_nrmax(w);
        extrapolate = 0;
        error_over_large_intervals = errsum;
        if iteration >= limit {
            break;
        }
    }

    *result = res_ext;
    *abserr = err_ext;

    if err_ext == ::DBL_MAX {
        return compute_result(f_w, result, abserr, errsum, error_type);
    }

    if error_type != 0 || error_type2 != 0 {
        if error_type2 != 0 {
            err_ext += correc;
        }

        if error_type == 0 {
            error_type = 3;
        }

        if *result != 0f64 && area != 0f64 {
            if err_ext / fabsf64(res_ext) > errsum / fabsf64(area) {
                return compute_result(f_w, result, abserr, errsum, error_type);
            }
        } else if err_ext > errsum {
            return compute_result(f_w, result, abserr, errsum, error_type);
        } else if area == 0f64 {
            return return_error(error_type);
        }
    }

    /*  Test on divergence. */
    {
        let max_area = fabsf64(res_ext).max(fabsf64(area));

        if !positive_integrand && max_area < 0.01f64 * resabs0 {
            return return_error(error_type);
        }
    }

    {
        let ratio = res_ext / area;

        if ratio < 0.01f64 || ratio > 100f64 || errsum > fabsf64(area) {
            error_type = 6;
        }
    }

    return_error(error_type)
}

unsafe fn qc25c<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, c: f64, result: &mut f64, abserr: &mut f64, err_reliable: &mut i32) {
    let cc = (2f64 * c - b - a) / (b - a);

    if fabsf64(cc) > 1.1f64 {
        let mut resabs = 0f64;
        let mut resasc = 0f64;

        let mut tmp_arg = InternParam{func: f, param: arg, p2: c};

        ::integration::qk15(fn_cauchy, &mut tmp_arg, a, b, result, abserr, &mut resabs, &mut resasc);

        if *abserr == resasc {
            *err_reliable = 0;
        } else  {
            *err_reliable = 1;
        }
    } else {
        let mut cheb12 : [f64, ..13] = [0f64, ..13];
        let mut cheb24 : [f64, ..25] = [0f64, ..25];
        let mut moment : [f64, ..25] = [0f64, ..25];
        let mut res12 = 0f64;
        let mut res24 = 0f64;

        gsl_integration_qcheb(f, arg, a, b, cheb12, cheb24);
        compute_moments(cc, moment);

        for i in range(0, 13) {
            res12 += cheb12[i] * moment[i];
        }

        for i in range(0, 25) {
            res24 += cheb24[i] * moment[i];
        }

        *result = res24;
        *abserr = fabsf64(res24 - res12) ;
        *err_reliable = 0;
    }
}

fn fn_cauchy<T>(x: f64, p: &mut InternParam<T>) -> f64 {
    let c = p.p2;
    let f = p.func;

    f(x, p.param) / (x - c)
}

unsafe fn compute_moments(cc: f64, moment: &mut [f64]) {
    let mut a0 = logf64(fabsf64((1f64 - cc) / (1f64 + cc)));
    let mut a1 = 2f64 + a0 * cc;

    moment[0] = a0;
    moment[1] = a1;

    for k in range(2u, 25u) {
        let a2 = if (k & 1) == 0 {
            2f64 * cc * a1 - a0
        } else {
            let km1 = k as f64 - 1f64;

            2f64 * cc * a1 - a0 - 4.0 / (km1 * km1 - 1f64)
        };

        moment[k] = a2;

        a0 = a1;
        a1 = a2;
    }
}

// This function computes the 12-th order and 24-th order Chebyshev approximations to f(x) on [a,b]
fn gsl_integration_qcheb<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, cheb12: &mut [f64], cheb24: &mut [f64]) {
    let mut fval : [f64, ..25] = [0f64, ..25];
    let mut v : [f64, ..12] = [0f64, ..12];

    /* These are the values of cos(pi*k/24) for k=1..11 needed for the Chebyshev expansion of f(x) */
    let x : [f64, ..11] = [
        0.9914448613738104f64,     
        0.9659258262890683f64,
        0.9238795325112868f64,     
        0.8660254037844386f64,
        0.7933533402912352f64,     
        0.7071067811865475f64,
        0.6087614290087206f64,     
        0.5000000000000000f64,
        0.3826834323650898f64,     
        0.2588190451025208f64,
        0.1305261922200516f64];
  
    let center = 0.5f64 * (b + a);
    let half_length =  0.5f64 * (b - a);
  
    fval[0] = 0.5f64 * f(b, arg);
    fval[12] = f(center, arg);
    fval[24] = 0.5f64 * f(a, arg);

    for i in range(1u, 12u) {
        let j = 24u - i;
        let u = half_length * x[i - 1u];

        fval[i] = f(center + u, arg);
        fval[j] = f(center - u, arg);
    }

    for i in range(1u, 12u) {
        let j = 24u - i;

        v[i] = fval[i] - fval[j];
        fval[i] = fval[i] + fval[j];
    }

    {
        let alam1 = v[0] - v[8];
        let alam2 = x[5] * (v[2] - v[6] - v[10]);

        cheb12[3] = alam1 + alam2;
        cheb12[9] = alam1 - alam2;
    }

    {
        let alam1 = v[1] - v[7] - v[9];
        let alam2 = v[3] - v[5] - v[11];
        {
            let alam = x[2] * alam1 + x[8] * alam2;

            cheb24[3] = cheb12[3] + alam;
            cheb24[21] = cheb12[3] - alam;
        }

        {
            let alam = x[8] * alam1 - x[2] * alam2;
            cheb24[9] = cheb12[9] + alam;
            cheb24[15] = cheb12[9] - alam;
        }
    }

    {
        let part1 = x[3] * v[4];
        let part2 = x[7] * v[8];
        let part3 = x[5] * v[6];
        
        {
            let alam1 = v[0] + part1 + part2;
            let alam2 = x[1] * v[2] + part3 + x[9] * v[10];
          
            cheb12[1] = alam1 + alam2;
            cheb12[11] = alam1 - alam2;
        }
        
        {
            let alam1 = v[0] - part1 + part2;
            let alam2 = x[9] * v[2] - part3 + x[1] * v[10];

            cheb12[5] = alam1 + alam2;
            cheb12[7] = alam1 - alam2;
        }
    }

    {
        let alam = x[0] * v[1] + x[2] * v[3] + x[4] * v[5] + x[6] * v[7] + x[8] * v[9] + x[10] * v[11];

        cheb24[1] = cheb12[1] + alam;
        cheb24[23] = cheb12[1] - alam;
    }

    {
        let alam = x[10] * v[1] - x[8] * v[3] + x[6] * v[5] - x[4] * v[7] + x[2] * v[9] - x[0] * v[11];

        cheb24[11] = cheb12[11] + alam;
        cheb24[13] = cheb12[11] - alam;
    }

    {
        let alam = x[4] * v[1] - x[8] * v[3] - x[0] * v[5] - x[10] * v[7] + x[2] * v[9] + x[6] * v[11];

        cheb24[5] = cheb12[5] + alam;
        cheb24[19] = cheb12[5] - alam;
    }

    {
        let alam = x[6] * v[1] - x[2] * v[3] - x[10] * v[5] + x[0] * v[7] - x[8] * v[9] - x[4] * v[11];

        cheb24[7] = cheb12[7] + alam;
        cheb24[17] = cheb12[7] - alam;
    }

    for i in range(0u, 6u) {
        let j = 12u - i;

        v[i] = fval[i] - fval[j];
        fval[i] = fval[i] + fval[j];
    }

    {
        let alam1 = v[0] + x[7] * v[4];
        let alam2 = x[3] * v[2];

        cheb12[2] = alam1 + alam2;
        cheb12[10] = alam1 - alam2;
    }

    cheb12[6] = v[0] - v[4];

    {
        let alam = x[1] * v[1] + x[5] * v[3] + x[9] * v[5];

        cheb24[2] = cheb12[2] + alam;
        cheb24[22] = cheb12[2] - alam;
    }

    {
        let alam = x[5] * (v[1] - v[3] - v[5]);

        cheb24[6] = cheb12[6] + alam;
        cheb24[18] = cheb12[6] - alam;
    }

    {
        let alam = x[9] * v[1] - x[5] * v[3] + x[1] * v[5];

        cheb24[10] = cheb12[10] + alam;
        cheb24[14] = cheb12[10] - alam;
    }

    for i in range(0u, 3u) {
        let j = 6 - i;

        v[i] = fval[i] - fval[j];
        fval[i] = fval[i] + fval[j];
    }

    cheb12[4] = v[0] + x[7] * v[2];
    cheb12[8] = fval[0] - x[7] * fval[2];

    {
        let alam = x[3] * v[1];

        cheb24[4] = cheb12[4] + alam;
        cheb24[20] = cheb12[4] - alam;
    }

    {
        let alam = x[7] * fval[1] - fval[3];

        cheb24[8] = cheb12[8] + alam;
        cheb24[16] = cheb12[8] - alam;
    }

    cheb12[0] = fval[0] + fval[2];

    {
        let alam = fval[1] + fval[3];

        cheb24[0] = cheb12[0] + alam;
        cheb24[24] = cheb12[0] - alam;
    }

    cheb12[12] = v[0] - v[2];
    cheb24[12] = cheb12[12];

    let mut tmp = 1f64 / 6f64;

    for i in range(0u, 12u) {
        cheb12[i] *= tmp;
    }

    tmp = 1f64 / 12f64;

    cheb12[0] *= tmp;
    cheb12[12] *= tmp;

    for i in range(1u, 24u) {
        cheb24[i] *= tmp;
    }

    tmp = 1f64 / 24f64;
    cheb24[0] *= tmp;
    cheb24[24] *= tmp;
}

struct fn_qaws_params<'r, T:'r> {
    function: ::function<T>,
    a: f64,
    b: f64,
    table: *mut ffi::gsl_integration_qaws_table,
    arg: &'r mut T
}

unsafe fn qc25s<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, a1: f64, b1: f64, t: *mut ffi::gsl_integration_qaws_table, result: &mut f64,
    abserr: &mut f64, err_reliable: &mut bool)
{
    let mut fn_params = fn_qaws_params{function: f, a: a, b: b, table: t, arg: arg};

    if a1 == a && ((*t).alpha != 0f64 || (*t).mu != 0) {
        let mut cheb12 : [f64, ..13] = [0f64, ..13];
        let mut cheb24 : [f64, ..25] = [0f64, ..25];

        let factor = powf64(0.5f64 * (b1 - a1), (*t).alpha + 1f64);

        gsl_integration_qcheb(fn_qaws_R, &mut fn_params, a1, b1, cheb12, cheb24);

        if (*t).mu == 0 {
            let mut res12 = 0f64;
            let mut res24 = 0f64;
            let u = factor;

            qc25s_compute_result((*t).ri, cheb12, cheb24, &mut res12, &mut res24);

            *result = u * res24;
            *abserr = fabsf64(u * (res24 - res12));
        } else {
            let mut res12a = 0f64;
            let mut res24a = 0f64;
            let mut res12b = 0f64;
            let mut res24b = 0f64;

            let u = factor * logf64(b1 - a1);
            let v = factor;

            qc25s_compute_result((*t).ri, cheb12, cheb24, &mut res12a, &mut res24a);
            qc25s_compute_result((*t).rg, cheb12, cheb24, &mut res12b, &mut res24b);

            *result = u * res24a + v * res24b;
            *abserr = fabsf64(u * (res24a - res12a)) + fabsf64(v * (res24b - res12b));
        }

        *err_reliable = false;
    } else if b1 == b && ((*t).beta != 0.0 || (*t).nu != 0) {
        let mut cheb12 : [f64, ..13] = [0f64, ..13];
        let mut cheb24 : [f64, ..25] = [0f64, ..25];
        let factor = powf64(0.5f64 * (b1 - a1), (*t).beta + 1f64);

        gsl_integration_qcheb(fn_qaws_L, &mut fn_params, a1, b1, cheb12, cheb24);

        if (*t).nu == 0 {
            let mut res12 = 0f64;
            let mut res24 = 0f64;
            let u = factor;

            qc25s_compute_result((*t).rj, cheb12, cheb24, &mut res12, &mut res24);

            *result = u * res24;
            *abserr = fabsf64(u * (res24 - res12));
        } else {
            let mut res12a = 0f64;
            let mut res24a = 0f64;
            let mut res12b = 0f64;
            let mut res24b = 0f64;

            let u = factor * logf64(b1 - a1);
            let v = factor;

            qc25s_compute_result((*t).rj, cheb12, cheb24, &mut res12a, &mut res24a);
            qc25s_compute_result((*t).rh, cheb12, cheb24, &mut res12b, &mut res24b);

            *result = u * res24a + v * res24b;
            *abserr = fabsf64(u * (res24a - res12a)) + fabsf64(v * (res24b - res12b));
        }
        *err_reliable = false;
    }
    else {
        let mut resabs = 0f64;
        let mut resasc = 0f64;

        ::integration::qk15(fn_qaws, &mut fn_params, a1, b1, result, abserr, &mut resabs, &mut resasc);

        if *abserr == resasc {
            *err_reliable = false;
        } else {
            *err_reliable = true;
        }
    }
}

fn fn_qaws<T>(x: f64, p: &mut fn_qaws_params<T>) -> f64 {
    let f = p.function;
    let t = p.table;

    let mut factor = 1f64;
  
    unsafe {
        if (*t).alpha != 0f64 {
            factor *= powf64(x - p.a, (*t).alpha);
        }

        if (*t).beta != 0f64 {
            factor *= powf64(p.b - x, (*t).beta);
        }

        if (*t).mu == 1 {
            factor *= logf64(x - p.a);
        }

        if (*t).nu == 1 {
            factor *= logf64(p.b - x);
        }

        factor * f(x, p.arg)
    }
}

fn fn_qaws_R<T>(x: f64, p: &mut fn_qaws_params<T>) -> f64 {
    let f = p.function;
    let t = p.table;

    let mut factor = 1f64;
  
    unsafe {
        if (*t).beta != 0f64 {
            factor *= powf64(p.b - x, (*t).beta);
        }

        if (*t).nu == 1 {
            factor *= logf64(p.b - x);
        }

        factor * f(x, p.arg)
    }
}

fn fn_qaws_L<T>(x: f64, p: &mut fn_qaws_params<T>) -> f64 {
    let f = p.function;
    let t = p.table;

    let mut factor = 1f64;
  
    unsafe {
        if (*t).alpha != 0f64 {
            factor *= powf64(x - p.a, (*t).alpha);
        }

        if (*t).mu == 1 {
            factor *= logf64(x - p.a);
        }

        factor * f(x, p.arg)
    }
}

unsafe fn qc25s_compute_result(r: &[f64], cheb12: &[f64], cheb24: &[f64], result12: &mut f64, result24: &mut f64) {
    let mut res12 = 0f64;
    let mut res24 = 0f64;

    for i in range(0u, 13u) {
        res12 += r[i] * cheb12[i];
    }

    for i in range(0u, 25u) {
        res24 += r[i] * cheb24[i];
    }

    *result12 = res12;
    *result24 = res24;
}