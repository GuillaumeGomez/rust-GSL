//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use enums;
use std::intrinsics::fabsf64;
use std::c_vec::CVec;

/// The QAG algorithm is a simple adaptive integration procedure. The integration region is divided into subintervals, and on each iteration
/// the subinterval with the largest estimated error is bisected. This reduces the overall error rapidly, as the subintervals become concentrated
/// around local difficulties in the integrand. These subintervals are managed by a gsl_integration_workspace struct, which handles the memory
/// for the subinterval ranges, results and error estimates.
pub struct IntegrationWorkspace {
    w: *mut ffi::gsl_integration_workspace
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
        intern_qag(f, arg, a, b, epsabs, epsrel, limit, self.w, result, abserr, integration_rule)
    }

    /// This function applies the Gauss-Kronrod 21-point integration rule adaptively until an estimate of the integral of f over (a,b) is achieved
    /// within the desired absolute and relative error limits, epsabs and epsrel. The results are extrapolated using the epsilon-algorithm, which
    /// accelerates the convergence of the integral in the presence of discontinuities and integrable singularities. The function returns the
    /// final approximation from the extrapolation, result, and an estimate of the absolute error, abserr. The subintervals and their results are
    /// stored in the memory provided by workspace. The maximum number of subintervals is given by limit, which may not exceed the allocated size
    /// of the workspace.
    pub fn qags<T>(&self, f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, result: &mut f64,
        abserr: &mut f64) -> enums::Value {
        unsafe { intern_qags(f, arg, a, b, epsabs, epsrel, limit, self.w, result, abserr, ::integration::qk21) }
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

fn initialise(w: *mut ffi::gsl_integration_workspace, a: f64, b: f64) {
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

fn set_initial_result(w: *mut ffi::gsl_integration_workspace, result: f64, error: f64) {
    unsafe {
        let mut rlist = CVec::new((*w).rlist, 1);
        let mut elist = CVec::new((*w).elist, 1);

        (*w).size = 1;
        rlist.as_mut_slice()[0] = result;
        elist.as_mut_slice()[0] = error;
    }
}

fn intern_qag<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, w: *mut ffi::gsl_integration_workspace,
    result: &mut f64, abserr: &mut f64, q: ::integration_function<T>) -> enums::Value {
    let mut roundoff_type1 = 0i32;
    let mut roundoff_type2 = 0i32;
    let mut error_type = 0i32;

    *result = 0f64;
    *abserr = 0f64;

    // Initialize results
    initialise(w, a, b);

    if unsafe { limit > (*w).limit } {
        let file = file!();
        "iteration limit exceeds available workspace".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::Inval as i32) }
            });
        });
    }
    if epsabs <= 0f64 && (epsrel < 50f64 * ::DBL_EPSILON || epsrel < 0.5e-28f64) {
        let file = file!();
        "tolerance cannot be achieved with given epsabs and epsrel".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::BadTol as i32) }
            });
        });
    }

    // perform the first integration
    let mut result0 = 0f64;
    let mut abserr0 = 0f64;
    let mut resabs0 = 0f64;
    let mut resasc0 = 0f64;
    q(f, arg, a, b, &mut result0, &mut abserr0, &mut resabs0, &mut resasc0);

    set_initial_result(w, result0, abserr0);

    // Test on accuracy
    let mut tolerance = unsafe { epsabs.max(epsrel * fabsf64(result0)) };

    // need IEEE rounding here to match original quadpack behavior
    let round_off = 50f64 * ::DBL_EPSILON * resabs0;

    if abserr0 <= round_off && abserr0 > tolerance {
        *result = result0;
        *abserr = abserr0;

        let file = file!();
        "cannot reach tolerance because of roundoff error on first attempt".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::Round as i32) }
            });
        });
    } else if (abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0f64 {
        *result = result0;
        *abserr = abserr0;

        return enums::Success;
    } else if limit == 1 {
        *result = result0;
        *abserr = abserr0;

        let file = file!();
        "ca maximum of one iteration was insufficient".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::MaxIter as i32) }
            });
        });
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
        ::util::retrieve(&ffi::FFI::wrap(w), &mut a_i, &mut b_i, &mut r_i, &mut e_i);

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
    *result = ::util::sum_results(&ffi::FFI::wrap(w));
    *abserr = errsum;

    if errsum <= tolerance {
        enums::Success
    } else if error_type == 2 {
        let file = file!();
        "roundoff error prevents tolerance from being achieved".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::Round as i32) }
            });
        });
        enums::Round
    } else if error_type == 3 {
        let file = file!();
        "bad integrand behavior found in the integration interval".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::Sing as i32) }
            });
        });
        enums::Sing
    } else if iteration == limit {
        let file = file!();
        "maximum number of subdivisions reached".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::MaxIter as i32) }
            });
        });
        enums::MaxIter
    } else {
        let file = file!();
        "could not integrate function".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::Failed as i32) }
            });
        });
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

unsafe fn compute_result(w: *mut ffi::gsl_integration_workspace, result: &mut f64, abserr: &mut f64, errsum: f64,
    error_type: i32) -> enums::Value {
    *result = ::util::sum_results(&ffi::FFI::wrap(w));
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
            let file = file!();
            "number of iterations was insufficient".with_c_str(|c_str|{
                file.with_c_str(|c_file|{
                    ffi::gsl_error(c_str, c_file, line!() as i32, enums::MaxIter as i32)
                });
            });
            enums::MaxIter
        }
        2 => {
            let file = file!();
            "cannot reach tolerance because of roundoff error".with_c_str(|c_str|{
                file.with_c_str(|c_file|{
                    ffi::gsl_error(c_str, c_file, line!() as i32, enums::Round as i32)
                });
            });
            enums::Round
        }
        3 => {
            let file = file!();
            "bad integrand behavior found in the integration interval".with_c_str(|c_str|{
                file.with_c_str(|c_file|{
                    ffi::gsl_error(c_str, c_file, line!() as i32, enums::Sing as i32)
                });
            });
            enums::Sing
        }
        4 => {
            let file = file!();
            "roundoff error detected in the extrapolation table".with_c_str(|c_str|{
                file.with_c_str(|c_file|{
                    ffi::gsl_error(c_str, c_file, line!() as i32, enums::Round as i32)
                });
            });
            enums::Round
        }
        5 => {
            let file = file!();
            "integral is divergent, or slowly convergent".with_c_str(|c_str|{
                file.with_c_str(|c_file|{
                    ffi::gsl_error(c_str, c_file, line!() as i32, enums::Round as i32)
                });
            });
            enums::Round
        }
        _ => {
            let file = file!();
            "integral is divergent, or slowly convergent".with_c_str(|c_str|{
                file.with_c_str(|c_file|{
                    ffi::gsl_error(c_str, c_file, line!() as i32, enums::Failed as i32)
                });
            });
            enums::Failed
        }
    }
}

unsafe fn intern_qags<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, w: *mut ffi::gsl_integration_workspace,
    result: &mut f64, abserr: &mut f64, q: ::integration_function<T>) -> enums::Value {
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
    initialise (w, a, b);

    *result = 0f64;
    *abserr = 0f64;

    if limit > (*w).limit {
        let file = file!();
        "iteration limit exceeds available workspace".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                ffi::gsl_error(c_str, c_file, line!() as i32, enums::Inval as i32)
            });
        });
    }

    /* Test on accuracy */
    if epsabs <= 0f64 && (epsrel < 50f64 * ::DBL_EPSILON || epsrel < 0.5e-28f64) {
        let file = file!();
        "tolerance cannot be acheived with given epsabs and epsrel".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                ffi::gsl_error(c_str, c_file, line!() as i32, enums::BadTol as i32)
            });
        });
    }

    /* Perform the first integration */
    q(f, arg, a, b, &mut result0, &mut abserr0, &mut resabs0, &mut resasc0);

    set_initial_result(w, result0, abserr0);

    let mut tolerance = epsabs.max(epsrel * fabsf64(result0));

    if abserr0 <= 100f64 * ::DBL_EPSILON * resabs0 && abserr0 > tolerance {
        *result = result0;
        *abserr = abserr0;

        let file = file!();
        "cannot reach tolerance because of roundoff error on first attempt".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                ffi::gsl_error(c_str, c_file, line!() as i32, enums::Round as i32)
            });
        });
    } else if (abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0f64 {
        *result = result0;
        *abserr = abserr0;

        return enums::Success;
    } else if limit == 1 {
        *result = result0;
        *abserr = abserr0;

        let file = file!();
        "a maximum of one iteration was insufficient".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                ffi::gsl_error(c_str, c_file, line!() as i32, enums::MaxIter as i32)
            });
        });
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
        ::util::retrieve(&ffi::FFI::wrap(w), &mut a_i, &mut b_i, &mut r_i, &mut e_i);

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
        ::util::update(&ffi::FFI::wrap(w), a1, b1, area1, error1, a2, b2, area2, error2);

        if errsum <= tolerance {
            return compute_result(w, result, abserr, errsum, error_type);
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
        return compute_result(w, result, abserr, errsum, error_type);
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
                return compute_result(w, result, abserr, errsum, error_type);
            }
        } else if err_ext > errsum {
            return compute_result(w, result, abserr, errsum, error_type);
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