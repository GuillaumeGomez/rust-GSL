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
    pub fn qags<T>(&self, f: ::function<T>, arg: &mut T, a: f64, b: f64, epsabs: f64, epsrel: f64, limit: u64, key: ::GaussKonrodRule,
        result: &mut f64, abserr: &mut f64) -> enums::Value {
        enums::Success
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