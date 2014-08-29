//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use std::intrinsics::fabsf64;
use std::c_vec::CVec;

pub fn subinterval_too_small(a1: f64, a2: f64, b2: f64) -> bool {
    let e = ::DBL_EPSILON;
    let u = ::DBL_MIN;

    let tmp = unsafe { (1f64 + 100f64 * e) * (fabsf64(a2) + 1000f64 * u) };

    unsafe { fabsf64(a1) <= tmp && fabsf64(b2) <= tmp }
}

pub fn sum_results(w: &::IntegrationWorkspace) -> f64 {
    unsafe {
        let f_w = ffi::FFI::unwrap(w);
        let mut result_sum = 0f64;
        let v_rlist = CVec::new((*f_w).rlist, (*f_w).size as uint);
        let rlist = v_rlist.as_slice();

        for k in range(0u, (*f_w).size as uint) {
            result_sum += rlist[k];
        }

        result_sum
    }
}

pub fn retrieve(f_w: &::IntegrationWorkspace, a: &mut f64, b: &mut f64, r: &mut f64, e: &mut f64) {
    unsafe {
        let w = ffi::FFI::unwrap(f_w);
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

pub fn update(f_w: &::IntegrationWorkspace, a1: f64, b1: f64, area1: f64, error1: f64, a2: f64, b2: f64, area2: f64, error2: f64) {
    let w = ffi::FFI::unwrap(f_w);

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

        qpsrt(f_w);
    }
}

pub fn qpsrt(f_w: &::IntegrationWorkspace) {
    let w = ffi::FFI::unwrap(f_w);

    unsafe {
        let mut order = CVec::new((*w).order, (*w).nrmax as uint + 1u);
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

        let elist = CVec::new((*w).elist, i_maxerr as uint + 1u);
        let errmax = elist.as_slice()[i_maxerr as uint];

        // This part of the routine is only executed if, due to a difficult integrand, subdivision increased the error estimate. In the normal
        // case the insert procedure should start after the nrmax-th largest error estimate.
        while i_nrmax > 0 && errmax > elist.as_slice()[order.as_slice()[i_nrmax as uint - 1u] as uint] {
            order.as_mut_slice()[i_nrmax as uint] = order.as_slice()[i_nrmax as uint - 1u];
            i_nrmax -= 1;
        }

        // Compute the number of elements in the list to be maintained in descending order. This number depends on the number of
        // subdivisions still allowed.
        let top =  if last < (limit / 2 + 2) {
            last
        } else {
            limit - last + 1
        };

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