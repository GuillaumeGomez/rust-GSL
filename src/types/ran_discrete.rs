//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

pub struct RanDiscrete {
    ran: *mut ffi::gsl_ran_discrete_t
}

impl RanDiscrete {
    pub fn new(P: &[f64]) -> Option<RanDiscrete> {
        let tmp = unsafe { ffi::gsl_ran_discrete_preproc(P.len() as u64, P.as_ptr()) };

        if tmp.is_null() {
            None
        } else {
            Some(RanDiscrete {
                ran: tmp
            })
        }
    }

    pub fn discrete(&self, r: &Rng) -> u64 {
        unsafe { ffi::gsl_ran_discrete(ffi::FFI::unwrap(r) as *const ffi::gsl_rng, self.ran as *const ffi::gsl_ran_discrete_t) }
    }

    pub fn discrete_pdf(&self, k: u64) -> f64 {
        unsafe { ffi::gsl_ran_discrete_pdf(k, self.ran as *const ffi::gsl_ran_discrete_t) }
    }
}

impl Drop for RanDiscrete {
    fn drop(&mut self) {
        unsafe { ffi::gsl_ran_discrete_free(self.ran) };
        self.ran = ::std::ptr::mut_null();
    }
}

impl ffi::FFI<ffi::gsl_ran_discrete_t> for RanDiscrete {
    fn wrap(r: *mut ffi::gsl_ran_discrete_t) -> RanDiscrete {
        RanDiscrete {
            ran: r
        }
    }

    fn unwrap(v: &RanDiscrete) -> *mut ffi::gsl_ran_discrete_t {
        v.ran
    }
}