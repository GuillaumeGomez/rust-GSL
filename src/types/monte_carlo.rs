//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Monte Carlo Integration

This chapter describes routines for multidimensional Monte Carlo integration. These include the traditional Monte Carlo method and adaptive 
algorithms such as VEGAS and MISER which use importance sampling and stratified sampling techniques. Each algorithm computes an estimate of 
a multidimensional definite integral of the form,

I = \int_xl^xu dx \int_yl^yu  dy ...  f(x, y, ...)

over a hypercubic region ((x_l,x_u), (y_l,y_u), ...) using a fixed number of function calls. The routines also provide a statistical estimate 
of the error on the result. This error estimate should be taken as a guide rather than as a strict error bound—random sampling of the region 
may not uncover all the important features of the function, resulting in an underestimate of the error.

##Interface

All of the Monte Carlo integration routines use the same general form of interface. There is an allocator to allocate memory for control 
variables and workspace, a routine to initialize those control variables, the integrator itself, and a function to free the space when done.

Each integration function requires a random number generator to be supplied, and returns an estimate of the integral and its standard deviation. 
The accuracy of the result is determined by the number of function calls specified by the user. If a known level of accuracy is required this 
can be achieved by calling the integrator several times and averaging the individual results until the desired accuracy is obtained.

Random sample points used within the Monte Carlo routines are always chosen strictly within the integration region, so that endpoint singularities 
are automatically avoided.

##PLAIN Monte Carlo

The plain Monte Carlo algorithm samples points randomly from the integration region to estimate the integral and its error. Using this algorithm 
the estimate of the integral E(f; N) for N randomly distributed points x_i is given by,

E(f; N) = =  V <f> = (V / N) \sum_i^N f(x_i)
where V is the volume of the integration region. The error on this estimate \sigma(E;N) is calculated from the estimated variance of the mean,

\sigma^2 (E; N) = (V^2 / N^2) \sum_i^N (f(x_i) -  <f>)^2.
For large N this variance decreases asymptotically as \Var(f)/N, where \Var(f) is the true variance of the function over the integration region. 
The error estimate itself should decrease as \sigma(f)/\sqrt{N}. The familiar law of errors decreasing as 1/\sqrt{N} applies—to reduce the 
error by a factor of 10 requires a 100-fold increase in the number of sample points.
!*/

use ffi;
use enums;
use std::f64::INFINITY;
use std::intrinsics::sqrtf64;
use std::c_vec::CVec;

pub struct PlainMonteCarlo {
    s: *mut ffi::gsl_monte_plain_state,
}

impl PlainMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions.
    pub fn new(dim: u64) -> Option<PlainMonteCarlo> {
        let tmp = unsafe { ffi::gsl_monte_plain_alloc(dim) };

        if tmp.is_null() {
            None
        } else {
            Some(PlainMonteCarlo {
                s: tmp
            })
        }
    }

    /// This function initializes a previously allocated integration state. This allows an existing workspace to be reused for different
    /// integrations.
    pub fn init(&self) -> enums::Value {
        unsafe { ffi::gsl_monte_plain_init(self.s) }
    }

    /// This routines uses the plain Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined
    /// by the lower and upper limits in the arrays xl and xu, each of the same size. The integration uses a fixed number of function calls
    /// calls, and obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied.
    /// The result of the integration is returned in result, with an estimated absolute error abserr.
    pub fn integrate<T>(&self, f: ::monte_function<T>, arg: &mut T, xl: &[f64], xu: &[f64], calls: u64, r: &::Rng, result: &mut f64,
        abserr: &mut f64) -> enums::Value {
        unsafe {
            let mut m = 0f64;
            let mut q = 0f64;
            let dim = (*self.s).dim as uint;
            let mut t_x = CVec::new((*self.s).x, dim);
            let x = t_x.as_mut_slice();

            if xl.len() != dim || xu.len() != dim {
                rgsl_error!("number of dimensions must match allocated size", enums::Inval);
            }

            for i in range(0u, dim as uint) {
                if xu[i] <= xl[i] {
                    rgsl_error!("xu must be greater than xl", enums::Inval);
                }

                if xu[i] - xl[i] > ::DBL_MAX {
                    rgsl_error!("Range of integration is too large, please rescale", enums::Inval);
                }
            }

            /* Compute the volume of the region */
            let mut vol = 1f64;

            for i in range(0u, dim) {
                vol *= xu[i] - xl[i];
            }

            for n in range(0u, calls as uint) {
                /* Choose a random point in the integration region */
                for i in range(0u, dim) {
                    x[i] = xl[i] + r.uniform_pos() * (xu[i] - xl[i]);
                }

                {
                    let fval = f(x, arg);

                    /* recurrence for mean and variance */
                    let d = fval - m;
                    m += d / (n as f64 + 1f64);
                    q += d * d * (n as f64 / (n as f64 + 1f64));
                }
            }

            *result = vol * m;

            *abserr = if calls < 2 {
                INFINITY
            } else {
                vol * sqrtf64(q / (calls as f64 * (calls as f64 - 1f64)))
            };

            enums::Success
        }
    }
}

impl Drop for PlainMonteCarlo {
    fn drop(&mut self) {
        unsafe { ffi::gsl_monte_plain_free(self.s) };
        self.s = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_monte_plain_state> for PlainMonteCarlo {
    fn wrap(s: *mut ffi::gsl_monte_plain_state) -> PlainMonteCarlo {
        PlainMonteCarlo {
            s: s
        }
    }

    fn unwrap(s: &PlainMonteCarlo) -> *mut ffi::gsl_monte_plain_state {
        s.s
    }
}