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

##VEGAS

The VEGAS algorithm of Lepage is based on importance sampling. It samples points from the probability distribution described by the function 
|f|, so that the points are concentrated in the regions that make the largest contribution to the integral.

In general, if the Monte Carlo integral of f is sampled with points distributed according to a probability distribution described by the function 
g, we obtain an estimate E_g(f; N),

E_g(f; N) = E(f/g; N)
with a corresponding variance,

\Var_g(f; N) = \Var(f/g; N).
If the probability distribution is chosen as g = |f|/I(|f|) then it can be shown that the variance V_g(f; N) vanishes, and the error in the 
estimate will be zero. In practice it is not possible to sample from the exact distribution g for an arbitrary function, so importance sampling 
algorithms aim to produce efficient approximations to the desired distribution.

The VEGAS algorithm approximates the exact distribution by making a number of passes over the integration region while histogramming the 
function f. Each histogram is used to define a sampling distribution for the next pass. Asymptotically this procedure converges to the desired 
distribution. In order to avoid the number of histogram bins growing like K^d the probability distribution is approximated by a separable 
function: g(x_1, x_2, ...) = g_1(x_1) g_2(x_2) ... so that the number of bins required is only Kd. This is equivalent to locating the 
peaks of the function from the projections of the integrand onto the coordinate axes. The efficiency of VEGAS depends on the validity of 
this assumption. It is most efficient when the peaks of the integrand are well-localized. If an integrand can be rewritten in a form which 
is approximately separable this will increase the efficiency of integration with VEGAS.

VEGAS incorporates a number of additional features, and combines both stratified sampling and importance sampling. The integration region 
is divided into a number of “boxes”, with each box getting a fixed number of points (the goal is 2). Each box can then have a fractional 
number of bins, but if the ratio of bins-per-box is less than two, Vegas switches to a kind variance reduction (rather than importance 
sampling).

The VEGAS algorithm computes a number of independent estimates of the integral internally, according to the iterations parameter described 
below, and returns their weighted average. Random sampling of the integrand can occasionally produce an estimate where the error is zero, 
particularly if the function is constant in some regions. An estimate with zero error causes the weighted average to break down and must 
be handled separately. In the original Fortran implementations of VEGAS the error estimate is made non-zero by substituting a small value 
(typically 1e-30). The implementation in GSL differs from this and avoids the use of an arbitrary constant—it either assigns the value a 
weight which is the average weight of the preceding estimates or discards it according to the following procedure,

current estimate has zero error, weighted average has finite error
The current estimate is assigned a weight which is the average weight of the preceding estimates.

current estimate has finite error, previous estimates had zero error
The previous estimates are discarded and the weighted averaging procedure begins with the current estimate.

current estimate has zero error, previous estimates had zero error
The estimates are averaged using the arithmetic mean, but no error is computed.
!*/

use ffi;
use libc::{c_void, c_double, size_t};
use std::slice;
use std::mem::transmute;

/// The plain Monte Carlo algorithm samples points randomly from the integration region to estimate the integral and its error. Using this algorithm 
/// the estimate of the integral E(f; N) for N randomly distributed points x_i is given by,
///
/// E(f; N) = =  V <f> = (V / N) \sum_i^N f(x_i)
/// where V is the volume of the integration region. The error on this estimate \sigma(E;N) is calculated from the estimated variance of the mean,
///
/// \sigma^2 (E; N) = (V^2 / N^2) \sum_i^N (f(x_i) -  <f>)^2.
///
/// For large N this variance decreases asymptotically as \Var(f)/N, where \Var(f) is the true variance of the function over the integration region. 
/// The error estimate itself should decrease as \sigma(f)/\sqrt{N}. The familiar law of errors decreasing as 1/\sqrt{N} applies—to reduce the 
/// error by a factor of 10 requires a 100-fold increase in the number of sample points.
pub struct PlainMonteCarlo {
    s: *mut ffi::gsl_monte_plain_state,
}

impl PlainMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions.
    pub fn new(dim: usize) -> Option<PlainMonteCarlo> {
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
    pub fn init(&self) -> ::Value {
        unsafe { ffi::gsl_monte_plain_init(self.s) }
    }

    /// This routines uses the plain Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined
    /// by the lower and upper limits in the arrays xl and xu, each of the same size. The integration uses a fixed number of function calls
    /// calls, and obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied.
    /// The result of the integration is returned in result, with an estimated absolute error abserr.
    ///
    /// It returns either Ok((result, abserr)) or Err(enums::Value).
    pub fn integrate<F: Fn(&[f64]) -> f64 + 'static>(&self, dim: usize, f: F, xl: &[f64], xu: &[f64],
                                                     t_calls: usize, r: &::Rng) -> Result<(f64, f64), ::Value> {
        unsafe {
            assert!(xl.len() == xu.len());
            let mut result = 0f64;
            let mut abserr = 0f64;
            let f: Box<Box<Fn(&[f64]) -> f64 + 'static>> = Box::new(Box::new(f));
            let mut func = ffi::gsl_monte_function {
                               f: transmute(monte_trampoline as usize),
                               dim: dim,
                               params: Box::into_raw(f) as *mut _,
                           };
            let ret = ffi::gsl_monte_plain_integrate(&mut func as *mut _ as *mut c_void,
                                                     xl.as_ptr(), xu.as_ptr(), xl.len(), t_calls,
                                                     ffi::FFI::unwrap(r), self.s,
                                                     (&mut result) as *mut c_double,
                                                     (&mut abserr) as *mut c_double);

            if ret == ::Value::Success {
                Ok((result, abserr))
            } else {
                Err(ret)
            }
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

/// The MISER algorithm of Press and Farrar is based on recursive stratified sampling. This technique aims to reduce the overall integration error 
/// by concentrating integration points in the regions of highest variance.
///
/// The idea of stratified sampling begins with the observation that for two disjoint regions a and b with Monte Carlo estimates of the integral 
/// E_a(f) and E_b(f) and variances \sigma_a^2(f) and \sigma_b^2(f), the variance \Var(f) of the combined estimate E(f) = (1/2) (E_a(f) + E_b(f)) 
/// is given by,
///
/// \Var(f) = (\sigma_a^2(f) / 4 N_a) + (\sigma_b^2(f) / 4 N_b).
///
/// It can be shown that this variance is minimized by distributing the points such that,
///
/// N_a / (N_a + N_b) = \sigma_a / (\sigma_a + \sigma_b).
///
/// Hence the smallest error estimate is obtained by allocating sample points in proportion to the standard deviation of the function in each 
/// sub-region.
///
/// The MISER algorithm proceeds by bisecting the integration region along one coordinate axis to give two sub-regions at each step. The direction 
/// is chosen by examining all d possible bisections and selecting the one which will minimize the combined variance of the two sub-regions. The 
/// variance in the sub-regions is estimated by sampling with a fraction of the total number of points available to the current step. The same 
/// procedure is then repeated recursively for each of the two half-spaces from the best bisection. The remaining sample points are allocated to 
/// the sub-regions using the formula for N_a and N_b. This recursive allocation of integration points continues down to a user-specified depth 
/// where each sub-region is integrated using a plain Monte Carlo estimate. These individual values and their error estimates are then combined 
/// upwards to give an overall result and an estimate of its error.
pub struct MiserMonteCarlo {
    s: *mut ffi::gsl_monte_miser_state,
}

impl MiserMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions. The workspace is used to maintain
    /// the state of the integration.
    pub fn new(dim: usize) -> Option<MiserMonteCarlo> {
        let tmp_pointer = unsafe { ffi::gsl_monte_miser_alloc(dim) };

        if tmp_pointer.is_null() {
            None
        } else {
            Some(MiserMonteCarlo {
                s: tmp_pointer
            })
        }
    }

    /// This function initializes a previously allocated integration state. This allows an existing workspace to be reused for different integrations.
    pub fn init(&self) -> ::Value {
        unsafe { ffi::gsl_monte_miser_init(self.s) }
    }

    /// This routines uses the MISER Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined by
    /// the lower and upper limits in the arrays xl and xu, each of size dim. The integration uses a fixed number of function calls calls,
    /// and obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied. The result
    /// of the integration is returned in result, with an estimated absolute error abserr.
    ///
    /// It returns either Ok((result, abserr)) or Err(enums::Value).
    pub fn integrate<F: Fn(&[f64]) -> f64 + 'static>(&self, dim: usize, f: F, xl: &[f64], xu: &[f64],
                                                     t_calls: usize, r: &::Rng) -> Result<(f64, f64), ::Value> {
        unsafe {
            assert!(xl.len() == xu.len());
            let mut result = 0f64;
            let mut abserr = 0f64;
            let f: Box<Box<Fn(&[f64]) -> f64 + 'static>> = Box::new(Box::new(f));
            let mut func = ffi::gsl_monte_function {
                               f: transmute(monte_trampoline as usize),
                               dim: dim,
                               params: Box::into_raw(f) as *mut _,
                           };
            let ret = ffi::gsl_monte_miser_integrate(&mut func as *mut _ as *mut c_void,
                                                     xl.as_ptr(), xu.as_ptr(), xl.len(), t_calls,
                                                     ffi::FFI::unwrap(r), self.s,
                                                     (&mut result) as *mut c_double,
                                                     (&mut abserr) as *mut c_double);

            if ret == ::Value::Success {
                Ok((result, abserr))
            } else {
                Err(ret)
            }
        }
    }

    /// This function copies the parameters of the integrator state into the user-supplied params structure.
    pub fn get_params(&self) -> MiserParams {
        let mut m = MiserParams {
            estimate_frac: 0f64,
            min_calls: 0,
            min_calls_per_bisection: 0,
            alpha: 0f64,
            dither: 0f64,
        };

        unsafe {
            ffi::gsl_monte_miser_params_get(self.s, &mut m as *mut MiserParams);
        }
        m
    }

    /// This function sets the integrator parameters based on values provided in the params structure.
    pub fn set_params(&self, params: &MiserParams) {
        unsafe {
            ffi::gsl_monte_miser_params_set(self.s, params as *const MiserParams);
        }
    }
}

impl Drop for MiserMonteCarlo {
    fn drop(&mut self) {
        unsafe { ffi::gsl_monte_miser_free(self.s) };
        self.s = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_monte_miser_state> for MiserMonteCarlo {
    fn wrap(s: *mut ffi::gsl_monte_miser_state) -> MiserMonteCarlo {
        MiserMonteCarlo {
            s: s
        }
    }

    fn unwrap(s: &MiserMonteCarlo) -> *mut ffi::gsl_monte_miser_state {
        s.s
    }
}

#[repr(C)]
pub struct MiserParams {
    /// This parameter specifies the fraction of the currently available number of function calls which
    /// are allocated to estimating the variance at each recursive step. The default value is 0.1.
    pub estimate_frac: f64,
    /// This parameter specifies the minimum number of function calls required for each estimate of the
    /// variance. If the number of function calls allocated to the estimate using estimate_frac falls
    /// below min_calls then min_calls are used instead. This ensures that each estimate maintains a
    /// reasonable level of accuracy. The default value of min_calls is 16 * dim.
    pub min_calls: usize,
    /// This parameter specifies the minimum number of function calls required to proceed with a bisection
    /// step. When a recursive step has fewer calls available than min_calls_per_bisection it performs
    /// a plain Monte Carlo estimate of the current sub-region and terminates its branch of the recursion.
    /// The default value of this parameter is 32 * min_calls.
    pub min_calls_per_bisection: usize,
    /// This parameter controls how the estimated variances for the two sub-regions of a bisection are
    /// combined when allocating points. With recursive sampling the overall variance should scale better
    /// than 1/N, since the values from the sub-regions will be obtained using a procedure which explicitly
    /// minimizes their variance. To accommodate this behavior the MISER algorithm allows the total variance
    /// to depend on a scaling parameter \alpha,
    ///
    /// \Var(f) = {\sigma_a \over N_a^\alpha} + {\sigma_b \over N_b^\alpha}.
    /// 
    /// The authors of the original paper describing MISER recommend the value \alpha = 2 as a good choice,
    /// obtained from numerical experiments, and this is used as the default value in this implementation.
    pub alpha: f64,
    /// This parameter introduces a random fractional variation of size dither into each bisection, which
    /// can be used to break the symmetry of integrands which are concentrated near the exact center of
    /// the hypercubic integration region. The default value of dither is zero, so no variation is introduced.
    /// If needed, a typical value of dither is 0.1.
    pub dither: f64,
}

/// The VEGAS algorithm of Lepage is based on importance sampling. It samples points from the probability
/// distribution described by the function |f|, so that the points are concentrated in the regions that
/// make the largest contribution to the integral.
///
/// In general, if the Monte Carlo integral of f is sampled with points distributed according to a
/// probability distribution described by the function g, we obtain an estimate E_g(f; N),
///
/// E_g(f; N) = E(f/g; N)
///
/// with a corresponding variance,
///
/// \Var_g(f; N) = \Var(f/g; N).
///
/// If the probability distribution is chosen as g = |f|/I(|f|) then it can be shown that the variance
/// V_g(f; N) vanishes, and the error in the estimate will be zero. In practice it is not possible to
/// sample from the exact distribution g for an arbitrary function, so importance sampling algorithms
/// aim to produce efficient approximations to the desired distribution.
///
/// The VEGAS algorithm approximates the exact distribution by making a number of passes over the
/// integration region while histogramming the function f. Each histogram is used to define a sampling
/// distribution for the next pass. Asymptotically this procedure converges to the desired distribution.
/// In order to avoid the number of histogram bins growing like K^d the probability distribution is
/// approximated by a separable function: g(x_1, x_2, ...) = g_1(x_1) g_2(x_2) ... so that the number
/// of bins required is only Kd. This is equivalent to locating the peaks of the function from the
/// projections of the integrand onto the coordinate axes. The efficiency of VEGAS depends on the
/// validity of this assumption. It is most efficient when the peaks of the integrand are well-localized.
/// If an integrand can be rewritten in a form which is approximately separable this will increase
/// the efficiency of integration with VEGAS.
///
/// VEGAS incorporates a number of additional features, and combines both stratified sampling and
/// importance sampling. The integration region is divided into a number of “boxes”, with each box
/// getting a fixed number of points (the goal is 2). Each box can then have a fractional number of
/// bins, but if the ratio of bins-per-box is less than two, Vegas switches to a kind variance reduction
/// (rather than importance sampling).
///
/// The VEGAS algorithm computes a number of independent estimates of the integral internally, according
/// to the iterations parameter described below, and returns their weighted average. Random sampling of
/// the integrand can occasionally produce an estimate where the error is zero, particularly if the function
/// is constant in some regions. An estimate with zero error causes the weighted average to break down and
/// must be handled separately. In the original Fortran implementations of VEGAS the error estimate is made
/// non-zero by substituting a small value (typically 1e-30). The implementation in GSL differs from this
/// and avoids the use of an arbitrary constant—it either assigns the value a weight which is the average
/// weight of the preceding estimates or discards it according to the following procedure,
///
/// current estimate has zero error, weighted average has finite error
///
///     The current estimate is assigned a weight which is the average weight of the preceding estimates.
/// current estimate has finite error, previous estimates had zero error
///
///     The previous estimates are discarded and the weighted averaging procedure begins with the current estimate.
/// current estimate has zero error, previous estimates had zero error
///
///     The estimates are averaged using the arithmetic mean, but no error is computed.
pub struct VegasMonteCarlo {
    s: *mut ffi::gsl_monte_vegas_state,
}

impl VegasMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions.
    /// The workspace is used to maintain the state of the integration.
    pub fn new(dim: usize) -> Option<VegasMonteCarlo> {
        let tmp_pointer = unsafe { ffi::gsl_monte_vegas_alloc(dim) };

        if tmp_pointer.is_null() {
            None
        } else {
            Some(VegasMonteCarlo {
                s: tmp_pointer
            })
        }
    }

    /// This function initializes a previously allocated integration state. This allows an existing workspace
    /// to be reused for different integrations.
    pub fn init(&self) -> ::Value {
        unsafe { ffi::gsl_monte_vegas_init(self.s) }
    }

    /// This routines uses the VEGAS Monte Carlo algorithm to integrate the function f over the dim-dimensional
    /// hypercubic region defined by the lower and upper limits in the arrays xl and xu, each of size dim.
    /// The integration uses a fixed number of function calls calls, and obtains random sampling points using
    /// the random number generator r. A previously allocated workspace s must be supplied. The result of the
    /// integration is returned in result, with an estimated absolute error abserr. The result and its error
    /// estimate are based on a weighted average of independent samples. The chi-squared per degree of freedom
    /// for the weighted average is returned via the state struct component, s->chisq, and must be consistent
    /// with 1 for the weighted average to be reliable.
    ///
    /// It returns either Ok((result, abserr)) or Err(enums::Value).
    pub fn integrate<F: Fn(&[f64]) -> f64 + 'static>(&self, dim: usize, f: F, xl: &[f64], xu: &[f64],
                                                     t_calls: usize, r: &::Rng) -> Result<(f64, f64), ::Value> {
        unsafe {
            assert!(xl.len() == xu.len());
            let mut result = 0f64;
            let mut abserr = 0f64;
            let f: Box<Box<Fn(&[f64]) -> f64 + 'static>> = Box::new(Box::new(f));
            let mut func = ffi::gsl_monte_function {
                               f: transmute(monte_trampoline as usize),
                               dim: dim,
                               params: Box::into_raw(f) as *mut _,
                           };
            let ret = ffi::gsl_monte_vegas_integrate(&mut func as *mut _ as *mut c_void,
                                                     xl.as_ptr(), xu.as_ptr(), xl.len(), t_calls,
                                                     ffi::FFI::unwrap(r), self.s,
                                                     (&mut result) as *mut c_double,
                                                     (&mut abserr) as *mut c_double);

            if ret == ::Value::Success {
                Ok((result, abserr))
            } else {
                Err(ret)
            }
        }
    }

    /// This function returns the chi-squared per degree of freedom for the weighted estimate of the integral.
    /// The returned value should be close to 1. A value which differs significantly from 1 indicates that
    /// the values from different iterations are inconsistent. In this case the weighted error will be
    /// under-estimated, and further iterations of the algorithm are needed to obtain reliable results.
    pub fn chisq(&self) -> f64 {
        unsafe {
            ffi::gsl_monte_vegas_chisq(self.s)
        }
    }

    /// This function returns the raw (unaveraged) values of the integral result and its error sigma from
    /// the most recent iteration of the algorithm.
    pub fn runval(&self, result: &mut f64, sigma: &mut f64) {
        unsafe {
            ffi::gsl_monte_vegas_runval(self.s, result as *mut c_double, sigma as *mut c_double)
        }
    }

    /*pub params_get(&self) -> VegasParams {
        ;
    }

    pub params_set(&self, params: &VegasParams) {
        ;
    }*/
}

impl Drop for VegasMonteCarlo {
    fn drop(&mut self) {
        unsafe { ffi::gsl_monte_vegas_free(self.s) };
        self.s = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_monte_vegas_state> for VegasMonteCarlo {
    fn wrap(s: *mut ffi::gsl_monte_vegas_state) -> VegasMonteCarlo {
        VegasMonteCarlo {
            s: s
        }
    }

    fn unwrap(s: &VegasMonteCarlo) -> *mut ffi::gsl_monte_vegas_state {
        s.s
    }
}

unsafe extern "C" fn monte_trampoline(x: *mut c_double, dim: size_t, param: *mut c_void) -> c_double {
    let f: &Box<Fn(&[f64]) -> f64 + 'static> = transmute(param);
    f(slice::from_raw_parts(x, dim as usize))
}
