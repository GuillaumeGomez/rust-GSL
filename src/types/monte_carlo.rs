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

##MISER

The MISER algorithm of Press and Farrar is based on recursive stratified sampling. This technique aims to reduce the overall integration error 
by concentrating integration points in the regions of highest variance.

The idea of stratified sampling begins with the observation that for two disjoint regions a and b with Monte Carlo estimates of the integral 
E_a(f) and E_b(f) and variances \sigma_a^2(f) and \sigma_b^2(f), the variance \Var(f) of the combined estimate E(f) = (1/2) (E_a(f) + E_b(f)) 
is given by,

\Var(f) = (\sigma_a^2(f) / 4 N_a) + (\sigma_b^2(f) / 4 N_b).
It can be shown that this variance is minimized by distributing the points such that,

N_a / (N_a + N_b) = \sigma_a / (\sigma_a + \sigma_b).
Hence the smallest error estimate is obtained by allocating sample points in proportion to the standard deviation of the function in each 
sub-region.

The MISER algorithm proceeds by bisecting the integration region along one coordinate axis to give two sub-regions at each step. The direction 
is chosen by examining all d possible bisections and selecting the one which will minimize the combined variance of the two sub-regions. The 
variance in the sub-regions is estimated by sampling with a fraction of the total number of points available to the current step. The same 
procedure is then repeated recursively for each of the two half-spaces from the best bisection. The remaining sample points are allocated to 
the sub-regions using the formula for N_a and N_b. This recursive allocation of integration points continues down to a user-specified depth 
where each sub-region is integrated using a plain Monte Carlo estimate. These individual values and their error estimates are then combined 
upwards to give an overall result and an estimate of its error.

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
use std::f64::INFINITY;
use std::intrinsics::{sqrtf64, powf64, fabsf64, floorf64, logf64};
use std::default::Default;
use libc::c_void;
use std::num::Int;
use c_str::ToCStr;
use c_vec::CVec;

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
    pub fn init(&self) -> ::Value {
        unsafe { ffi::gsl_monte_plain_init(self.s) }
    }

    /// This routines uses the plain Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined
    /// by the lower and upper limits in the arrays xl and xu, each of the same size. The integration uses a fixed number of function calls
    /// calls, and obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied.
    /// The result of the integration is returned in result, with an estimated absolute error abserr.
    pub fn integrate<T>(&self, f: ::monte_function<T>, arg: &mut T, xl: &[f64], xu: &[f64], calls: u64, r: &::Rng, result: &mut f64,
        abserr: &mut f64) -> ::Value {
        unsafe {
            let mut m = 0f64;
            let mut q = 0f64;
            let dim = (*self.s).dim as usize;
            let mut t_x = CVec::new((*self.s).x, dim);
            let x = t_x.as_mut_slice();

            if xl.len() != dim || xu.len() != dim {
                rgsl_error!("number of dimensions must match allocated size", ::Value::Inval);
            }

            for i in range(0us, dim as usize) {
                if xu[i] <= xl[i] {
                    rgsl_error!("xu must be greater than xl", ::Value::Inval);
                }

                if xu[i] - xl[i] > ::DBL_MAX {
                    rgsl_error!("Range of integration is too large, please rescale", ::Value::Inval);
                }
            }

            /* Compute the volume of the region */
            let mut vol = 1f64;

            for i in range(0us, dim) {
                vol *= xu[i] - xl[i];
            }

            for n in range(0us, calls as usize) {
                /* Choose a random point in the integration region */
                for i in range(0us, dim) {
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

            ::Value::Success
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

pub struct MiserMonteCarlo {
    s: *mut ffi::gsl_monte_miser_state
}

impl MiserMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions. The workspace is used to maintain
    /// the state of the integration.
    pub fn new(dim: u64) -> Option<MiserMonteCarlo> {
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
    #[allow(unused_assignments)]
    pub fn integrate<T>(&self, f: ::monte_function<T>, arg: &mut T, xl: &[f64], xu: &[f64], t_calls: u64, r: &::Rng, result: &mut f64,
        abserr: &mut f64) -> ::Value {
        unsafe {
            let mut calls = t_calls;
            let mut calls_l = 0u64;
            let mut calls_r = 0u64;
            let min_calls = (*self.s).min_calls;
            let mut i_bisect = 0us;
            let dim = (*self.s).dim as usize;
            let mut found_best = false;

            let mut res_est = 0f64;
            let mut err_est = 0f64;
            let mut res_r = 0f64;
            let mut err_r = 0f64;
            let mut res_l = 0f64;
            let mut err_l = 0f64;

            let mut weight_l = 0f64;
            let mut weight_r = 0f64;

            let mut t_x = CVec::new((*self.s).x, dim);
            let x = t_x.as_mut_slice();
            let mut t_xmid = CVec::new((*self.s).xmid, dim);
            let xmid = t_xmid.as_mut_slice();
            let mut t_sigma_l = CVec::new((*self.s).sigma_l, dim);
            let sigma_l = t_sigma_l.as_mut_slice();
            let mut t_sigma_r = CVec::new((*self.s).sigma_r, dim);
            let sigma_r = t_sigma_r.as_mut_slice();

            if dim != xl.len() || dim != xu.len() {
                rgsl_error!("number of dimensions must match allocated size", ::Value::Inval);
            }

            for i in range(0us, dim) {
                if xu[i] <= xl[i] {
                    rgsl_error!("xu must be greater than xl", ::Value::Inval);
                }

                if xu[i] - xl[i] > ::DBL_MAX {
                    rgsl_error!("Range of integration is too large, please rescale", ::Value::Inval);
                }
            }

            if (*self.s).alpha < 0f64 {
                rgsl_error!("alpha must be non-negative", ::Value::Inval);
            }

            /* Compute volume */
            let mut vol = 1f64;

            for i in range(0us, dim) {
                vol *= xu[i] - xl[i];
            }

            if calls < (*self.s).min_calls_per_bisection {
                let mut m = 0f64;
                let mut q = 0f64;

                if calls < 2 {
                    rgsl_error!("insufficient calls for subvolume", ::Value::Failed);
                }

                for n in range(0us, calls as usize) {
                    /* Choose a random point in the integration region */
                    for i in range(0us, dim) {
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

                *abserr = vol * sqrtf64(q / (calls as f64 * (calls as f64 - 1f64)));

                return ::Value::Success;
            }

            let tmp = calls as f64 * (*self.s).estimate_frac;
            let estimate_calls = if min_calls as f64 > tmp {
                min_calls as u64
            } else {
                tmp as u64
            };

            if estimate_calls < dim as u64 * 4u64 {
                rgsl_error!("insufficient calls to sample all halfspaces", ::Value::Sanity);
            }

            /* Flip coins to bisect the integration region with some fuzz */
            for i in range(0us, dim) {
                let s = if (r.uniform() - 0.5f64) >= 0f64 { (*self.s).dither } else { -(*self.s).dither };
                
                xmid[i] = (0.5f64 + s) * xl[i] + (0.5f64 - s) * xu[i];
            }

            /* The idea is to chose the direction to bisect based on which will
             give the smallest total variance.  We could (and may do so later)
             use MC to compute these variances.  But the NR guys simply estimate
             the variances by finding the min and max function values 
             for each half-region for each bisection. */
            estimate_corrmc(f, arg, xl, xu, estimate_calls, r, self, &mut res_est, &mut err_est, xmid, sigma_l, sigma_r);

            /* We have now used up some calls for the estimation */

            calls -= estimate_calls;

            /* Now find direction with the smallest total "variance" */

            {
                let mut best_var = ::DBL_MAX;
                let beta = 2f64 / (1f64 + (*self.s).alpha);
                weight_r = 1f64;
                weight_l = weight_r;

                for i in range(0us, dim) {
                    if sigma_l[i] >= 0f64 && sigma_r[i] >= 0f64 {
                        /* estimates are okay */
                        let var = powf64(sigma_l[i], beta) + powf64(sigma_r[i], beta);

                        if var <= best_var {
                            found_best = true;
                            best_var = var;
                            i_bisect = i;
                            weight_l = powf64(sigma_l[i], beta);
                            weight_r = powf64(sigma_r[i], beta);

                            if weight_l == 0f64 && weight_r == 0f64 {
                                weight_l = 1f64;
                                weight_r = 1f64;
                            }
                        }
                    } else {
                        if sigma_l[i] < 0f64 {
                            rgsl_error!("no points in left-half space!", ::Value::Sanity);
                        }
                        if sigma_r[i] < 0f64 {
                            rgsl_error!("no points in right-half space!", ::Value::Sanity);
                        }
                    }
                }
            }

            if !found_best {
                /* All estimates were the same, so chose a direction at random */
                i_bisect = r.uniform_int(dim as u64) as usize;
            }

            let xbi_l = xl[i_bisect];
            let xbi_m = xmid[i_bisect];
            let xbi_r = xu[i_bisect];

            /* Get the actual fractional sizes of the two "halves", and
             distribute the remaining calls among them */
            {
                let fraction_l = fabsf64((xbi_m - xbi_l) / (xbi_r - xbi_l));
                let fraction_r = 1f64 - fraction_l;

                let a = fraction_l * weight_l;
                let b = fraction_r * weight_r;

                calls_l = (min_calls as f64 + (calls as f64 - 2f64 * min_calls as f64) * (a as f64 / (a as f64 + b as f64))) as u64;
                calls_r = (min_calls as f64 + (calls as f64 - 2f64 * min_calls as f64) * (b as f64 / (a as f64 + b as f64))) as u64;
            }

            /* Compute the integral for the left hand side of the bisection */

            /* Due to the recursive nature of the algorithm we must allocate
             some new memory for each recursive call */

            {
                let mut xu_tmp : Vec<f64> = Vec::with_capacity(dim);

                // Useless in Rust...
                /*if xu_tmp == 0 {
                    rgsl_error!("out of memory for left workspace", ::NoMem);
                }*/

                for it in xu.iter() {
                    xu_tmp.push(*it);
                }

                xu_tmp[i_bisect] = xbi_m;

                let status = self.integrate(f, arg, xl.as_slice(), xu_tmp.as_slice(), calls_l, r, &mut res_l, &mut err_l);

                if status != ::Value::Success {
                    return status;
                }
            }

            /* Compute the integral for the right hand side of the bisection */

            {
                let mut xl_tmp : Vec<f64> = Vec::with_capacity(dim);

                // Useless in Rust...
                /*if xl_tmp == 0 {
                    rgsl_error!("out of memory for right workspace", ::NoMem);
                }*/

                for it in xl.iter() {
                    xl_tmp.push(*it);
                }

                xl_tmp[i_bisect] = xbi_m;

                let status = self.integrate(f, arg, xl_tmp.as_slice(), xu, calls_r, r, &mut res_r, &mut err_r);

                if status != ::Value::Success {
                    return status;
                }
            }

            *result = res_l + res_r;
            *abserr = sqrtf64(err_l * err_l + err_r * err_r);

            ::Value::Success
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

fn estimate_corrmc<T>(f: ::monte_function<T>, arg: &mut T, xl: &[f64], xu: &[f64], calls: u64, r: &::Rng, state: &MiserMonteCarlo,
    result: &mut f64, abserr: &mut f64, xmid: &[f64], sigma_l: &mut [f64], sigma_r: &mut [f64]) -> ::Value {
    unsafe {
        let dim = (*state.s).dim as usize;
        let mut t_x = CVec::new((*state.s).x, dim);
        let x = t_x.as_mut_slice();
        let mut t_fsum_l = CVec::new((*state.s).fsum_l, dim);
        let fsum_l = t_fsum_l.as_mut_slice();
        let mut t_fsum_r = CVec::new((*state.s).fsum_r, dim);
        let fsum_r = t_fsum_r.as_mut_slice();
        let mut t_fsum2_l = CVec::new((*state.s).fsum2_l, dim);
        let fsum2_l = t_fsum2_l.as_mut_slice();
        let mut t_fsum2_r = CVec::new((*state.s).fsum2_r, dim);
        let fsum2_r = t_fsum2_r.as_mut_slice();
        let mut t_hits_l = CVec::new((*state.s).hits_l, dim);
        let hits_l = t_hits_l.as_mut_slice();
        let mut t_hits_r = CVec::new((*state.s).hits_r, dim);
        let hits_r = t_hits_r.as_mut_slice();

        let mut m = 0f64;
        let mut q = 0f64; 
        let mut vol = 1f64;

        for i in range(0us, dim) {
            vol *= xu[i] - xl[i];
            hits_r[i] = 0;
            hits_l[i] = hits_r[i];
            fsum_r[i] = 0f64;
            fsum_l[i] = fsum_r[i];
            fsum2_r[i] = 0f64;
            fsum2_l[i] = fsum2_r[i];
            sigma_r[i] = -1f64;
            sigma_l[i] = sigma_r[i];
        }

        for n in range(0us, calls as usize) {
            let j = (n / 2us) % dim;
            let side = n % 2us;

            for i in range(0us, dim) {
                let z = r.uniform_pos();

                if i != j {
                    x[i] = xl[i] + z * (xu[i] - xl[i]);
                } else {
                    if side == 0 {
                        x[i] = xmid[i] + z * (xu[i] - xmid[i]);
                    } else {
                        x[i] = xl[i] + z * (xmid[i] - xl[i]);
                    }
                }
            }

            let fval = f(x, arg);

            /* recurrence for mean and variance */
            {
                let d = fval - m;

                m += d / (n as f64 + 1f64);
                q += d * d * (n as f64 / (n as f64 + 1f64));
            }

            /* compute the variances on each side of the bisection */
            for i in range(0us, dim) {
                if x[i] <= xmid[i] {
                    fsum_l[i] += fval;
                    fsum2_l[i] += fval * fval;
                    hits_l[i] += 1;
                } else {
                    fsum_r[i] += fval;
                    fsum2_r[i] += fval * fval;
                    hits_r[i] += 1;
                }
            }
        }

        for i in range(0us, dim) {
            let fraction_l = (xmid[i] - xl[i]) / (xu[i] - xl[i]);

            if hits_l[i] > 0 {
                fsum_l[i] /= hits_l[i] as f64;
                sigma_l[i] = sqrtf64(fsum2_l[i] - fsum_l[i] * fsum_l[i] / hits_l[i] as f64);
                sigma_l[i] *= fraction_l * vol / hits_l[i] as f64;
            }

            if hits_r[i] > 0 {
                fsum_r[i] /= hits_r[i] as f64;
                sigma_r[i] = sqrtf64(fsum2_r[i] - fsum_r[i] * fsum_r[i] / hits_r[i] as f64);
                sigma_r[i] *= (1f64 - fraction_l) * vol / hits_r[i] as f64;
            }
        }

        *result = vol * m;

        if calls < 2 {
            *abserr = INFINITY;
        } else {
            *abserr = vol * sqrtf64(q / (calls as f64 * (calls as f64 - 1f64)));
        }

        ::Value::Success
    }
}

#[repr(C)]
#[derive(Copy)]
pub struct VegasParams {
    /// The parameter alpha controls the stiffness of the rebinning algorithm. It is typically set between one and two. A value of zero prevents
    /// rebinning of the grid. The default value is 1.5.
    pub alpha: f64,
    /// The number of iterations to perform for each call to the routine. The default value is 5 iterations.
    pub iterations: u64,
    /// Setting this determines the stage of the calculation. Normally, stage = 0 which begins with a new uniform grid and empty weighted average.
    /// Calling VEGAS with stage = 1 retains the grid from the previous run but discards the weighted average, so that one can “tune” the grid
    /// using a relatively small number of points and then do a large run with stage = 1 on the optimized grid. Setting stage = 2 keeps the grid
    /// and the weighted average from the previous run, but may increase (or decrease) the number of histogram bins in the grid depending on the
    /// number of calls available. Choosing stage = 3 enters at the main loop, so that nothing is changed, and is equivalent to performing
    /// additional iterations in a previous call.
    pub stage: i32,
    /// The possible choices are GSL_VEGAS_MODE_IMPORTANCE, GSL_VEGAS_MODE_STRATIFIED, GSL_VEGAS_MODE_IMPORTANCE_ONLY. This determines whether
    /// VEGAS will use importance sampling or stratified sampling, or whether it can pick on its own. In low dimensions VEGAS uses strict
    /// stratified sampling (more precisely, stratified sampling is chosen if there are fewer than 2 bins per box).
    pub mode: ::VegasMode,
    /// These parameters set the level of information printed by VEGAS. All information is written to the stream ostream. The default setting
    /// of verbose is -1, which turns off all output. A verbose value of 0 prints summary information about the weighted average and final
    /// result, while a value of 1 also displays the grid coordinates. A value of 2 prints information from the rebinning procedure for each
    /// iteration.
    pub verbose: i32,
    ostream: *mut c_void
}

impl Default for VegasParams {
    fn default() -> VegasParams {
        VegasParams {
            alpha: 1.5f64,
            iterations: 5u64,
            stage: 0i32,
            mode: ::VegasMode::Importance,
            verbose: 0i32,
            ostream: ::std::ptr::null_mut()
        }
    }
}

pub struct VegasMonteCarlo {
    s: *mut ffi::gsl_monte_vegas_state
}

impl VegasMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions. The workspace is used to maintain
    /// the state of the integration.
    pub fn new(dim: u64) -> Option<VegasMonteCarlo> {
        let tmp_pointer = unsafe { ffi::gsl_monte_vegas_alloc(dim) };

        if tmp_pointer.is_null() {
            None
        } else {
            Some(VegasMonteCarlo {
                s: tmp_pointer
            })
        }
    }

    /// This function initializes a previously allocated integration state. This allows an existing workspace to be reused for different integrations.
    pub fn init(&self) -> ::Value {
        unsafe { ffi::gsl_monte_vegas_init(self.s) }
    }

    /// This routines uses the VEGAS Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined by
    /// the lower and upper limits in the arrays xl and xu, each of size dim. The integration uses a fixed number of function calls calls, and
    /// obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied. The result of
    /// the integration is returned in result, with an estimated absolute error abserr. The result and its error estimate are based on a weighted
    /// average of independent samples. The chi-squared per degree of freedom for the weighted average is returned via the state struct
    /// component, s->chisq, and must be consistent with 1 for the weighted average to be reliable.
    pub fn integrate<T>(&self, f: ::monte_function<T>, arg: &mut T, xl: &[f64], xu: &[f64], t_calls: u64, r: &::Rng, result: &mut f64,
        abserr: &mut f64) -> ::Value {
        unsafe {
            let mut calls = t_calls;
            let dim = (*self.s).dim as usize;

            if xl.len() != dim || xu.len() != dim {
                rgsl_error!("number of dimensions must match allocated size", ::Value::Inval);
            }

            for i in range(0us, dim) {
                if xu[i] <= xl[i] {
                    rgsl_error!("xu must be greater than xl", ::Value::Inval);
                }

                if xu[i] - xl[i] > ::DBL_MAX {
                    rgsl_error!("Range of integration is too large, please rescale", ::Value::Inval);
                }
            }

            if (*self.s).stage == 0 {
                init_grid(self.s, xl, xu);

                // FIXME : allow verbose mode
                /*if (*self.s).verbose >= 0 {
                    print_lim(self.s, xl, xu, dim);
                }*/
            }

            if (*self.s).stage <= 1 {
                (*self.s).wtd_int_sum = 0f64;
                (*self.s).sum_wgts = 0f64;
                (*self.s).chi_sum = 0f64;
                (*self.s).it_num = 1;
                (*self.s).samples = 0;
                (*self.s).chisq = 0f64;
            }

            if (*self.s).stage <= 2 {
                let mut bins = (*self.s).bins_max;
                let mut boxes = 1u64;

                if (*self.s).mode != ::VegasMode::ImportanceOnly {
                    /* shooting for 2 calls/box */
                    boxes = floorf64(powf64(calls as f64 / 2f64, 1f64 / dim as f64)) as u64;
                    (*self.s).mode = ::VegasMode::Importance;

                    if 2 * boxes >= (*self.s).bins_max {
                        /* if bins/box < 2 */
                        let tmp = boxes / (*self.s).bins_max;
                        let box_per_bin = if tmp > 1 { tmp } else { 1 };

                        let tmp2 = boxes / box_per_bin;
                        bins = if tmp2 < (*self.s).bins_max { tmp2 } else { (*self.s).bins_max };
                        boxes = box_per_bin * bins;

                        (*self.s).mode = ::VegasMode::Stratified;
                    }
                }

                {
                    let tot_boxes = boxes.pow(dim);
                    let tmp = calls / tot_boxes;

                    (*self.s).calls_per_box = if tmp > 2 { tmp as u32 } else { 2 };
                    calls = ((*self.s).calls_per_box * tot_boxes as u32) as u64;
                }

                /* total volume of x-space/(avg num of calls/bin) */
                (*self.s).jac = (*self.s).vol * powf64(bins as f64, dim as f64) / (calls as f64);

                (*self.s).boxes = boxes as u32;

                /* If the number of bins changes from the previous invocation, bins
                 are expanded or contracted accordingly, while preserving bin density */

                if bins != (*self.s).bins as u64 {
                    resize_grid(self.s, bins as usize);

                    // FIXME: allow verbose mode
                    /*if (*self.s).verbose > 1 {
                        print_grid(self.s, dim);
                    }*/
                }

                // FIXME: allow verbose mode
                /*if (*self.s).verbose >= 0 {
                    print_head(self.s, dim, calls, (*self.s).it_num, (*self.s).bins, (*self.s).boxes);
                }*/
            }

            (*self.s).it_start = (*self.s).it_num;

            let mut cum_int = 0f64;
            let mut cum_sig = 0f64;

            for it in range(0us, (*self.s).iterations as usize) {
                let mut intgrl = 0f64;
                let mut tss = 0f64;
                let calls_per_box = (*self.s).calls_per_box;
                let jacbin = (*self.s).jac;

                let mut x = CVec::new((*self.s).x, dim);
                let mut bin = CVec::new((*self.s).bin, dim);

                (*self.s).it_num = (*self.s).it_start + it as u32;

                reset_grid_values(self.s);
                let mut t_box = CVec::new((*self.s).box_, dim);
                init_box_coord(self.s, t_box.as_mut_slice());

                //let mut c = 0u;
                loop {
                    let mut m = 0f64;
                    let mut q = 0f64;

                    for k in range(0us, calls_per_box as usize) {
                        let mut bin_vol = 0f64;

                        random_point(x.as_mut_slice(), bin.as_mut_slice(), &mut bin_vol, t_box.as_mut_slice(), xl, xu, self.s, r/*, c == 4640*/);

                        let fval = jacbin * bin_vol * f(x.as_mut_slice(), arg);

                        /* recurrence for mean and variance (sum of squares) */
                        {
                            let d = fval - m;

                            m += d / (k as f64 + 1f64);
                            q += d * d * (k as f64 / (k as f64 + 1f64));
                        }

                        if (*self.s).mode != ::VegasMode::Stratified {
                            let f_sq = fval * fval;

                            accumulate_distribution(self.s, bin.as_mut_slice(), f_sq);
                        }

                    }

                    intgrl += m * calls_per_box as f64;

                    let f_sq_sum = q * calls_per_box as f64;

                    //c += 1;
                    tss += f_sq_sum;

                    if (*self.s).mode == ::VegasMode::Stratified {
                        accumulate_distribution(self.s, bin.as_mut_slice(), f_sq_sum);
                    }
                    if !change_box_coord(self.s, t_box.as_mut_slice()) {
                        break;
                    }
                }

                /* Compute final results for this iteration   */
                let var = tss / (calls_per_box as f64 - 1f64);

                let wgt = if var > 0f64 {
                    1f64 / var
                } else if (*self.s).sum_wgts > 0f64 {
                    (*self.s).sum_wgts / (*self.s).samples as f64
                } else {
                    0f64
                };
                
                let intgrl_sq = intgrl * intgrl;

                let sig = sqrtf64(var);

                (*self.s).result = intgrl;
                (*self.s).sigma  = sig;

                if wgt > 0f64 {
                    let sum_wgts = (*self.s).sum_wgts;
                    let wtd_int_sum = (*self.s).wtd_int_sum;
                    let m = if sum_wgts > 0f64 {
                        wtd_int_sum / sum_wgts
                    } else {
                        0f64
                    };
                    let q = intgrl - m;

                    (*self.s).samples += 1;
                    (*self.s).sum_wgts += wgt;
                    (*self.s).wtd_int_sum += intgrl * wgt;
                    (*self.s).chi_sum += intgrl_sq * wgt;

                    cum_int = (*self.s).wtd_int_sum / (*self.s).sum_wgts;
                    cum_sig = sqrtf64(1f64 / (*self.s).sum_wgts);

            /*#if USE_ORIGINAL_CHISQ_FORMULA
            // This is the chisq formula from the original Lepage paper.  It
            // computes the variance from <x^2> - <x>^2 and can suffer from
            // catastrophic cancellations, e.g. returning negative chisq.
                 if (state->samples > 1)
                   {
                     state->chisq = (state->chi_sum - state->wtd_int_sum * cum_int) /
                       (state->samples - 1.0);
                   }
            #else*/
            // The new formula below computes exactly the same quantity as above but using a stable recurrence
                    if (*self.s).samples == 1 {
                        (*self.s).chisq = 0f64;
                    } else {
                        (*self.s).chisq *= (*self.s).samples as f64 - 2f64;
                        (*self.s).chisq += (wgt / (1f64 + (wgt / sum_wgts))) * q * q;
                        (*self.s).chisq /= (*self.s).samples as f64 - 1f64;
                    }
            //#endif
                } else {
                    cum_int += (intgrl - cum_int) / (it as f64 + 1f64);
                    cum_sig = 0f64;
                }         


                // FIXME: allow verbose mode
                /*if (*self.s).verbose >= 0 {
                    print_res(self.s, (*self.s).it_num, intgrl, sig, cum_int, cum_sig, (*self.s).chisq);
                    if it + 1 == (*self.s).iterations && (*self.s).verbose > 0 {
                        print_grid(self.s, dim);
                    }
                }*/

                // FIXME: allow verbose mode
                /*if (*self.s).verbose > 1 {
                    print_dist(self.s, dim);
                }*/

                refine_grid(self.s);

                // FIXME: allow verbose mode
                /*if (*self.s).verbose > 1 {
                    print_grid(self.s, dim);
                }*/
            }

            /* By setting stage to 1 further calls will generate independent
             estimates based on the same grid, although it may be rebinned. */

            (*self.s).stage = 1;  

            *result = cum_int;
            *abserr = cum_sig;

            ::Value::Success
        }
    }

    /// This function returns the chi-squared per degree of freedom for the weighted estimate of the integral. The returned value should be
    /// close to 1. A value which differs significantly from 1 indicates that the values from different iterations are inconsistent. In this
    /// case the weighted error will be under-estimated, and further iterations of the algorithm are needed to obtain reliable results.
    pub fn chisq(&self) -> f64 {
        unsafe { ffi::gsl_monte_vegas_chisq(self.s as *const ffi::gsl_monte_vegas_state) }
    }

    /// This function returns the raw (unaveraged) values of the integral result and its error sigma from the most recent iteration of the
    /// algorithm.
    pub fn runval(&self, result: &mut f64, sigma: &mut f64) {
        unsafe { ffi::gsl_monte_vegas_runval(self.s as *const ffi::gsl_monte_vegas_state, result, sigma) }
    }

    /// This function copies the parameters of the integrator state into the returned params structure.
    pub fn params_get(&self) -> VegasParams {
        let mut tmp : VegasParams = Default::default();

        unsafe { ffi::gsl_monte_vegas_params_get(self.s as *const ffi::gsl_monte_vegas_state, &mut tmp) };
        tmp
    }

    /// This function sets the integrator parameters based on values provided in the params structure.
    pub fn params_set(&self, params: &VegasParams) {
        unsafe { ffi::gsl_monte_vegas_params_set(self.s, params as *const VegasParams) }
    }
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

unsafe fn init_grid(s: *mut ffi::gsl_monte_vegas_state, xl: &[f64], xu: &[f64]) {
    let mut vol = 1f64;
    let mut t_delx = CVec::new((*s).delx, xl.len());
    let delx = t_delx.as_mut_slice();

    (*s).bins = 1;

    for j in range(0us, xl.len()) {
        let dx = xu[j] - xl[j];
        
        delx[j] = dx;
        vol *= dx;

        *(*s).xi.offset(j as isize) = 0f64;
        *(*s).xi.offset((*s).dim as isize + j as isize) = 1f64;
    }

    (*s).vol = vol;
}

unsafe fn reset_grid_values(s: *mut ffi::gsl_monte_vegas_state) {
    let dim = (*s).dim as isize;
    let bins = (*s).bins as isize;

    for i in range(0is, bins) {
        for j in range(0is, dim) {
            *(*s).d.offset(i * dim + j) = 0f64;
        }
    }
}

unsafe fn accumulate_distribution(s: *mut ffi::gsl_monte_vegas_state, bin: &[i32], y: f64) {
    let dim = (*s).dim as usize;

    for j in range(0is, dim as isize) {
        let i = bin[j as usize] as isize;

        *(*s).d.offset(i * dim as isize + j) += y;
    }
}

#[allow(unused_variables)]
#[allow(unused_assignments)]
unsafe fn random_point(x: &mut [f64], bin: &mut [i32], bin_vol: &mut f64, box_: &[i32], xl: &[f64], xu: &[f64], s: *mut ffi::gsl_monte_vegas_state,
    r: &::Rng) {
    /* Use the random number generator r to return a random position x
       in a given box. The value of bin gives the bin location of the
       random position (there may be several bins within a given box) */

    let mut vol = 1f64;

    let dim = (*s).dim as usize;
    let bins = (*s).bins as usize;
    let boxes = (*s).boxes as usize;
    let mut t_delx = CVec::new((*s).delx, xl.len());
    let delx = t_delx.as_mut_slice();

    for j in range(0us, dim) {
        /* box[j] + ran gives the position in the box units, while z
           is the position in bin units.  */

        let z = ((box_[j] as f64 + r.uniform_pos()) / boxes as f64) * bins as f64;
        let k = z as isize;

        let mut y = 0f64;
        let mut bin_width = 0f64;

        bin[j] = k as i32;

        if k == 0 {
            bin_width = *(*s).xi.offset(dim as isize + j as isize);
            y = z * bin_width;
        } else {
            // FIXME: problem is here -> xi hasn't the good values, when c == 4640 on second turn
            bin_width = *(*s).xi.offset((k + 1) * (*s).dim as isize + j as isize) - *(*s).xi.offset(k * (*s).dim as isize + j as isize);
            y = *(*s).xi.offset(k as isize * (*s).dim as isize + j as isize) + (z - k as f64) * bin_width;
        }

        x[j] = xl[j] + y * delx[j];

        vol *= bin_width;
    }

    *bin_vol = vol;
}

unsafe fn resize_grid(s: *mut ffi::gsl_monte_vegas_state, bins: usize) {
    let dim = (*s).dim as usize;

    /* weight is ratio of bin sizes */
    let pts_per_bin = (*s).bins as f64 / bins as f64;

    for j in range(0us, dim) {
        let mut xnew = 0f64;
        let mut dw = 0f64;
        let mut i = 1us;

        for k in range(1us, (*s).bins as usize + 1us) {
            dw += 1f64;
            let xold = xnew;
            xnew = *(*s).xi.offset(k as isize * dim as isize + j as isize);

            while dw > pts_per_bin {
                dw -= pts_per_bin;
                *(*s).xin.offset(i as isize) = xnew - (xnew - xold) * dw;
                i += 1;
            }
        }

        for k in range(1us, bins) {
            *(*s).xi.offset(k as isize * dim as isize + j as isize) = *(*s).xin.offset(k as isize);
        }

        *(*s).xi.offset(bins as isize * dim as isize + j as isize) = 1f64;
    }

    (*s).bins = bins as u32;
}

unsafe fn refine_grid(s: *mut ffi::gsl_monte_vegas_state) {
    let dim = (*s).dim as usize;
    let bins = (*s).bins as usize;

    for j in range(0us, dim) {
        let mut t_weight = CVec::new((*s).weight, bins);
        let weight = t_weight.as_mut_slice();

        let mut oldg = *(*s).d.offset(j as isize);
        let mut newg = *(*s).d.offset(dim as isize + j as isize);

        *(*s).d.offset(j as isize) = (oldg + newg) / 2f64;
        let mut grid_tot_j = *(*s).xi.offset(j as isize);

        /* This implements gs[i][j] = (gs[i-1][j]+gs[i][j]+gs[i+1][j])/3 */

        for i in range(1us, bins - 1us) {
            let rc = oldg + newg;

            oldg = newg;
            newg = *(*s).d.offset((i as isize + 1) * dim as isize + j as isize);
            *(*s).d.offset(i as isize * dim as isize + j as isize) = (rc + newg) / 3f64;
            grid_tot_j += *(*s).d.offset(i as isize * dim as isize + j as isize);
        }
        *(*s).d.offset((bins as isize - 1is) * dim as isize + j as isize) = (newg + oldg) / 2f64;

        grid_tot_j += *(*s).d.offset((bins as isize - 1is) * dim as isize + j as isize);

        let mut tot_weight = 0f64;

        for i in range(0us, bins) {
            weight[i] = 0f64;

            if *(*s).d.offset(i as isize * dim as isize + j as isize) > 0f64 {
                oldg = grid_tot_j / *(*s).d.offset(i as isize * dim as isize + j as isize);
                /* damped change */
                weight[i] = powf64(((oldg - 1f64) / oldg / logf64(oldg)), (*s).alpha);
            }

          tot_weight += weight[i];

/*#ifdef DEBUG
          printf("weight[%d] = %g\n", i, weight[i]);
#endif*/
        }

        {
            let pts_per_bin = tot_weight / bins as f64;

            let mut xnew = 0f64;
            let mut dw = 0f64;
            let mut i = 1us;

            for k in range(0us, bins) {
                dw += weight[k];
                let xold = xnew;
                xnew = *(*s).xi.offset((k as isize + 1is) * dim as isize + j as isize);

                while dw > pts_per_bin {
                    dw -= pts_per_bin;
                    *(*s).xin.offset(i as isize) = xnew - (xnew - xold) * dw / weight[k];
                    i += 1;
                }
            }

            for k in range(1us, bins) {
                *(*s).xi.offset(k as isize * dim as isize + j as isize) = *(*s).xin.offset(k as isize);
            }

            *(*s).xi.offset(bins as isize * dim as isize + j as isize) = 1f64;
        }
    }
}

/*unsafe fn print_lim(state: *mut ffi::gsl_monte_vegas_state, xl: &[f64], xu : &[f64]) {
    if *(state).ostream.is_null() {
        return;
    }
    let t_output_stream : Option<Writer> = ::std::mem::transmute((*state).ostream);
    let output_stream = t_output_stream.unwrap();

    write!(output_stream, "The limits of integration are:\n");
    for j in range(0us, xl.len()) {
        write!(output_stream, "\nxl[{}]={}    xu[{}]={}", j, xl[j], j, xu[j]);
    }
    write!(output_stream, "\n");
    output_stream.flush();
}

unsafe fn print_head(state: *mut ffi::gsl_monte_vegas_state, num_dim: u32, calls: u32, it_num: u32, bins: u32, boxes: u32) {
    let t_output_stream : Option<File> = ::std::mem::transmute((*state).ostream);
    let output_stream = match t_output_stream {
        Some(os) => os,
        None => return,
    };

    write!(output_stream, "\nnum_dim={}, calls={}, it_num={}, max_it_num={} ", num_dim, calls, it_num, (*state).iterations);
    write!(output_stream, "verb={}, alph={:.2},\nmode={}, bins={}, boxes={}\n", (*state).verbose, (*state).alpha, (*state).mode,
        bins, boxes);
    write!(output_stream, "\n       single.......iteration                   ");
    write!(output_stream, "accumulated......results   \n");

    write!(output_stream, "iteration     integral    sigma             integral   ");
    write!(output_stream, "      sigma     chi-sq/it\n\n");
    output_stream.flush();
}

unsafe fn print_res(state: *mut ffi::gsl_monte_vegas_state, itr: u32, res: f64, err: f64, cum_res: f64, cum_err: f64, chi_sq: f64) {
    let t_output_stream : Option<File> = ::std::mem::transmute((*state).ostream);
    let output_stream = match t_output_stream {
        Some(os) => os,
        None => return,
    };

    write!(output_stream, "{:4}        {:6.4}e {:10.2}e          {:6.4}e      {:8.2}e  {:10.2}e\n",
        itr, res, err, cum_res, cum_err, chi_sq);
    output_stream.flush();
}

unsafe fn print_dist(state: *mut ffi::gsl_monte_vegas_state, dim: u32) {
    let t_output_stream : Option<File> = ::std::mem::transmute((*state).ostream);
    let output_stream = match t_output_stream {
        Some(os) => os,
        None => return,
    };
    let p = (*state).verbose;
    
    if p < 1 {
        return;
    }

    for j in range(0i, dim as isize) {
        write!(output_stream, "\n axis {} \n", j);
        write!(output_stream, "      x   g\n");

        for i in range(0i, (*state).bins as isize) {
            write!(output_stream, "weight [{:11.2}e , {:11.2}e] = ",
                *(*state).xi.offset(i * (*state).dim as isize + j),
                *(*state).xi.offset((i + 1is) * (*state).dim as isize + j));
            write!(output_stream, " {:11.2}e\n", *(*state).d.offset(i * (*state).dim as isize + j));
        }
        write!(output_stream, "\n");
    }
    write!(output_stream, "\n");
    output_stream.flush();
}

unsafe fn print_grid(state: *mut ffi::gsl_monte_vegas_state, dim: u32) {
    let t_output_stream : Option<File> = ::std::mem::transmute((*state).ostream);
    let output_stream = match t_output_stream {
        Some(os) => os,
        None => return,
    };
    let p = (*state).verbose;

    if p < 1 {
        return;
    }

    for j in range(0us, dim as usize) {
        write!(output_stream, "\n axis {} \n", j);
        write!(output_stream, "      x   \n");

        for i in range(0us, (*state).bins as usize) {
            write!(output_stream, "{:11.2}e", *(*state).xi.offset(i as isize * (*state).dim as isize + j as isize));
            
            if i % 5u == 4u {
                write!(output_stream, "\n");
            }
        }
        write!(output_stream, "\n");
    }
    write!(output_stream, "\n");
    output_stream.flush();
}*/

unsafe fn init_box_coord(s: *mut ffi::gsl_monte_vegas_state, box_: &mut [ffi::coord]) {
    let dim = (*s).dim as usize;

    for i in range(0us, dim) {
        box_[i] = 0;
    }
}

/* change_box_coord steps through the box coord like
   {0,0}, {0, 1}, {0, 2}, {0, 3}, {1, 0}, {1, 1}, {1, 2}, ...
*/
unsafe fn change_box_coord(s: *mut ffi::gsl_monte_vegas_state, box_: &mut [ffi::coord]) -> bool {
    let mut j = (*s).dim as isize - 1;
    let ng = (*s).boxes;

    while j >= 0 {
        box_[j as usize] = (box_[j as usize] + 1) % ng as i32;

        if box_[j as usize] != 0 {
            return true;
        }

        j -= 1;
    }

    false
}