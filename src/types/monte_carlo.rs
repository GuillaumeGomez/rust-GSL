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
!*/

use ffi;
use enums;
use std::f64::INFINITY;
use std::intrinsics::{sqrtf64, powf64, fabsf64};
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
    pub fn init(&self) -> enums::Value {
        unsafe { ffi::gsl_monte_miser_init(self.s) }
    }

    /// This routines uses the MISER Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined by
    /// the lower and upper limits in the arrays xl and xu, each of size dim. The integration uses a fixed number of function calls calls,
    /// and obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied. The result
    /// of the integration is returned in result, with an estimated absolute error abserr.
    #[allow(dead_assignment)]
    pub fn integrate<T>(&self, f: ::monte_function<T>, arg: &mut T, xl: &[f64], xu: &[f64], calls: u64, r: &::Rng, result: &mut f64,
        abserr: &mut f64) -> enums::Value {
        unsafe {
            let mut calls_l = 0u64;
            let mut calls_r = 0u64;
            let min_calls = (*self.s).min_calls;
            let mut i_bisect = 0u;
            let dim = (*self.s).dim as uint;
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
                rgsl_error!("number of dimensions must match allocated size", enums::Inval);
            }

            for i in range(0u, dim) {
                if xu[i] <= xl[i] {
                    rgsl_error!("xu must be greater than xl", enums::Inval);
                }

                if xu[i] - xl[i] > ::DBL_MAX {
                    rgsl_error!("Range of integration is too large, please rescale", enums::Inval);
                }
            }

            if (*self.s).alpha < 0f64 {
                rgsl_error!("alpha must be non-negative", enums::Inval);
            }

            /* Compute volume */
            let mut vol = 1f64;

            for i in range(0u, dim) {
                vol *= xu[i] - xl[i];
            }

            if calls < (*self.s).min_calls_per_bisection {
                let mut m = 0f64;
                let mut q = 0f64;

                if calls < 2 {
                    rgsl_error!("insufficient calls for subvolume", enums::Failed);
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

                *abserr = vol * sqrtf64(q / (calls as f64 * (calls as f64 - 1f64)));

                return enums::Success;
            }

            let tmp = calls as f64 * (*self.s).estimate_frac;
            let estimate_calls = if min_calls as f64 > tmp {
                min_calls as u64
            } else {
                tmp as u64
            };

            if estimate_calls < dim as u64 * 4u64 {
                rgsl_error!("insufficient calls to sample all halfspaces", enums::Sanity);
            }

            /* Flip coins to bisect the integration region with some fuzz */
            for i in range(0u, dim) {
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

            let t_calls = calls - estimate_calls;

            /* Now find direction with the smallest total "variance" */

            {
                let mut best_var = ::DBL_MAX;
                let beta = 2f64 / (1f64 + (*self.s).alpha);
                weight_r = 1f64;
                weight_l = weight_r;

                for i in range(0u, dim) {
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
                            rgsl_error!("no points in left-half space!", enums::Sanity);
                        }
                        if sigma_r[i] < 0f64 {
                            rgsl_error!("no points in right-half space!", enums::Sanity);
                        }
                    }
                }
            }

            if !found_best {
                /* All estimates were the same, so chose a direction at random */
                i_bisect = r.uniform_int(dim as u64) as uint;
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

                calls_l = min_calls + (t_calls - 2 * min_calls) * (a / (a + b)) as u64;
                calls_r = min_calls + (t_calls - 2 * min_calls) * (b / (a + b)) as u64;
            }

            /* Compute the integral for the left hand side of the bisection */

            /* Due to the recursive nature of the algorithm we must allocate
             some new memory for each recursive call */

            {
                let mut t_xu_tmp : Vec<f64> = Vec::with_capacity(dim);
                let xu_tmp = t_xu_tmp.as_mut_slice();

                // Useless in Rust...
                /*if xu_tmp == 0 {
                    rgsl_error!("out of memory for left workspace", enums::NoMem);
                }*/

                for i in range(0u, dim) {
                    xu_tmp[i] = xu[i];
                }

                xu_tmp[i_bisect] = xbi_m;

                let status = self.integrate(f, arg, xl.as_slice(), xu_tmp.as_slice(), calls_l, r, &mut res_l, &mut err_l);

                if status != enums::Success {
                    return status;
                }
            }

            /* Compute the integral for the right hand side of the bisection */

            {
                let mut t_xl_tmp : Vec<f64> = Vec::with_capacity(dim);
                let xl_tmp = t_xl_tmp.as_mut_slice();

                // Useless in Rust...
                /*if xl_tmp == 0 {
                    rgsl_error!("out of memory for right workspace", enums::NoMem);
                }*/

                for i in range(0u, dim) {
                    xl_tmp[i] = xl[i];
                }

                xl_tmp[i_bisect] = xbi_m;

                let status = self.integrate(f, arg, xl_tmp.as_slice(), xu, calls_r, r, &mut res_r, &mut err_r);

                if status != enums::Success {
                    return status;
                }
            }

            *result = res_l + res_r;
            *abserr = sqrtf64(err_l * err_l + err_r * err_r);

            enums::Success
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
    result: &mut f64, abserr: &mut f64, xmid: &[f64], sigma_l: &mut [f64], sigma_r: &mut [f64]) -> enums::Value {
    unsafe {
        let dim = (*state.s).dim as uint;
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

        for i in range(0u, dim) {
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

        for n in range(0u, calls as uint) {
            let j = (n / 2u) % dim;
            let side = n % 2u;

            for i in range(0u, dim) {
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
            for i in range(0u, dim) {
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

        for i in range(0u, dim) {
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

        enums::Success
    }
}