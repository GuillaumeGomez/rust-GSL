//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Monte Carlo Integration

This chapter describes routines for multidimensional Monte Carlo integration. These include the traditional Monte Carlo method and adaptive
algorithms such as VEGAS and MISER which use importance sampling and stratified sampling techniques. Each algorithm computes an estimate of
a multidimensional definite integral of the form,

I = \int_xl^xu dx \int_yl^yu  dy ...  f(x, y, ...)

over a hypercubic region ((x_l,x_u), (y_l,y_u), ...) using a fixed number of function calls. The routines also provide a statistical estimate
of the error on the result. This error estimate should be taken as a guide rather than as a strict error bound—random sampling of the region
may not uncover all the important features of the function, resulting in an underestimate of the error.

## Interface

All of the Monte Carlo integration routines use the same general form of interface. There is an allocator to allocate memory for control
variables and workspace, a routine to initialize those control variables, the integrator itself, and a function to free the space when done.

Each integration function requires a random number generator to be supplied, and returns an estimate of the integral and its standard deviation.
The accuracy of the result is determined by the number of function calls specified by the user. If a known level of accuracy is required this
can be achieved by calling the integrator several times and averaging the individual results until the desired accuracy is obtained.

Random sample points used within the Monte Carlo routines are always chosen strictly within the integration region, so that endpoint singularities
are automatically avoided.

## VEGAS

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

use ffi::FFI;
use std::marker::PhantomData;
use std::mem::transmute;
use std::os::raw::c_void;
use std::slice;

ffi_wrapper!(PlainMonteCarlo, *mut sys::gsl_monte_plain_state, gsl_monte_plain_free,
"The plain Monte Carlo algorithm samples points randomly from the integration region to estimate
the integral and its error. Using this algorithm the estimate of the integral E(f; N) for N
randomly distributed points x_i is given by,

```text
E(f; N) = =  V <f> = (V / N) sum_i^N f(x_i)
```

where V is the volume of the integration region. The error on this estimate `sigma(E;N)` is
calculated from the estimated variance of the mean,

```text
sigma^2 (E; N) = (V^2 / N^2) sum_i^N (f(x_i) -  <f>)^2.
```

For large N this variance decreases asymptotically as `Var(f)/N`, where `Var(f)` is the true
variance of the function over the integration region. The error estimate itself should decrease as
`sigma(f)/sqrt{N}`. The familiar law of errors decreasing as `1/sqrt{N}` applies-to reduce the
error by a factor of 10 requires a 100-fold increase in the number of sample points.");

impl PlainMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions.
    #[doc(alias = "gsl_monte_plain_alloc")]
    pub fn new(dim: usize) -> Option<PlainMonteCarlo> {
        let tmp = unsafe { sys::gsl_monte_plain_alloc(dim) };

        if tmp.is_null() {
            None
        } else {
            Some(PlainMonteCarlo::wrap(tmp))
        }
    }

    /// This function initializes a previously allocated integration state. This allows an existing workspace to be reused for different
    /// integrations.
    #[doc(alias = "gsl_monte_plain_init")]
    pub fn init(&mut self) -> ::Value {
        ::Value::from(unsafe { sys::gsl_monte_plain_init(self.unwrap_unique()) })
    }

    /// This routines uses the plain Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined
    /// by the lower and upper limits in the arrays xl and xu, each of the same size. The integration uses a fixed number of function calls
    /// calls, and obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied.
    /// The result of the integration is returned in result, with an estimated absolute error abserr.
    ///
    /// In C, the function takes a `gsl_monte_function` as first argument. In here, you have to
    /// pass the `dim` argument and the function pointer (which became a closure) directly to the
    /// function.
    ///
    /// It returns either Ok((result, abserr)) or Err(Value).
    #[doc(alias = "gsl_monte_plain_integrate")]
    pub fn integrate<F: FnMut(&[f64]) -> f64>(
        &mut self,
        f: F,
        xl: &[f64],
        xu: &[f64],
        t_calls: usize,
        r: &mut ::Rng,
    ) -> (::Value, f64, f64) {
        assert!(xl.len() == xu.len());
        let mut result = 0f64;
        let mut abserr = 0f64;
        let f: Box<F> = Box::new(f);
        let ret = unsafe {
            let func = sys::gsl_monte_function {
                f: transmute(monte_trampoline::<F> as usize),
                dim: xl.len() as _,
                params: Box::into_raw(f) as *mut _,
            };
            sys::gsl_monte_plain_integrate(
                &func,
                xl.as_ptr(),
                xu.as_ptr(),
                xl.len() as _,
                t_calls,
                r.unwrap_unique(),
                self.unwrap_unique(),
                &mut result,
                &mut abserr,
            )
        };

        (::Value::from(ret), result, abserr)
    }
}

ffi_wrapper!(MiserMonteCarlo, *mut sys::gsl_monte_miser_state, gsl_monte_miser_free,
"The MISER algorithm of Press and Farrar is based on recursive stratified sampling. This technique
aims to reduce the overall integration error by concentrating integration points in the regions of
highest variance.

The idea of stratified sampling begins with the observation that for two disjoint regions a and b
with Monte Carlo estimates of the integral E_a(f) and E_b(f) and variances `sigma_a^2(f)` and
`sigma_b^2(f)`, the variance `Var(f)` of the combined estimate `E(f) = (1/2) (E_a(f) + E_b(f))`
is given by,

```text
Var(f) = (sigma_a^2(f) / 4 N_a) + (sigma_b^2(f) / 4 N_b).
```

It can be shown that this variance is minimized by distributing the points such that,

```text
N_a / (N_a + N_b) = sigma_a / (sigma_a + sigma_b).
```

Hence the smallest error estimate is obtained by allocating sample points in proportion to the
standard deviation of the function in each sub-region.

The MISER algorithm proceeds by bisecting the integration region along one coordinate axis to give
two sub-regions at each step. The direction is chosen by examining all d possible bisections and
selecting the one which will minimize the combined variance of the two sub-regions. The variance in
the sub-regions is estimated by sampling with a fraction of the total number of points available to
the current step. The same procedure is then repeated recursively for each of the two half-spaces
from the best bisection. The remaining sample points are allocated to the sub-regions using the
formula for N_a and N_b. This recursive allocation of integration points continues down to a
user-specified depth where each sub-region is integrated using a plain Monte Carlo estimate. These
individual values and their error estimates are then combined upwards to give an overall result and
an estimate of its error.");

impl MiserMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions. The workspace is used to maintain
    /// the state of the integration.
    #[doc(alias = "gsl_monte_miser_alloc")]
    pub fn new(dim: usize) -> Option<MiserMonteCarlo> {
        let tmp = unsafe { sys::gsl_monte_miser_alloc(dim) };

        if tmp.is_null() {
            None
        } else {
            Some(MiserMonteCarlo::wrap(tmp))
        }
    }

    /// This function initializes a previously allocated integration state. This allows an existing workspace to be reused for different integrations.
    #[doc(alias = "gsl_monte_miser_init")]
    pub fn init(&mut self) -> ::Value {
        ::Value::from(unsafe { sys::gsl_monte_miser_init(self.unwrap_unique()) })
    }

    /// This routines uses the MISER Monte Carlo algorithm to integrate the function f over the dim-dimensional hypercubic region defined by
    /// the lower and upper limits in the arrays xl and xu, each of size dim. The integration uses a fixed number of function calls calls,
    /// and obtains random sampling points using the random number generator r. A previously allocated workspace s must be supplied. The result
    /// of the integration is returned in result, with an estimated absolute error abserr.
    ///
    /// In C, the function takes a `gsl_monte_function` as first argument. In here, you have to
    /// pass the `dim` argument and the function pointer (which became a closure) directly to the
    /// function.
    ///
    /// It returns either Ok((result, abserr)) or Err(Value).
    #[doc(alias = "gsl_monte_miser_integrate")]
    pub fn integrate<F: FnMut(&[f64]) -> f64>(
        &mut self,
        f: F,
        xl: &[f64],
        xu: &[f64],
        t_calls: usize,
        r: &mut ::Rng,
    ) -> (::Value, f64, f64) {
        assert!(xl.len() == xu.len());
        let mut result = 0f64;
        let mut abserr = 0f64;
        let f: Box<F> = Box::new(f);
        let ret = unsafe {
            let mut func = sys::gsl_monte_function {
                f: transmute(monte_trampoline::<F> as usize),
                dim: xl.len() as _,
                params: Box::into_raw(f) as *mut _,
            };
            sys::gsl_monte_miser_integrate(
                &mut func,
                xl.as_ptr(),
                xu.as_ptr(),
                xl.len() as _,
                t_calls,
                r.unwrap_unique(),
                self.unwrap_unique(),
                &mut result,
                &mut abserr,
            )
        };
        (::Value::from(ret), result, abserr)
    }

    /// This function copies the parameters of the integrator state into the user-supplied params structure.
    #[doc(alias = "gsl_monte_miser_params_get")]
    pub fn get_params(&self) -> MiserParams {
        let mut m = sys::gsl_monte_miser_params {
            estimate_frac: 0f64,
            min_calls: 0,
            min_calls_per_bisection: 0,
            alpha: 0f64,
            dither: 0f64,
        };

        unsafe {
            sys::gsl_monte_miser_params_get(self.unwrap_shared(), &mut m);
        }
        MiserParams(m)
    }

    /// This function sets the integrator parameters based on values provided in the params structure.
    #[doc(alias = "gsl_monte_miser_params_set")]
    pub fn set_params(&mut self, params: &MiserParams) {
        unsafe {
            sys::gsl_monte_miser_params_set(self.unwrap_unique(), &params.0 as *const _);
        }
    }
}

#[derive(Debug, Clone)]
#[repr(C)]
pub struct MiserParams(pub sys::gsl_monte_miser_params);

ffi_wrapper!(VegasMonteCarlo, *mut sys::gsl_monte_vegas_state, gsl_monte_vegas_free,
"The VEGAS algorithm of Lepage is based on importance sampling. It samples points from the probability
distribution described by the function |f|, so that the points are concentrated in the regions that
make the largest contribution to the integral.

In general, if the Monte Carlo integral of f is sampled with points distributed according to a
probability distribution described by the function g, we obtain an estimate E_g(f; N),

```text
E_g(f; N) = E(f/g; N)
```

with a corresponding variance,

```text
Var_g(f; N) = Var(f/g; N).
```

If the probability distribution is chosen as g = |f|/I(|f|) then it can be shown that the variance
V_g(f; N) vanishes, and the error in the estimate will be zero. In practice it is not possible to
sample from the exact distribution g for an arbitrary function, so importance sampling algorithms
aim to produce efficient approximations to the desired distribution.

The VEGAS algorithm approximates the exact distribution by making a number of passes over the
integration region while histogramming the function f. Each histogram is used to define a sampling
distribution for the next pass. Asymptotically this procedure converges to the desired distribution.
In order to avoid the number of histogram bins growing like K^d the probability distribution is
approximated by a separable function: g(x_1, x_2, ...) = g_1(x_1) g_2(x_2) ... so that the number
of bins required is only Kd. This is equivalent to locating the peaks of the function from the
projections of the integrand onto the coordinate axes. The efficiency of VEGAS depends on the
validity of this assumption. It is most efficient when the peaks of the integrand are
well-localized. If an integrand can be rewritten in a form which is approximately separable this
will increase the efficiency of integration with VEGAS.

VEGAS incorporates a number of additional features, and combines both stratified sampling and
importance sampling. The integration region is divided into a number of “boxes”, with each box
getting a fixed number of points (the goal is 2). Each box can then have a fractional number of
bins, but if the ratio of bins-per-box is less than two, Vegas switches to a kind variance reduction
(rather than importance sampling).

The VEGAS algorithm computes a number of independent estimates of the integral internally, according
to the iterations parameter described below, and returns their weighted average. Random sampling of
the integrand can occasionally produce an estimate where the error is zero, particularly if the
function is constant in some regions. An estimate with zero error causes the weighted average to
break down and must be handled separately. In the original Fortran implementations of VEGAS the
error estimate is made non-zero by substituting a small value (typically 1e-30). The implementation
in GSL differs from this and avoids the use of an arbitrary constant—it either assigns the value a
weight which is the average weight of the preceding estimates or discards it according to the
following procedure,

current estimate has zero error, weighted average has finite error

* The current estimate is assigned a weight which is the average weight of the preceding estimates.
current estimate has finite error, previous estimates had zero error

* The previous estimates are discarded and the weighted averaging procedure begins with the current estimate.
current estimate has zero error, previous estimates had zero error

* The estimates are averaged using the arithmetic mean, but no error is computed.");

impl VegasMonteCarlo {
    /// This function allocates and initializes a workspace for Monte Carlo integration in dim dimensions.
    /// The workspace is used to maintain the state of the integration.
    #[doc(alias = "gsl_monte_vegas_alloc")]
    pub fn new(dim: usize) -> Option<VegasMonteCarlo> {
        let tmp = unsafe { sys::gsl_monte_vegas_alloc(dim) };

        if tmp.is_null() {
            None
        } else {
            Some(VegasMonteCarlo::wrap(tmp))
        }
    }

    /// This function initializes a previously allocated integration state. This allows an existing workspace
    /// to be reused for different integrations.
    #[doc(alias = "gsl_monte_vegas_init")]
    pub fn init(&mut self) -> ::Value {
        ::Value::from(unsafe { sys::gsl_monte_vegas_init(self.unwrap_unique()) })
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
    /// In C, the function takes a `gsl_monte_function` as first argument. In here, you have to
    /// pass the `dim` argument and the function pointer (which became a closure) directly to the
    /// function.
    ///
    /// It returns either Ok((result, abserr)) or Err(Value).
    #[doc(alias = "gsl_monte_vegas_integrate")]
    pub fn integrate<F: FnMut(&[f64]) -> f64>(
        &mut self,
        f: F,
        xl: &[f64],
        xu: &[f64],
        t_calls: usize,
        r: &mut ::Rng,
    ) -> (::Value, f64, f64) {
        assert!(xl.len() == xu.len());
        let mut result = 0f64;
        let mut abserr = 0f64;
        let f: Box<F> = Box::new(f);
        let ret = unsafe {
            let mut func = sys::gsl_monte_function {
                f: transmute(monte_trampoline::<F> as usize),
                dim: xl.len() as _,
                params: Box::into_raw(f) as *mut _,
            };
            sys::gsl_monte_vegas_integrate(
                &mut func,
                xl.as_ptr() as usize as *mut _,
                xu.as_ptr() as usize as *mut _,
                xl.len() as _,
                t_calls,
                r.unwrap_unique(),
                self.unwrap_unique(),
                &mut result,
                &mut abserr,
            )
        };
        (::Value::from(ret), result, abserr)
    }

    /// This function returns the chi-squared per degree of freedom for the weighted estimate of the integral.
    /// The returned value should be close to 1. A value which differs significantly from 1 indicates that
    /// the values from different iterations are inconsistent. In this case the weighted error will be
    /// under-estimated, and further iterations of the algorithm are needed to obtain reliable results.
    #[doc(alias = "gsl_monte_vegas_chisq")]
    pub fn chisq(&mut self) -> f64 {
        unsafe { sys::gsl_monte_vegas_chisq(self.unwrap_unique()) }
    }

    /// This function returns the raw (unaveraged) values of the integral result and its error sigma
    /// from the most recent iteration of the algorithm.
    ///
    /// Returns `(result, sigma)`.
    #[doc(alias = "gsl_monte_vegas_runval")]
    pub fn runval(&mut self) -> (f64, f64) {
        let mut result = 0.;
        let mut sigma = 0.;
        unsafe { sys::gsl_monte_vegas_runval(self.unwrap_unique(), &mut result, &mut sigma) };
        (result, sigma)
    }

    #[doc(alias = "gsl_monte_vegas_params_get")]
    pub fn get_params(&self) -> VegasParams {
        let mut params = VegasParams::default();
        unsafe {
            sys::gsl_monte_vegas_params_get(self.unwrap_shared(), &mut params.inner as *mut _);
        }
        params
    }

    #[doc(alias = "gsl_monte_vegas_params_set")]
    pub fn set_params(&mut self, params: &VegasParams) {
        unsafe {
            sys::gsl_monte_vegas_params_set(self.unwrap_unique(), &params.inner as *const _);
        }
    }
}

pub struct VegasParams<'a> {
    inner: sys::gsl_monte_vegas_params,
    lt: PhantomData<&'a ()>,
}

impl<'a> VegasParams<'a> {
    /// alpha: The parameter alpha controls the stiffness of the rebinning algorithm. It is typically
    /// set between one and two. A value of zero prevents rebinning of the grid. The default
    /// value is 1.5.
    ///
    /// iterations: The number of iterations to perform for each call to the routine. The default value
    /// is 5 iterations.
    ///
    /// stage: Setting this determines the stage of the calculation. Normally, stage = 0 which begins
    /// with a new uniform grid and empty weighted average. Calling vegas with stage =
    /// 1 retains the grid from the previous run but discards the weighted average, so that
    /// one can “tune” the grid using a relatively small number of points and then do a large
    /// run with stage = 1 on the optimized grid. Setting stage = 2 keeps the grid and the
    /// weighted average from the previous run, but may increase (or decrease) the number
    /// of histogram bins in the grid depending on the number of calls available. Choosing
    /// stage = 3 enters at the main loop, so that nothing is changed, and is equivalent to
    /// performing additional iterations in a previous call.
    ///
    /// mode: The possible choices are GSL_VEGAS_MODE_IMPORTANCE, GSL_VEGAS_MODE_
    /// STRATIFIED, GSL_VEGAS_MODE_IMPORTANCE_ONLY. This determines whether vegas
    /// will use importance sampling or stratified sampling, or whether it can pick on
    /// its own. In low dimensions vegas uses strict stratified sampling (more precisely,
    /// stratified sampling is chosen if there are fewer than 2 bins per box).
    ///
    /// verbosity + stream: These parameters set the level of information printed by vegas.
    pub fn new(
        alpha: f64,
        iterations: usize,
        stage: i32,
        mode: ::VegasMode,
        verbosity: VegasVerbosity,
        stream: Option<&'a mut ::IOStream>,
    ) -> Result<VegasParams, String> {
        if !verbosity.is_off() && stream.is_none() {
            return Err(
                "rust-GSL: need to provide an input stream for Vegas Monte Carlo \
                        integration if verbosity is not 'Off'"
                    .to_string(),
            );
        } else if verbosity.is_off() && stream.is_some() {
            return Err(
                "rust-GSL: need to provide the verbosity flag for Vegas Monta Carlo \
                        integration, currently set to 'Off'"
                    .to_string(),
            );
        }

        let stream = if let Some(stream) = stream {
            if !stream.write_mode() {
                return Err("rust-GSL: input stream not flagged as 'write' mode".to_string());
            }
            stream.as_raw()
        } else {
            ::std::ptr::null_mut()
        };
        Ok(VegasParams {
            inner: sys::gsl_monte_vegas_params {
                alpha,
                iterations,
                stage,
                mode: mode.into(),
                verbose: verbosity.to_int(),
                ostream: stream,
            },
            lt: PhantomData,
        })
    }
}

impl<'a> ::std::default::Default for VegasParams<'a> {
    fn default() -> VegasParams<'a> {
        VegasParams {
            inner: sys::gsl_monte_vegas_params {
                alpha: 1.5,
                iterations: 5,
                stage: 0,
                mode: ::VegasMode::ImportanceOnly.into(),
                verbose: -1,
                ostream: ::std::ptr::null_mut(),
            },
            lt: PhantomData,
        }
    }
}

/// The default setting of verbose is `Off`, which turns off all output.
/// A verbose value of `Summary` prints summary information about the weighted average
/// and final result, while a value of `Grid` also displays the grid coordinates.
/// A value of 'Rebinning' prints information from the rebinning procedure for each iteration.
#[derive(Clone, Copy)]
pub enum VegasVerbosity {
    Off,       // -1
    Summary,   // 0
    Grid,      // 1
    Rebinning, // 2
}

impl VegasVerbosity {
    fn to_int(self) -> i32 {
        match self {
            VegasVerbosity::Off => -1,
            VegasVerbosity::Summary => 0,
            VegasVerbosity::Grid => 1,
            VegasVerbosity::Rebinning => 2,
        }
    }

    fn is_off(self) -> bool {
        matches!(self, VegasVerbosity::Off)
    }
}

unsafe extern "C" fn monte_trampoline<F: FnMut(&[f64]) -> f64>(
    x: *mut f64,
    dim: usize,
    param: *mut c_void,
) -> f64 {
    let f: &mut F = &mut *(param as *mut F);
    f(slice::from_raw_parts(x, dim))
}

// The following tests have been made and tested against the following C code:
//
// ```ignore
// double
// g (double *k, size_t dim, void *params)
// {
//   (void)(dim);
//   (void)(params);
//   double A = 1.0 / (M_PI * M_PI * M_PI);
//   return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
// }
//
// void
// display_results (char *title, double result, double error)
// {
//   printf ("%s ==================\n", title);
//   printf ("result = % .6f\n", result);
//   printf ("sigma  = % .6f\n", error);
//   printf ("exact  = % .6f\n", exact);
//   printf ("error  = % .6f = %.2g sigma\n", result - exact,
//       fabs (result - exact) / error);
// }
//
// int main(void) {
//   double res, err;
//
//   double xl[3] = { 0, 0, 0 };
//   double xu[3] = { M_PI, M_PI, M_PI };
//
//   const gsl_rng_type *T;
//   gsl_rng *r;
//
//   gsl_monte_function G = { &g, 3, 0 };
//
//   size_t calls = 500000;
//
//   gsl_rng_env_setup ();
//
//   T = gsl_rng_default;
//   r = gsl_rng_alloc (T);
//
//   {
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);
//
//     gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
//                    &res, &err);
//     display_results ("vegas warm-up", res, err);
//
//     printf ("converging...\n");
//
//         do
//       {
//         gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
//                        &res, &err);
//         printf ("result = % .6f sigma = % .6f "
//             "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
//       }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
//
//     display_results ("vegas final", res, err);
//
//     gsl_monte_vegas_free (s);
//   }
//
//   gsl_rng_free (r);
//
//   return 0;
// }
// ```
#[test]
fn plain() {
    use std::f64::consts::PI;

    fn g(k: &[f64]) -> f64 {
        let a = 1f64 / (PI * PI * PI);

        a / (1.0 - k[0].cos() * k[1].cos() * k[2].cos())
    }

    let xl: [f64; 3] = [0f64; 3];
    let xu: [f64; 3] = [PI, PI, PI];

    let calls = 500000;

    ::RngType::env_setup();
    let mut r = ::Rng::new(::RngType::default()).unwrap();

    {
        let mut s = PlainMonteCarlo::new(3).unwrap();

        let (_, res, err) = s.integrate(g, &xl, &xu, calls, &mut r);
        assert_eq!(&format!("{:.6}", res), "1.412209");
        assert_eq!(&format!("{:.6}", err), "0.013436");
    }
}

#[test]
fn miser() {
    use std::f64::consts::PI;
    fn g(k: &[f64]) -> f64 {
        let a = 1f64 / (PI * PI * PI);

        a / (1.0 - k[0].cos() * k[1].cos() * k[2].cos())
    }

    let xl: [f64; 3] = [0f64; 3];
    let xu: [f64; 3] = [PI, PI, PI];

    let calls = 500000;

    ::RngType::env_setup();
    let mut r = ::Rng::new(::RngType::default()).unwrap();

    {
        let mut s = MiserMonteCarlo::new(3).unwrap();

        let (_, res, err) = s.integrate(g, &xl, &xu, calls, &mut r);
        assert_eq!(&format!("{:.6}", res), "1.389530");
        assert_eq!(&format!("{:.6}", err), "0.005011");
    }
}

#[test]
fn miser_closure() {
    use std::f64::consts::PI;

    let xl: [f64; 3] = [0f64; 3];
    let xu: [f64; 3] = [PI, PI, PI];

    let calls = 500000;

    ::RngType::env_setup();
    let mut r = ::Rng::new(::RngType::default()).unwrap();

    {
        let mut s = MiserMonteCarlo::new(3).unwrap();

        let (_, res, err) = s.integrate(
            |k| {
                let a = 1f64 / (PI * PI * PI);

                a / (1.0 - k[0].cos() * k[1].cos() * k[2].cos())
            },
            &xl,
            &xu,
            calls,
            &mut r,
        );
        assert_eq!(&format!("{:.6}", res), "1.389530");
        assert_eq!(&format!("{:.6}", err), "0.005011");
    }
}

#[test]
fn vegas_warm_up() {
    use std::f64::consts::PI;
    fn g(k: &[f64]) -> f64 {
        let a = 1f64 / (PI * PI * PI);

        a / (1.0 - k[0].cos() * k[1].cos() * k[2].cos())
    }

    let xl: [f64; 3] = [0f64; 3];
    let xu: [f64; 3] = [PI, PI, PI];

    ::RngType::env_setup();
    let mut r = ::Rng::new(::RngType::default()).unwrap();

    {
        let mut s = VegasMonteCarlo::new(3).unwrap();

        let (_, res, err) = s.integrate(g, &xl, &xu, 10000, &mut r);
        assert_eq!(&format!("{:.6}", res), "1.385603");
        assert_eq!(&format!("{:.6}", err), "0.002212");
    }
}

#[test]
fn vegas() {
    use std::f64::consts::PI;
    fn g(k: &[f64]) -> f64 {
        let a = 1. / (PI * PI * PI);

        a / (1. - k[0].cos() * k[1].cos() * k[2].cos())
    }

    let calls = 500000;

    let xl: [f64; 3] = [0f64; 3];
    let xu: [f64; 3] = [PI, PI, PI];

    ::RngType::env_setup();
    let mut r = ::Rng::new(::RngType::default()).unwrap();

    {
        let mut s = VegasMonteCarlo::new(3).unwrap();

        s.integrate(g, &xl, &xu, 10000, &mut r);
        let mut res;
        let mut err;
        loop {
            let (_, _res, _err) = s.integrate(g, &xl, &xu, calls / 5, &mut r);
            res = _res;
            err = _err;
            println!(
                "result = {:.6} sigma = {:.6} chisq/dof = {:.1}",
                res,
                err,
                s.chisq()
            );
            if (s.chisq() - 1f64).abs() <= 0.5f64 {
                break;
            }
        }
        assert_eq!(&format!("{:.6}", res), "1.393307");
        assert_eq!(&format!("{:.6}", err), "0.000335");
    }
}
