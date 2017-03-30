//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function returns a random integer from the hypergeometric distribution. The probability distribution for hypergeometric random variates is,
/// 
/// p(k) =  C(n_1, k) C(n_2, t - k) / C(n_1 + n_2, t)
/// 
/// where C(a,b) = a!/(b!(a-b)!) and t <= n_1 + n_2. The domain of k is max(0,t-n_2), ..., min(t,n_1).
/// 
/// If a population contains n_1 elements of “type 1” and n_2 elements of “type 2” then the hypergeometric distribution gives the probability of obtaining
/// k elements of “type 1” in t samples from the population without replacement.
pub fn hypergeometric(r: &mut Rng, n1: u32, n2: u32, t: u32) -> u32 {
    unsafe { ffi::gsl_ran_hypergeometric(ffi::FFI::unwrap(r), n1, n2, t) }
}

/// This function computes the probability p(k) of obtaining k from a hypergeometric distribution with parameters n1, n2, t, using the formula given above.
pub fn hypergeometric_pdf(k: u32, n1: u32, n2: u32, t: u32) -> f64 {
    unsafe { ffi::gsl_ran_hypergeometric_pdf(k, n1, n2, t) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the hypergeometric distribution with parameters n1, n2 and t.
pub fn hypergeometric_P(k: u32, n1: u32, n2: u32, t: u32) -> f64 {
    unsafe { ffi::gsl_cdf_hypergeometric_P(k, n1, n2, t) }
}

/// This function computes the cumulative distribution functions P(k), Q(k) for the hypergeometric distribution with parameters n1, n2 and t.
pub fn hypergeometric_Q(k: u32, n1: u32, n2: u32, t: u32) -> f64 {
    unsafe { ffi::gsl_cdf_hypergeometric_Q(k, n1, n2, t) }
}