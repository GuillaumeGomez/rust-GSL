//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use ffi;
use types::Rng;

/// This function computes a random sample n[] from the multinomial distribution formed by N trials from an underlying distribution p[K]. The distribution function for n[] is,
/// 
/// P(n_1, n_2, ..., n_K) = 
///   (N!/(n_1! n_2! ... n_K!)) p_1^n_1 p_2^n_2 ... p_K^n_K
/// 
/// where (n_1, n_2, ..., n_K) are nonnegative integers with sum_{k=1}^K n_k = N, and (p_1, p_2, ..., p_K) is a probability distribution with \sum p_i = 1. If the array p[K] is not normalized then its
/// entries will be treated as weights and normalized appropriately. The arrays n[] and p[] must both be of length K.
/// 
/// Random variates are generated using the conditional binomial method (see C.S. Davis, The computer generation of multinomial random variates, Comp. Stat. Data Anal. 16 (1993) 205–217 for details).
pub fn multinomial(r: &Rng, N: u32, p: &[f64], n: &mut [u32]) {
    unsafe { ffi::gsl_ran_multinomial(ffi::FFI::unwrap(r), p.len() as u64, N, p.as_ptr(), n.as_mut_ptr()) }
}

/// This function computes the probability P(n_1, n_2, ..., n_K) of sampling n[K] from a multinomial distribution with parameters p[K], using the formula given above.
pub fn multinomial_pdf(p: &[f64], n: &[u32]) -> f64 {
    unsafe { ffi::gsl_ran_multinomial_pdf(p.len() as u64, p.as_ptr(), n.as_ptr()) }
}

/// This function returns the logarithm of the probability for the multinomial distribution P(n_1, n_2, ..., n_K) with parameters p[K].
pub fn multinomial_lnpdf(p: &[f64], n: &[u32]) -> f64 {
    unsafe { ffi::gsl_ran_multinomial_lnpdf(p.len() as u64, p.as_ptr(), n.as_ptr()) }
}