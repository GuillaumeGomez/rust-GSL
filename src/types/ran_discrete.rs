//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
Given K discrete events with different probabilities `P[k]`, produce a random value k consistent with its probability.

The obvious way to do this is to preprocess the probability list by generating a cumulative probability array with K+1 elements:

```text
  C[0] = 0
C[k+1] = C[k]+P[k].
```

Note that this construction produces `C[K]=1`. Now choose a uniform deviate u between 0 and 1, and find the value of k such that `C[k] <= u < C[k+1]`. Although this in principle requires of order \log K steps per random number generation, they are fast steps, and if you use something like \lfloor uK \rfloor as a starting point, you can often do pretty well.

But faster methods have been devised. Again, the idea is to preprocess the probability list, and save the result in some form of lookup table; then the individual calls for a random discrete event can go rapidly. An approach invented by G. Marsaglia (Generating discrete random variables in a computer, Comm ACM 6, 37–38 (1963)) is very clever, and readers interested in examples of good algorithm design are directed to this short and well-written paper. Unfortunately, for large K, Marsaglia’s lookup table can be quite large.

A much better approach is due to Alastair J. Walker (An efficient method for generating discrete random variables with general distributions, ACM Trans on Mathematical Software 3, 253–256 (1977); see also Knuth, v2, 3rd ed, p120–121,139). This requires two lookup tables, one floating point and one integer, but both only of size K. After preprocessing, the random numbers are generated in O(1) time, even for large K. The preprocessing suggested by Walker requires O(K^2) effort, but that is not actually necessary, and the implementation provided here only takes O(K) effort. In general, more preprocessing leads to faster generation of the individual random numbers, but a diminishing return is reached pretty early. Knuth points out that the optimal preprocessing is combinatorially difficult for large K.

This method can be used to speed up some of the discrete random number generators below, such as the binomial distribution. To use it for something like the Poisson Distribution, a modification would have to be made, since it only takes a finite set of K outcomes.
!*/

use ffi::FFI;
use types::Rng;

ffi_wrapper!(
    RanDiscrete,
    *mut sys::gsl_ran_discrete_t,
    gsl_ran_discrete_free
);

impl RanDiscrete {
    /// This function returns a pointer to a structure that contains the lookup table for the discrete random number generator. The array P[] contains the probabilities of the discrete events;
    /// these array elements must all be positive, but they needn’t add up to one (so you can think of them more generally as “weights”)—the preprocessor will normalize appropriately.
    /// This return value is used as an argument for the gsl_ran_discrete function below.
    pub fn new(P: &[f64]) -> Option<RanDiscrete> {
        let tmp = unsafe { sys::gsl_ran_discrete_preproc(P.len() as _, P.as_ptr()) };

        if tmp.is_null() {
            None
        } else {
            Some(RanDiscrete::wrap(tmp))
        }
    }

    /// After the new, above, has been called, you use this function to get the discrete random numbers.
    pub fn discrete(&self, r: &mut Rng) -> usize {
        unsafe { sys::gsl_ran_discrete(r.unwrap_unique(), self.unwrap_shared()) }
    }

    /// Returns the probability `P[k]` of observing the variable k. Since `P[k]` is not
    /// stored as part of the lookup table, it must be recomputed; this computation takes O(K),
    /// so if K is large and you care about the original array `P[k]` used to create the lookup
    /// table, then you should just keep this original array `P[k]` around.
    pub fn discrete_pdf(&self, k: usize) -> f64 {
        unsafe { sys::gsl_ran_discrete_pdf(k, self.unwrap_shared()) }
    }
}
