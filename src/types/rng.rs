//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Random Number Generation

The library provides a large collection of random number generators which can be accessed through a uniform interface. 
Environment variables allow you to select different generators and seeds at runtime, so that you can easily switch between generators without needing to recompile your program.
Each instance of a generator keeps track of its own state, allowing the generators to be used in multi-threaded programs.
Additional functions are available for transforming uniform random numbers into samples from continuous or discrete probability distributions such as the Gaussian, log-normal or Poisson distributions.

##General comments on random numbers

In 1988, Park and Miller wrote a paper entitled “Random number generators: good ones are hard to find.” [Commun. ACM, 31, 1192–1201]. Fortunately, some excellent random number generators are available, though poor ones are still in common use. You may be happy with the system-supplied random number generator on your computer, but you should be aware that as computers get faster, requirements on random number generators increase. Nowadays, a simulation that calls a random number generator millions of times can often finish before you can make it down the hall to the coffee machine and back.

A very nice review of random number generators was written by Pierre L’Ecuyer, as Chapter 4 of the book: Handbook on Simulation, Jerry Banks, ed. (Wiley, 1997). The chapter is available in postscript from L’Ecuyer’s ftp site (see references). Knuth’s volume on Seminumerical Algorithms (originally published in 1968) devotes 170 pages to random number generators, and has recently been updated in its 3rd edition (1997). It is brilliant, a classic. If you don’t own it, you should stop reading right now, run to the nearest bookstore, and buy it.

A good random number generator will satisfy both theoretical and statistical properties. Theoretical properties are often hard to obtain (they require real math!), but one prefers a random number generator with a long period, low serial correlation, and a tendency not to “fall mainly on the planes.” Statistical tests are performed with numerical simulations. Generally, a random number generator is used to estimate some quantity for which the theory of probability provides an exact answer. Comparison to this exact answer provides a measure of “randomness”.

##The Random Number Generator Interface

It is important to remember that a random number generator is not a “real” function like sine or cosine. Unlike real functions, successive calls to a random number generator yield different return values. Of course that is just what you want for a random number generator, but to achieve this effect, the generator must keep track of some kind of “state” variable.
Sometimes this state is just an integer (sometimes just the value of the previously generated random number), but often it is more complicated than that and may involve a whole array of numbers, possibly with some indices thrown in. To use the random number generators, you do not need to know the details of what comprises the state, and besides that varies from algorithm to algorithm.

The random number generator library uses two special structs, RngType which holds static information about each type of generator and Rng which describes an instance of a generator created from a given RngType.

##Performance

The following table shows the relative performance of a selection the available random number generators. The fastest simulation quality generators are taus, gfsr4 and mt19937. The generators which offer the best mathematically-proven quality are those based on the RANLUX algorithm.

 * 1754 k ints/sec,    870 k doubles/sec, taus
 * 1613 k ints/sec,    855 k doubles/sec, gfsr4
 * 1370 k ints/sec,    769 k doubles/sec, mt19937
 *  565 k ints/sec,    571 k doubles/sec, ranlxs0
 *  400 k ints/sec,    405 k doubles/sec, ranlxs1
 *  490 k ints/sec,    389 k doubles/sec, mrg
 *  407 k ints/sec,    297 k doubles/sec, ranlux
 *  243 k ints/sec,    254 k doubles/sec, ranlxd1
 *  251 k ints/sec,    253 k doubles/sec, ranlxs2
 *  238 k ints/sec,    215 k doubles/sec, cmrg
 *  247 k ints/sec,    198 k doubles/sec, ranlux389
 *  141 k ints/sec,    140 k doubles/sec, ranlxd2

##Random number environment variables

The library allows you to choose a default generator and seed from the environment variables GSL_RNG_TYPE and GSL_RNG_SEED and the function gsl_rng_env_setup. This makes it easy try out different generators and seeds without having to recompile your program.

##References and Further Reading

The subject of random number generation and testing is reviewed extensively in Knuth’s Seminumerical Algorithms.

Donald E. Knuth, The Art of Computer Programming: Seminumerical Algorithms (Vol 2, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896842.
Further information is available in the review paper written by Pierre L’Ecuyer,

P. L’Ecuyer, “Random Number Generation”, Chapter 4 of the Handbook on Simulation, Jerry Banks Ed., Wiley, 1998, 93–137.
http://www.iro.umontreal.ca/~lecuyer/papers.html

The source code for the DIEHARD random number generator tests is also available online,

DIEHARD source code G. Marsaglia,
http://stat.fsu.edu/pub/diehard/
A comprehensive set of random number generator tests is available from NIST,

NIST Special Publication 800-22, “A Statistical Test Suite for the Validation of Random Number Generators and Pseudo Random Number Generators for Cryptographic Applications”.
http://csrc.nist.gov/rng/

##Acknowledgements

Thanks to Makoto Matsumoto, Takuji Nishimura and Yoshiharu Kurita for making the source code to their generators (MT19937, MM&TN; TT800, MM&YK) available under the GNU General Public License. Thanks to Martin Lüscher for providing notes and source code for the RANLXS and RANLXD generators.
!*/

use ffi;
use std::string;
use std::default::Default;
use enums;

pub struct Rng {
    r: *mut ffi::gsl_rng
}

impl Rng {
    /// This function returns a pointer to a newly-created instance of a random number generator of type T. For example, the following code creates an instance of the Tausworthe generator,
    /// 
    /// ```Rust
    /// let r = Rng::new(gsl_rng_taus);
    /// ```
    /// 
    /// If there is insufficient memory to create the generator then the function returns a null pointer and the error handler is invoked with an error code of GSL_ENOMEM.
    /// 
    /// The generator is automatically initialized with the default seed, gsl_rng_default_seed. This is zero by default but can be changed either directly or by using the environment variable
    /// GSL_RNG_SEED (see [`Random number environment variables`](https://www.gnu.org/software/gsl/manual/html_node/Random-number-environment-variables.html#Random-number-environment-variables)).
    pub fn new(T: &RngType) -> Option<Rng> {
        let tmp = unsafe { ffi::gsl_rng_alloc(ffi::FFI::unwrap(T) as *const ffi::gsl_rng_type) };

        if tmp.is_null() {
            None
        } else {
            Some(Rng {
                r: tmp
            })
        }
    }

    /// This function initializes (or ‘seeds’) the random number generator. If the generator is seeded with the same value of s on two different runs, the same stream of random numbers will be generated by successive calls to the routines below.
    /// If different values of s >= 1 are supplied, then the generated streams of random numbers should be completely different. If the seed s is zero then the standard seed from the original implementation is used instead.
    /// For example, the original Fortran source code for the ranlux generator used a seed of 314159265, and so choosing s equal to zero reproduces this when using gsl_rng_ranlux.
    /// 
    /// When using multiple seeds with the same generator, choose seed values greater than zero to avoid collisions with the default setting.
    /// 
    /// Note that the most generators only accept 32-bit seeds, with higher values being reduced modulo 2^32. For generators with smaller ranges the maximum seed value will typically be lower.
    pub fn set(&self, s: u64) {
        unsafe { ffi::gsl_rng_set(self.r as *const ffi::gsl_rng, s) }
    }

    /// This function returns a random integer from the generator r. The minimum and maximum values depend on the algorithm used, but all integers in the range [min,max] are equally likely.
    /// The values of min and max can be determined using the auxiliary functions gsl_rng_max (r) and gsl_rng_min (r).
    pub fn get(&self) -> u64 {
        unsafe { ffi::gsl_rng_get(self.r as *const ffi::gsl_rng) }
    }

    /// This function returns a double precision floating point number uniformly distributed in the range [0,1). The range includes 0.0 but excludes 1.0.
    /// The value is typically obtained by dividing the result of gsl_rng_get(r) by gsl_rng_max(r) + 1.0 in double precision.
    /// Some generators compute this ratio internally so that they can provide floating point numbers with more than 32 bits of randomness (the maximum number of bits that can be portably represented in a single unsigned long int).
    pub fn uniform(&self) -> f64 {
        unsafe { ffi::gsl_rng_uniform(self.r as *const ffi::gsl_rng) }
    }

    /// This function returns a positive double precision floating point number uniformly distributed in the range (0,1), excluding both 0.0 and 1.0.
    /// The number is obtained by sampling the generator with the algorithm of gsl_rng_uniform until a non-zero value is obtained.
    /// You can use this function if you need to avoid a singularity at 0.0.
    pub fn uniform_pos(&self) -> f64 {
        unsafe { ffi::gsl_rng_uniform_pos(self.r as *const ffi::gsl_rng) }
    }

    /// This function returns a random integer from 0 to n-1 inclusive by scaling down and/or discarding samples from the generator r.
    /// All integers in the range [0,n-1] are produced with equal probability. For generators with a non-zero minimum value an offset is applied so that zero is returned with the correct probability.
    /// 
    /// Note that this function is designed for sampling from ranges smaller than the range of the underlying generator. The parameter n must be less than or equal to the range of the generator r.
    /// If n is larger than the range of the generator then the function calls the error handler with an error code of GSL_EINVAL and returns zero.
    /// 
    /// In particular, this function is not intended for generating the full range of unsigned integer values [0,2^32-1].
    /// Instead choose a generator with the maximal integer range and zero minimum value, such as gsl_rng_ranlxd1, gsl_rng_mt19937 or gsl_rng_taus, and sample it directly using gsl_rng_get. The range of each generator can be found using the auxiliary functions described in the next section.
    pub fn uniform_int(&self, n: u64) -> u64 {
        unsafe { ffi::gsl_rng_uniform_int(self.r as *const ffi::gsl_rng, n) }
    }

    /// This function returns a pointer to the name of the generator. For example,
    /// 
    /// ```Rust
    /// println!("r is a '{}' generator", r.get_name());
    /// ```
    /// 
    /// would print something like "r is a 'taus' generator".
    pub fn get_name(&self) -> String {
        unsafe { string::raw::from_buf(ffi::gsl_rng_name(self.r as *const ffi::gsl_rng) as *const u8) }
    }

    /// This function returns the largest value that the get function can return.
    pub fn max(&self) -> u64 {
        unsafe { ffi::gsl_rng_max(self.r as *const ffi::gsl_rng) }
    }

    /// This function returns the smallest value that gsl_rng_get can return. Usually this value is zero.
    /// There are some generators with algorithms that cannot return zero, and for these generators the minimum value is 1.
    pub fn min(&self) -> u64 {
        unsafe { ffi::gsl_rng_min(self.r as *const ffi::gsl_rng) }
    }

    /// This function returns a pointer to the state of generator r. You can use this information to access the state directly. For example, the following code will write the state of a generator to a stream,
    /// 
    /// ```C
    /// void * state = gsl_rng_state (r);
    /// size_t n = gsl_rng_size (r);
    /// fwrite (state, n, 1, stream);
    /// ```
    pub fn state<'r, T>(&self) -> &'r mut T {
        unsafe { ::std::mem::transmute(ffi::gsl_rng_state(self.r as *const ffi::gsl_rng)) }
    }

    /// This function copies the random number generator src into the pre-existing generator dest, making dest into an exact copy of src. The two generators must be of the same type.
    pub fn copy(&self, other: &Rng) -> enums::Value {
        unsafe { ffi::gsl_rng_memcpy(other.r, self.r as *const ffi::gsl_rng) }
    }

    /// This function returns the size of the state of generator r. You can use this information to access the state directly. For example, the following code will write the state of a generator to a stream,
    /// 
    /// ```C
    /// void * state = gsl_rng_state (r);
    /// size_t n = gsl_rng_size (r);
    /// fwrite (state, n, 1, stream);
    /// ```
    pub fn size(&self) -> u64 {
        unsafe { ffi::gsl_rng_size(self.r as *const ffi::gsl_rng) }
    }

    /// Equivalent to DefaultRngSeed
    pub fn default_seed() -> u64 {
        ffi::gsl_rng_default_seed
    }
}

impl Clone for Rng {
    /// This function returns a pointer to a newly created generator which is an exact copy of the generator r.
    fn clone(&self) -> Rng {
        unsafe { ffi::FFI::wrap(ffi::gsl_rng_clone(self.r as *const ffi::gsl_rng)) }
    }
}

impl Drop for Rng {
    fn drop(&mut self) {
        unsafe { ffi::gsl_rng_free(self.r as *const ffi::gsl_rng) };
        self.r = ::std::ptr::mut_null();
    }
}

impl ffi::FFI<ffi::gsl_rng> for Rng {
    fn wrap(r: *mut ffi::gsl_rng) -> Rng {
        Rng {
            r: r
        }
    }

    fn unwrap(m: &Rng) -> *mut ffi::gsl_rng {
        m.r
    }
}

pub struct RngType {
    ptr: *mut ffi::gsl_rng_type
}

impl RngType {
    /// wrapper for name element
    pub fn name(&self) -> String {
        if self.ptr.is_null() {
            String::new()
        } else {
            unsafe { ::std::string::raw::from_buf((*self.ptr).name as *const u8) }
        }
    }

    /// wrapper for max element
    pub fn max(&self) -> u64 {
        if self.ptr.is_null() {
            0u64
        } else {
            unsafe { (*self.ptr).max }
        }
    }

    /// wrapper for min element
    pub fn min(&self) -> u64 {
        if self.ptr.is_null() {
            0u64
        } else {
            unsafe { (*self.ptr).min }
        }
    }

    /// wrapper for size element
    pub fn size(&self) -> u64 {
        if self.ptr.is_null() {
            0u64
        } else {
            unsafe { (*self.ptr).size }
        }
    }

    /// This function returns a pointer to an array of all the available generator types, terminated by a null pointer.
    /// The function should be called once at the start of the program, if needed. The following code fragment shows how to iterate over the array of generator types to print the names of the available algorithms,
    /// 
    /// ```Rust
    /// let t = RngType::types_setup ();
    /// 
    /// println!("Available generators:");
    /// for tmp in t.iter() {
    ///     println!("{}", tmp.name);
    /// }
    /// ```
    pub fn types_setup() -> Vec<RngType> {
        let ptr = unsafe { ffi::gsl_rng_types_setup() };
        let mut ret = Vec::new();

        if ptr.is_not_null() {
            unsafe {
                let mut it = 0;
                loop {
                    let tmp = ptr.offset(it);

                    if (*tmp).is_null() {
                        break;
                    }
                    ret.push(ffi::FFI::wrap(*tmp));
                    it += 1;
                }
            }
        }
        ret
    }

    /// This function reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED and uses their values to set the corresponding library variables gsl_rng_default and gsl_rng_default_seed. These global variables are defined as follows,
    /// 
    /// ```C
    /// extern const gsl_rng_type *gsl_rng_default
    /// extern unsigned long int gsl_rng_default_seed
    /// ```
    /// 
    /// The environment variable GSL_RNG_TYPE should be the name of a generator, such as taus or mt19937. The environment variable GSL_RNG_SEED should contain the desired seed value.
    /// It is converted to an unsigned long int using the C library function strtoul.
    /// 
    /// If you don’t specify a generator for GSL_RNG_TYPE then gsl_rng_mt19937 is used as the default. The initial value of gsl_rng_default_seed is zero.
    /// See rng example in examples folder for more details.
    pub fn env_setup() -> Option<RngType> {
        let tmp = unsafe { ffi::gsl_rng_env_setup() };

        if tmp.is_null() {
            None
        } else {
            Some(ffi::FFI::wrap(tmp as *mut ffi::gsl_rng_type))
        }
    }
}

impl Default for RngType {
    fn default() -> RngType {
        ffi::FFI::wrap(ffi::gsl_rng_default as *mut ffi::gsl_rng_type)
    }
}

impl ffi::FFI<ffi::gsl_rng_type> for RngType {
    fn wrap(r: *mut ffi::gsl_rng_type) -> RngType {
        RngType {
            ptr: r
        }
    }

    fn unwrap(m: &RngType) -> *mut ffi::gsl_rng_type {
        m.ptr
    }
}