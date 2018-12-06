//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The following functions allow the shuffling and sampling of a set of objects.
The algorithms rely on a random number generator as a source of randomness and a poor quality generator can lead to correlations in the output.
In particular it is important to avoid generators with a short period. For more information see Knuth, v2, 3rd ed, Section 3.4.2, “Random Sampling and Shuffling”.
!*/

use ffi;
use types::Rng;
use enums;
use libc::c_void;

/// This function randomly shuffles the order of n objects, each of size size, stored in the array base[0..n-1]. The output of the random number generator r is used to
/// produce the permutation. The algorithm generates all possible n! permutations with equal probability, assuming a perfect source of random numbers.
/// 
/// The following code shows how to shuffle the numbers from 0 to 51,
/// 
/// ```C
/// int a[52];
/// 
/// for (i = 0; i < 52; i++)
///   {
///     a[i] = i;
///   }
/// 
/// gsl_ran_shuffle (r, a, 52, sizeof (int));
/// ```
pub fn shuffle<T>(r: &mut Rng, base: &mut [T]) {
    unsafe { ffi::gsl_ran_shuffle(ffi::FFI::unwrap_unique(r),
        base.as_mut_ptr() as *mut c_void,
        base.len() as usize,
        ::std::mem::size_of::<T>() as usize) }
}

/// This function fills the array dest[k] with k objects taken randomly from the n elements of the array src[0..n-1]. The objects are each of size size.
/// The output of the random number generator r is used to make the selection. The algorithm ensures all possible samples are equally likely, assuming a perfect source of randomness.
/// 
/// The objects are sampled without replacement, thus each object can only appear once in dest[k]. It is required that k be less than or equal to n.
/// The objects in dest will be in the same relative order as those in src. You will need to call gsl_ran_shuffle(r, dest, n, size) if you want to randomize the order.
/// 
/// The following code shows how to select a random sample of three unique numbers from the set 0 to 99,
/// 
/// ```C
/// double a[3], b[100];
/// 
/// for (i = 0; i < 100; i++)
///   {
///     b[i] = (double) i;
///   }
/// 
/// gsl_ran_choose (r, a, 3, b, 100, sizeof (double));
/// ```
pub fn choose<T>(r: &mut Rng, dest: &mut [T], src: &[T]) -> enums::Value {
    enums::Value::from(unsafe { ffi::gsl_ran_choose(ffi::FFI::unwrap_unique(r),
                                                    dest.as_mut_ptr() as *mut c_void,
                                                    dest.len() as usize,
                                                    src.as_ptr() as *mut c_void,
                                                    src.len() as usize,
                                                    ::std::mem::size_of::<T>() as usize)
    })
}

/// This function is like gsl_ran_choose but samples k items from the original array of n items src with replacement, so the same object can appear more
/// than once in the output sequence dest. There is no requirement that k be less than n in this case.
pub fn sample<T>(r: &mut Rng, dest: &mut [T], src: &[T]) -> enums::Value {
    enums::Value::from(unsafe { ffi::gsl_ran_sample(ffi::FFI::unwrap_unique(r),
                                                    dest.as_mut_ptr() as *mut c_void,
                                                    dest.len() as usize,
                                                    src.as_ptr() as *mut c_void,
                                                    src.len() as usize,
                                                    ::std::mem::size_of::<T>() as usize)
    })
}
