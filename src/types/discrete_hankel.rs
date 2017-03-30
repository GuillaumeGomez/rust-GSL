//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
#Discrete Hankel Transforms

This chapter describes functions for performing Discrete Hankel Transforms (DHTs).

##Definitions

The discrete Hankel transform acts on a vector of sampled data, where the samples are assumed to
have been taken at points related to the zeroes of a Bessel function of fixed order; compare this to
the case of the discrete Fourier transform, where samples are taken at points related to the zeroes
of the sine or cosine function.

Specifically, let f(t) be a function on the unit interval and j_(\nu,m) the m-th zero of the Bessel
function J_\nu(x). Then the finite \nu-Hankel transform of f(t) is defined to be the set of numbers
g_m given by,

g_m = \int_0^1 t dt J_\nu(j_(\nu,m)t) f(t),

so that,

f(t) = \sum_{m=1}^\infty (2 J_\nu(j_(\nu,m)t) / J_(\nu+1)(j_(\nu,m))^2) g_m.

Suppose that f is band-limited in the sense that g_m=0 for m > M. Then we have the following
fundamental sampling theorem.

g_m = (2 / j_(\nu,M)^2)
      \sum_{k=1}^{M-1} f(j_(\nu,k)/j_(\nu,M))
          (J_\nu(j_(\nu,m) j_(\nu,k) / j_(\nu,M)) / J_(\nu+1)(j_(\nu,k))^2).

It is this discrete expression which defines the discrete Hankel transform. The kernel in the
summation above defines the matrix of the \nu-Hankel transform of size M-1. The coefficients of this
matrix, being dependent on \nu and M, must be precomputed and stored; the gsl_dht object
encapsulates this data. The allocation function gsl_dht_alloc returns a gsl_dht object which must be
properly initialized with gsl_dht_init before it can be used to perform transforms on data sample
vectors, for fixed \nu and M, using the gsl_dht_apply function. The implementation allows a scaling
of the fundamental interval, for convenience, so that one can assume the function is defined on the
interval [0,X], rather than the unit interval.

Notice that by assumption f(t) vanishes at the endpoints of the interval, consistent with the
inversion formula and the sampling formula given above. Therefore, this transform corresponds to an
orthogonal expansion in eigenfunctions of the Dirichlet problem for the Bessel differential
equation.

##References and Further Reading

The algorithms used by these functions are described in the following papers,

H. Fisk Johnson, Comp. Phys. Comm. 43, 181 (1987).
D. Lemoine, J. Chem. Phys. 101, 3936 (1994).
!*/

use ffi;
use enums;

pub struct DiscreteHankel {
    t: *mut ffi::gsl_dht,
}

impl DiscreteHankel {
    /// This function allocates a Discrete Hankel transform object of size `size`.
    pub fn new(size: usize) -> Option<DiscreteHankel> {
        let tmp = unsafe { ffi::gsl_dht_alloc(size) };

        if tmp.is_null() {
            None
        } else {
            Some(DiscreteHankel {
                t: tmp,
            })
        }
    }

    /// This function allocates a Discrete Hankel transform object of size `size` and initializes it
    /// for the given values of `nu` and `xmax`.
    pub fn new_with_init(size: usize, nu: f64, xmax: f64) -> Option<DiscreteHankel> {
        let tmp = unsafe { ffi::gsl_dht_new(size, nu, xmax) };

        if tmp.is_null() {
            None
        } else {
            Some(DiscreteHankel {
                t: tmp,
            })
        }
    }

    /// This function initializes the transform `self` for the given values of `nu` and `xmax`.
    pub fn init(&mut self, nu: f64, xmax: f64) -> enums::Value {
        unsafe { ffi::gsl_dht_init(self.t, nu, xmax) }
    }

    /// This function applies the transform t to the array f_in whose size is equal to the size of
    /// the transform. The result is stored in the array `f_out` which must be of the same length.
    ///
    /// Applying this function to its output gives the original data multiplied by (1/j_(\nu,M))^2,
    /// up to numerical errors.
    pub fn apply(&self, f_in: &[f64]) -> Result<Vec<f64>, enums::Value> {
        unsafe {
            assert!((*self.t).size == f_in.len(),
                    "f_in and f_out must have the same length as this struct");
            let mut f_out: Vec<f64> = ::std::iter::repeat(0.).take(f_in.len()).collect();
            match ffi::gsl_dht_apply(self.t, f_in.as_ptr(), f_out.as_mut_ptr()) {
                enums::Value::Success => Ok(f_out),
                err => Err(err),
            }
        }
    }

    /// This function returns the value of the n-th sample point in the unit interval,
    /// (j_{\nu,n+1}/j_{\nu,M}) X. These are the points where the function f(t) is assumed to be
    /// sampled.
    pub fn x_sample(&self, n: i32) -> f64 {
        unsafe { ffi::gsl_dht_x_sample(self.t, n) }
    }

    /// This function returns the value of the n-th sample point in “k-space”, j_{\nu,n+1}/X.
    pub fn k_sample(&self, n: i32) -> f64 {
        unsafe { ffi::gsl_dht_k_sample(self.t, n) }
    }
}

impl Drop for DiscreteHankel {
    fn drop(&mut self) {
        unsafe { ffi::gsl_dht_free(self.t) };
        self.t = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_dht> for DiscreteHankel {
    fn wrap(t: *mut ffi::gsl_dht) -> DiscreteHankel {
        DiscreteHankel {
            t: t
        }
    }

    fn soft_wrap(t: *mut ffi::gsl_dht) -> DiscreteHankel {
        Self::wrap(t)
    }

    fn unwrap_shared(t: &DiscreteHankel) -> *const ffi::gsl_dht {
        t.t as *const _
    }

    fn unwrap_unique(t: &mut DiscreteHankel) -> *mut ffi::gsl_dht {
        t.t
    }
}

// The following tests have been made and tested against the following C code:
//
// ```ignore
// #include <gsl/gsl_dht.h>
//
// int main() {
//     gsl_dht *t = gsl_dht_alloc(3);
//     printf("%d\n", gsl_dht_init(t, 3., 2.));
//     printf("%f %f\n", gsl_dht_x_sample(t, 1), gsl_dht_k_sample(t, 1));
//     double in[] = {100., 2., 3.};
//     double out[] = {0., 0., 0.};
//     gsl_dht_apply(t, in, out);
//     printf("%f %f %f\n", out[0], out[1], out[2]);
//     gsl_dht_free(t);
//     return 0;
// }
// ```
#[test]
fn discrete_hankel() {
    let mut d = DiscreteHankel::new(3).unwrap();
    assert_eq!(d.init(3., 2.), ::Value::Success);
    assert_eq!(&format!("{:.4} {:.4}", d.x_sample(1), d.k_sample(1)), "1.2033 4.8805");
    let v = d.apply(&[100., 2., 3.]);
    assert_eq!(true, v.is_ok());
    let v = v.unwrap();
    assert_eq!(&format!("{:.4} {:.4} {:.4}", v[0], v[1], v[2]), "8.5259 13.9819 11.7320");
}
