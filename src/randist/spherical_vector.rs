//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The spherical distributions generate random vectors, located on a spherical surface.
They can be used as random directions, for example in the steps of a random walk.
!*/

use ffi;
use types::Rng;

/// This function returns a random direction vector v = (x,y) in two dimensions. The vector is normalized such that |v|^2 = x^2 + y^2 = 1.
/// The obvious way to do this is to take a uniform random number between 0 and 2\pi and let x and y be the sine and cosine respectively.
/// Two trig functions would have been expensive in the old days, but with modern hardware implementations, this is sometimes the fastest way to go.
/// This is the case for the Pentium (but not the case for the Sun Sparcstation).
/// One can avoid the trig evaluations by choosing x and y in the interior of a unit circle (choose them at random from the interior of the enclosing square,
/// and then reject those that are outside the unit circle), and then dividing by \sqrt{x^2 + y^2}. A much cleverer approach, attributed to von Neumann
/// (See Knuth, v2, 3rd ed, p140, exercise 23), requires neither trig nor a square root. In this approach, u and v are chosen at random from the interior of
/// a unit circle, and then x=(u^2-v^2)/(u^2+v^2) and y=2uv/(u^2+v^2).
pub fn dir_2d(r: &mut Rng, x: &mut f64, y: &mut f64) {
    unsafe { ffi::randist::gsl_ran_dir_2d(ffi::FFI::unwrap_unique(r), x, y) }
}

/// This function returns a random direction vector v = (x,y) in two dimensions. The vector is normalized such that |v|^2 = x^2 + y^2 = 1.
/// The obvious way to do this is to take a uniform random number between 0 and 2\pi and let x and y be the sine and cosine respectively.
/// Two trig functions would have been expensive in the old days, but with modern hardware implementations, this is sometimes the fastest way to go.
/// This is the case for the Pentium (but not the case for the Sun Sparcstation).
/// One can avoid the trig evaluations by choosing x and y in the interior of a unit circle (choose them at random from the interior of the enclosing square,
/// and then reject those that are outside the unit circle), and then dividing by \sqrt{x^2 + y^2}. A much cleverer approach, attributed to von Neumann
/// (See Knuth, v2, 3rd ed, p140, exercise 23), requires neither trig nor a square root. In this approach, u and v are chosen at random from the interior of
/// a unit circle, and then x=(u^2-v^2)/(u^2+v^2) and y=2uv/(u^2+v^2).
pub fn dir_2d_trig_method(r: &mut Rng, x: &mut f64, y: &mut f64) {
    unsafe { ffi::randist::gsl_ran_dir_2d_trig_method(ffi::FFI::unwrap_unique(r), x, y) }
}

/// This function returns a random direction vector v = (x,y,z) in three dimensions. The vector is normalized such that |v|^2 = x^2 + y^2 + z^2 = 1.
/// The method employed is due to Robert E. Knop (CACM 13, 326 (1970)), and explained in Knuth, v2, 3rd ed, p136. It uses the surprising fact that the
/// distribution projected along any axis is actually uniform (this is only true for 3 dimensions).
pub fn dir_3d(r: &mut Rng, x: &mut f64, y: &mut f64, z: &mut f64) {
    unsafe { ffi::randist::gsl_ran_dir_3d(ffi::FFI::unwrap_unique(r), x, y, z) }
}

/// This function returns a random direction vector v = (x_1,x_2,...,x_n) in n dimensions. The vector is normalized such that |v|^2 = x_1^2 + x_2^2 + ... + x_n^2 = 1.
/// The method uses the fact that a multivariate Gaussian distribution is spherically symmetric. Each component is generated to have a Gaussian distribution, and then
/// the components are normalized. The method is described by Knuth, v2, 3rd ed, p135–136, and attributed to G. W. Brown, Modern Mathematics for the Engineer (1956).
pub fn dir_nd(r: &mut Rng, x: &mut [f64]) {
    unsafe {
        ffi::randist::gsl_ran_dir_nd(ffi::FFI::unwrap_unique(r), x.len() as usize, x.as_mut_ptr())
    }
}
