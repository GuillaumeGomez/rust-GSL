//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*
The following code calculates an estimate of \zeta(2) = \pi^2 / 6 using the series,

\zeta(2) = 1 + 1/2^2 + 1/3^2 + 1/4^2 + ...
After N terms the error in the sum is O(1/N), making direct summation of the series converge slowly.
*/

extern crate rgsl;

pub static N : usize = 20usize;

fn main() {
    let mut t : [f64; 20] = [0f64; 20];
    let mut sum_accel = 0f64;
    let mut err = 0f64;
    let mut sum = 0f64;

    let w = rgsl::LevinUWorkspace::new(N).unwrap();

    let zeta_2 = ::std::f64::consts::PI * ::std::f64::consts::PI / 6f64;

    /* terms for zeta(2) = \sum_{n=1}^{\infty} 1/n^2 */

    for n in 0usize..N {
        let np1 = n as f64 + 1f64;

        t[n] = 1f64 / (np1 * np1);
        sum += t[n];
    }

    w.accel(&t, &mut sum_accel, &mut err);

    println!("term-by-term sum = {:.16} using {} terms", sum, N);

    println!("term-by-term sum = {:.16} using {} terms", w.sum_plain(), w.terms_used());

    println!("exact value      = {:.16}", zeta_2);
    println!("accelerated sum  = {:.16} using {} terms", sum_accel, w.terms_used());

    println!("estimated error  = {:.16}", err);
    println!("actual error     = {:.16}", sum_accel - zeta_2);
}