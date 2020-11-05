//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::randist::gaussian;

fn main() {
    let x = 2.;

    let p = gaussian::ugaussian_P(x);
    println!("prob(x < {}) = {}", x, p);

    let q = gaussian::ugaussian_Q(x);
    println!("prob(x > {}) = {}", x, q);

    let x = gaussian::ugaussian_Pinv(p);
    println!("Pinv({}) = {}", p, x);

    let x = gaussian::ugaussian_Qinv(q);
    println!("Pinv({}) = {}", q, x);
}
