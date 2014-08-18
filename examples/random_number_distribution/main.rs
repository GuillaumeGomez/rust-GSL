//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::Rng;

fn main() {
    rgsl::RngType::env_setup();

    let mu = 3f64;
    let t = rgsl::rng::default();
    let r = Rng::new(&t).unwrap();

    // print n random variates chosen from the poisson distribution with mean parameter mu
    for tmp in range(0i, 10i) {
        print!("{} ", rgsl::randist::poisson::poisson(&r, mu));
    }
    println!("");
}