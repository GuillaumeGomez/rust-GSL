//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::Rng;

fn main() {
    rgsl::RngType::env_setup();

    let mu = 3f64;
    let t = rgsl::rng::default();
    let mut r = Rng::new(&t).unwrap();

    // print n random variates chosen from the poisson distribution with mean parameter mu
    for _ in 0u8..10u8 {
        print!("{} ", rgsl::randist::poisson::poisson(&mut r, mu));
    }
    println!("");
}