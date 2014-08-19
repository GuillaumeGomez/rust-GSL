//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{Permutation, RngType};

fn main() {
    rgsl::RngType::env_setup();

    let mut p = Permutation::new(10).unwrap();
    let q = Permutation::new(10).unwrap();
    let t : RngType = rgsl::rng::default();
    let r = rgsl::Rng::new(&t).unwrap();

    println!("initial permutation :");
    p.init();
    println!("{}\n", p);

    println!("random permutation :");
    rgsl::randist::shuffling_sampling::shuffle(&r, p.data());
    println!("{}\n", p);

    println!("inverse permutation :");
    p.inverse(&q);
    println!("{}\n", q);
}