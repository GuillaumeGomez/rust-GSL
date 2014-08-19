//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::RngType;

fn main() {
    rgsl::RngType::env_setup();

    let t : RngType = rgsl::rng::default();
    let r = rgsl::Rng::new(&t).unwrap();
    let k = 5;
    let n = 100000;
    let mut x : [f64, ..100000] = [0f64, ..100000];
    let mut small : [f64, ..5] = [0f64, ..5];

    for tmp in range(0, n) {
        x[tmp] = r.uniform();
    }

    rgsl::sort::select::sort_smallest(small, k, x, 1);
    println!("{} smallest values from {}", k, n);
    for tmp in range(0, k as uint) {
        println!("{}: {}", k, small[tmp]);
    }
}