//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following program prints the first 1024 points of the 2-dimensional Sobol sequence.

#![feature(core)]

extern crate rgsl;

fn main() {
    let q = rgsl::QRng::new(&rgsl::QRngType::sobol(), 2).unwrap();

    for i in range(0usize, 1024usize) {
        let mut v : [f64; 2] = [0f64, 0f64];

        q.get(&mut v);
        println!("{}: {:.5} {:.5}", i, v[0], v[1]);
    }
}