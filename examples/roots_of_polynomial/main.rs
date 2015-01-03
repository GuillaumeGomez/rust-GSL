//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// To demonstrate the use of the general polynomial solver we will take the polynomial P(x) = x^5 - 1 which has the following roots,
// 
// 1, e^{2\pi i /5}, e^{4\pi i /5}, e^{6\pi i /5}, e^{8\pi i /5}
// The following program will find these roots.

extern crate rgsl;

use rgsl::PolyComplex;

fn main() {
    let a : [f64; 6] = [-1f64, 0f64, 0f64, 0f64, 0f64, 1f64];
    let mut z : [f64; 10] = [0f64; 10];
    let w = PolyComplex::new(6).unwrap();

    w.solve(&a, &mut z);
    for i in range(0, 5) {
        println!("z{} = {} {}", i, z[2 * i], z[2 * i + 1]);
    }
}