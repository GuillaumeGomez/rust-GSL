//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::bessel;

fn main() {
    let x = 5.;
    let y = bessel::J0(x);
    println!("J0({}) = {:.18}", x, y);
}
