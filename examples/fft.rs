//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use num_complex::c64;
use rgsl::fft;

const N: usize = 128;

fn main() {
    let data = &mut [c64(0., 0.); N];

    data[0].re = 1.;

    for i in 1..=10 {
        data[i].re = 1.;
        data[N - i].re = 1.;
    }

    for i in 0..N {
        println!("{} {}", i, data[i]);
    }
    println!();
    println!();

    fft::radix2::forward(data).unwrap();

    for i in 0..N {
        println!("{} {}", i, data[i] / 128f64.sqrt());
    }
}
