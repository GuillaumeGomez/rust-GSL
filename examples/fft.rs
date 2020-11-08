//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::fft;

macro_rules! real {
    ($z:ident, $i:expr) => {
        $z[2 * ($i)]
    };
}
macro_rules! imag {
    ($z:ident, $i:expr) => {
        $z[2 * ($i) + 1]
    };
}

const N: usize = 128;

fn main() {
    let data = &mut [0.; 2 * N];

    real!(data, 0) = 1.;

    for i in 1..=10 {
        real!(data, i) = 1.;
        real!(data, N - i) = 1.;
    }

    for i in 0..N {
        println!("{} {} {}", i, real!(data, i), imag!(data, i));
    }
    println!("");
    println!("");

    fft::radix2::forward(data, 1, N);

    for i in 0..N {
        println!(
            "{} {} {}",
            i,
            real!(data, i) / 128f64.sqrt(),
            imag!(data, i) / 128f64.sqrt()
        );
    }
}
