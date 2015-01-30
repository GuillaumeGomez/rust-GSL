//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// Here is an example program which computes the FFT of a short pulse in a sample of length 128. To make the resulting Fourier transform real
// the pulse is defined for equal positive and negative times (-10 … 10), where the negative times wrap around the end of the array.
//
// The second part is Here is an example which computes the FFT of a short pulse in a sample of length 630 (=2*3*3*5*7) using the mixed-radix algorithm.

#![feature(core)]

extern crate rgsl;

use std::num::Float;

fn main() {
    /* Part 1 */
    println!("=== PART 1 ===");
    let mut data : [f64; 256] = [0f64; 256];

    data[0] = 1f64;
    for i in range(1us, 11us) {
        data[2 * i] = 1f64;
        data[(128 - i) * 2] = 1f64;
    }
    for i in range(0us, 128us) {
        println!("{} {:.6} {:.6}", i, data[2 * i], data[2 * i + 1]);
    }
    println!("");
    rgsl::fft::radix2::forward(&mut data, 1, 128);
    for i in range(0us, 128us) {
        println!("{} {:.6} {:.6}", i, data[2 * i] / 128f64.sqrt(), data[2 * i + 1] / 128f64.sqrt());
    }

    /* Part 2 */
    println!("\n=== PART 2 ===");
    let mut data2 : [f64; 1260] = [0f64; 1260];
    let n = 630u64;

    for i in range(1us, 11us) {
        data2[2 * i] = 1f64;
        data2[(128 - i) * 2] = 1f64;
    }
    for i in range(0us, n as usize) {
        println!("{} {:.6} {:.6}", i, data2[2 * i], data2[2 * i + 1]);
    }
    let mut wavetable = rgsl::FftComplexWaveTable::new(n).unwrap();
    let workspace = rgsl::FftComplexWorkspace::new(n).unwrap();

    for i in range(0, wavetable.factor().len()) {
        println!("# factor {}: {}", i, wavetable.factor()[i]);
    }

    rgsl::fft::mixed_radix::forward(&mut data2, 1, n, &wavetable, &workspace);

    println!("");
    for i in range(0us, n as usize) {
        println!("{}: {:.6} {:.6}", i, data2[2 * i], data2[2 * i + 1]);
    }
}