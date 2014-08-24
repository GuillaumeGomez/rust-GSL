//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following example program computes Chebyshev approximations to a step function. This is an extremely difficult approximation to make,
// due to the discontinuity, and was chosen as an example where approximation error is visible. For smooth functions the Chebyshev
// approximation converges extremely rapidly and errors would not be visible.

extern crate rgsl;

fn main() {
    let mut data : [f64, ..256] = [0f64, ..256];

    data[0] = 1f64;
    for i in range(1u, 11u) {
        data[2 * i] = 1f64;
        data[(128 - i) * 2] = 1f64;
    }
    for i in range(0u, 128u) {
        println!("{} {:.6} {:.6}", i, data[2 * i], data[2 * i + 1]);
    }
    println!("");
    rgsl::radix2_fft::forward(data, 1, 128);
    for i in range(0u, 128u) {
        println!("{} {:.6} {:.6}", i, data[2 * i] / 128f64.sqrt(), data[2 * i + 1] / 128f64.sqrt());
    }
}