//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following example program computes Chebyshev approximations to a step function. This is an extremely difficult approximation to make,
// due to the discontinuity, and was chosen as an example where approximation error is visible. For smooth functions the Chebyshev
// approximation converges extremely rapidly and errors would not be visible.

extern crate rgsl;

use rgsl::ChebSeries;

#[allow(unused_variables)]
fn f(x: f64, param: &mut i32) -> f64 {
    if x < 0.5 {
        0.25
    } else {
        0.75
    }
}

fn main() {
    let n = 10000i32;
    let mut cs = ChebSeries::new(40).unwrap();

    cs.init(f, 0f64, 1f64, &mut 1i32);
    for i in range(0, n) {
        let x = i as f64 / n as f64;
        let r10 = cs.eval_n(10, x);
        let r40 = cs.eval(x);
        println!("{} {} {} {}", x, f(x, &mut 0i32), r10, r40);
    }
}