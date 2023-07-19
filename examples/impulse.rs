//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{FilterEnd, FilterImpulseWorkspace, FilterScale, Rng, RngType, VectorF64, VectorI32};
use std::f64::consts::PI;

const N: usize = 1000; // length of time series
const K: usize = 25; // window size
const T: f64 = 4.; // number of scale factors for outlier detection

fn main() {
    // input vector
    let mut x = VectorF64::new(N).expect("VectorF64::new failed");
    // filtered output vector
    let mut y = VectorF64::new(N).expect("VectorF64::new failed");
    // window medians
    let mut xmedian = VectorF64::new(N).expect("VectorF64::new failed");
    // window scale estimates
    let mut xsigma = VectorF64::new(N).expect("VectorF64::new failed");
    // outlier detected?
    let mut ioutlier = VectorI32::new(N).expect("VectorF64::new failed");
    let mut w = FilterImpulseWorkspace::new(K).expect("FilterImpulseWorkspace::new failed");
    let mut r = Rng::new(RngType::default()).expect("Rng::new failed");

    // generate input signal
    for i in 0..N {
        let xi = 10. * (2. * PI * i as f64 / N as f64).sin();
        let ei = r.gaussian(2.);
        let u = r.uniform();
        let outlier = if u < 0.001 {
            15. * if ei >= 0. { 1. } else { -1. }
        } else {
            0.
        };

        x.set(i, xi + ei + outlier);
    }

    // apply impulse detection filter
    w.impulse(
        FilterEnd::Truncate,
        FilterScale::QN,
        T,
        &x,
        &mut y,
        &mut xmedian,
        &mut xsigma,
        &mut ioutlier,
    )
    .unwrap();

    for i in 0..N {
        let xi = x.get(i);
        let yi = y.get(i);
        let xmedi = xmedian.get(i);
        let xsigmai = xsigma.get(i);
        let outlier = ioutlier.get(i);

        println!(
            "{} {} {} {} {} {}",
            i,
            xi,
            yi,
            xmedi + T * xsigmai,
            xmedi - T * xsigmai,
            outlier
        );
    }
}
