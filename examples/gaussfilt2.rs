//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{FilterEnd, FilterGaussianWorkspace, Rng, RngType, VectorF64};

const N: usize = 1000; // length of time series
const K: usize = 61; // window size
const ALPHA: f64 = 3.; // Gaussian kernel has +/- 3 standard deviations

fn main() {
    // input vector
    let mut x = VectorF64::new(N).expect("VectorF64::new failed");
    // filtered output vector
    let mut y = VectorF64::new(N).expect("VectorF64::new failed");
    // first derivative filtered vector
    let mut dy = VectorF64::new(N).expect("VectorF64::new failed");
    // second derivative filtered vector
    let mut d2y = VectorF64::new(N).expect("VectorF64::new failed");
    let mut r = Rng::new(RngType::default()).expect("Rng::new failed");
    let mut gauss_p = FilterGaussianWorkspace::new(K).expect("FilterGaussianWorkspace::new failed");

    // generate input signal
    for i in 0..N {
        let xi = if i > N / 2 { 0.5 } else { 0. };
        let ei = r.gaussian(0.1);

        x.set(i, xi + ei);
    }

    // apply filters
    gauss_p
        .gaussian(FilterEnd::PadValue, ALPHA, 0, &x, &mut y)
        .unwrap();
    gauss_p
        .gaussian(FilterEnd::PadValue, ALPHA, 1, &x, &mut dy)
        .unwrap();
    gauss_p
        .gaussian(FilterEnd::PadValue, ALPHA, 2, &x, &mut d2y)
        .unwrap();

    // print results
    for i in 0..N {
        let xi = x.get(i);
        let yi = y.get(i);
        let dyi = dy.get(i);
        let d2yi = d2y.get(i);
        let dxi = if i == 0 {
            x.get(i + 1) - xi
        } else if i == N - 1 {
            x.get(i) - x.get(i - 1)
        } else {
            0.5 * (x.get(i + 1) - x.get(i - 1))
        };

        println!("{:.12} {:.12} {:.12} {:.12} {:.12}", xi, yi, dxi, dyi, d2yi);
    }
}
