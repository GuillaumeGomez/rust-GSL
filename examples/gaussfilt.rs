//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{filter, FilterEnd, FilterGaussianWorkspace, Rng, RngType, VectorF64};

const N: usize = 500; // length of time series
const K: usize = 51; // window size
const ALPHA: &[f64] = &[0.5, 3., 10.]; // alpha values

fn main() {
    // input vector
    let mut x = VectorF64::new(N).expect("VectorF64::new failed");
    // filtered output vector for alpha1
    let mut y1 = VectorF64::new(N).expect("VectorF64::new failed");
    // filtered output vector for alpha2
    let mut y2 = VectorF64::new(N).expect("VectorF64::new failed");
    // filtered output vector for alpha3
    let mut y3 = VectorF64::new(N).expect("VectorF64::new failed");
    // Gaussian kernel for alpha1
    let mut k1 = VectorF64::new(K).expect("VectorF64::new failed");
    // Gaussian kernel for alpha2
    let mut k2 = VectorF64::new(K).expect("VectorF64::new failed");
    // Gaussian kernel for alpha3
    let mut k3 = VectorF64::new(K).expect("VectorF64::new failed");

    let mut r = Rng::new(RngType::default()).expect("Rng::new failed");
    let mut gauss_p = FilterGaussianWorkspace::new(K).expect("FilterGaussianWorkspace::new failed");

    let mut sum = 0.;

    // generate input signal
    for i in 0..N {
        let ui = r.gaussian(1.);
        sum += ui;
        x.set(i, sum);
    }

    // compute kernels without normalization
    filter::gaussian_kernel(ALPHA[0], 0, false, &mut k1);
    filter::gaussian_kernel(ALPHA[1], 0, false, &mut k2);
    filter::gaussian_kernel(ALPHA[2], 0, false, &mut k3);

    // apply filters
    gauss_p.gaussian(FilterEnd::PadValue, ALPHA[0], 0, &x, &mut y1);
    gauss_p.gaussian(FilterEnd::PadValue, ALPHA[1], 0, &x, &mut y2);
    gauss_p.gaussian(FilterEnd::PadValue, ALPHA[2], 0, &x, &mut y3);

    // print kernels
    for i in 0..K {
        let k1i = k1.get(i);
        let k2i = k2.get(i);
        let k3i = k3.get(i);

        println!("{:.6} {:.6} {:.6}", k1i, k2i, k3i);
    }

    println!("");
    println!("");

    // print filter results
    for i in 0..N {
        let xi = x.get(i);
        let y1i = y1.get(i);
        let y2i = y2.get(i);
        let y3i = y3.get(i);

        println!("{:.6} {:.6} {:.6} {:.6}", xi, y1i, y2i, y3i);
    }
}
