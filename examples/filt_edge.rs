//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

#[cfg(feature = "v2_5")]
mod example {
    use rgsl::{FilterEnd, FilterMedianWorkspace, FilterRMedianWorkspace, Rng, RngType, VectorF64};
    use std::f64::consts::PI;

    const N: usize = 1000; // length of time series
    const K: usize = 7; // window size
    const F: f64 = 5.0; // frequency of square wave in Hz

    pub fn run() {
        let mut median_p = FilterMedianWorkspace::new(K).expect("FilterMedianWorkspace::new");
        let mut rmedian_p = FilterRMedianWorkspace::new(K).expect("FilterRMedianWorkspace::new");
        // time
        let mut t = VectorF64::new(N).expect("VectorF64::new failed");
        // input vector
        let mut x = VectorF64::new(N).expect("VectorF64::new failed");
        // median filtered output
        let mut y_median = VectorF64::new(N).expect("VectorF64::new failed");
        // recursive median filtered output
        let mut y_rmedian = VectorF64::new(N).expect("VectorF64::new failed");
        let mut r = Rng::new(RngType::default()).expect("Rng::new failed");

        // generate input signal
        for i in 0..N {
            let ti = i as f64 / (N as f64 - 1.);
            let tmp = (2. * PI * F * ti).sin();
            let xi = if tmp >= 0. { 1. } else { -1. };
            let ei = r.gaussian(0.1);

            t.set(i, ti);
            x.set(i, xi + ei);
        }

        median_p
            .median(FilterEnd::PadValue, &x, &mut y_median)
            .unwrap();
        rmedian_p
            .rmedian(FilterEnd::PadValue, &x, &mut y_rmedian)
            .unwrap();

        // print results
        for i in 0..N {
            let ti = t.get(i);
            let xi = x.get(i);
            let medi = y_median.get(i);
            let rmedi = y_rmedian.get(i);

            println!("{:.6} {:.6} {:.6} {:.6}", ti, xi, medi, rmedi);
        }
    }
}

#[cfg(feature = "v2_5")]
fn main() {
    example::run();
}

#[cfg(not(feature = "v2_5"))]
fn main() {
    eprintln!("You need to enable the `v2_5` feature to be able to run this example");
}
