//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{Histogram2D, Histogram2DPdf, Rng, RngType};

fn main() {
    let mut h = Histogram2D::new(10, 10).expect("Histogram2D::new failed");

    h.set_ranges_uniform(0., 1., 0., 1.);

    h.accumulate(0.3, 0.3, 1.);
    h.accumulate(0.8, 0.1, 5.);
    h.accumulate(0.7, 0.9, 0.5);

    RngType::env_setup();

    let mut r = Rng::new(RngType::default()).expect("Rng::new failed");

    let mut p = Histogram2DPdf::new(h.nx(), h.ny()).expect("Histogram2DPdf::new failed");
    p.init(&mut h);

    for _ in 0..1000 {
        let u = r.uniform();
        let v = r.uniform();

        let (x, y, _) = p.sample(u, v);

        println!("{} {}", x, y);
    }
}
