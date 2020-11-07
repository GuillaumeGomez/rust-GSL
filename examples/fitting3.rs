//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{Rng, RngType};

fn main() {
    RngType::env_setup();

    let mut r = Rng::new(RngType::default()).expect("Rng::new failed");
    let mut x: f64 = 0.1;

    while x < 2. {
        let y0 = x.exp();
        let sigma = 0.1 * y0;
        let dy = r.gaussian(sigma);

        println!("{:.1} {:.5} {:.5}", x, y0 + dy, sigma);
        x += 0.1;
    }
}
