//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::fit;

const N: usize = 4;

fn main() {
    let x = &[1970., 1980., 1990., 2000.];
    let y = &[12., 11., 14., 13.];
    let w = &[0.1, 0.2, 0.3, 0.4];

    let (_, c0, c1, cov00, cov01, cov11, chisq) = fit::wlinear(x, 1, w, 1, y, 1, N);

    println!("# best fit: Y = {} + {} X", c0, c1);
    println!("# covariance matrix:");
    println!("# [ {}, {}\n#   {}, {}]", cov00, cov01, cov01, cov11);
    println!("# chisq = {}", chisq);

    for i in 0..N {
        println!("data: {} {} {}", x[i], y[i], 1. / w[i].sqrt());
    }

    println!("");

    for i in -30..130 {
        let xf = x[0] + (i as f64 / 100.) * (x[N - 1] - x[0]);
        let (_, yf, yf_err) = fit::linear_est(xf, c0, c1, cov00, cov01, cov11);

        println!("fit: {} {}", xf, yf);
        println!("hi : {} {}", xf, yf + yf_err);
        println!("lo : {} {}", xf, yf - yf_err);
    }
}
