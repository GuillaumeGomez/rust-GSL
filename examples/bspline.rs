//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{multilinear, stats};
use rgsl::{BSpLineWorkspace, MatrixF64, MultifitLinearWorkspace, Rng, RngType, VectorF64};

const N: usize = 200;
const NCOEFFS: usize = 12;

fn main() {
    let rng_ty = RngType::env_setup().expect("Failed to setup RngType...");
    let mut r = Rng::new(rng_ty).expect("Rng::new failed...");

    // allocate a cubic bspline workspace (k = 4)
    let mut bw = BSpLineWorkspace::new(4, NCOEFFS - 2).expect("BSpLineWorkspace::new failed...");
    let mut b = VectorF64::new(NCOEFFS).expect("VectorF64::new failed...");

    let mut x = VectorF64::new(N).expect("VectorF64::new failed...");
    let mut y = VectorF64::new(N).expect("VectorF64::new failed...");
    let mut mat_x = MatrixF64::new(N, NCOEFFS).expect("MatrixF64::new failed...");
    let mut c = VectorF64::new(NCOEFFS).expect("VectorF64::new failed...");
    let mut w = VectorF64::new(N).expect("VectorF64::new failed...");
    let mut cov = MatrixF64::new(NCOEFFS, NCOEFFS).expect("MatrixF64::new failed...");
    let mut mw =
        MultifitLinearWorkspace::new(N, NCOEFFS).expect("MultifitLinearWorkspace::new failed...");

    // this is the data to be fitted
    for i in 0..N {
        let xi = (15. / (N - 1) as f64) * i as f64;
        let mut yi = xi.cos() * (-0.1 * xi).exp();
        let sigma = 0.1 * yi;

        let dy = r.gaussian(sigma);
        yi += dy;

        x.set(i, xi);
        y.set(i, yi);
        w.set(i, 1. / (sigma * sigma));

        println!("{} {}", xi, yi);
    }

    // use uniform breakpoints on [0, 15]
    bw.knots_uniform(0., 15.);

    // construct the fit matrix X
    for i in 0..N {
        let xi = x.get(i);

        // compute B_j(xi) for all j
        bw.eval(xi, &mut b);

        // fill in row i of X
        for j in 0..NCOEFFS {
            let bj = b.get(j);
            mat_x.set(i, j, bj);
        }
    }

    // do the fit
    let (chisq, _) = mw.wlinear(&mat_x, &w, &y, &mut c, &mut cov);

    let dof = N - NCOEFFS;
    let tss = stats::wtss(
        w.as_slice().expect("as_slice failed"),
        1,
        y.as_slice().expect("as_slice failed"),
        1,
    );
    let rsq = 1. - chisq / tss;

    eprintln!("chisq/dof = {}, rsq = {}", chisq / dof as f64, rsq);

    println!("");
    println!("");

    // output the smoothed curve
    let mut xi = 0.;
    while xi < 15. {
        bw.eval(xi, &mut b);
        let (yi, _, _) = multilinear::linear_est(&b, &c, &cov);
        println!("{} {}", xi, yi);
        xi += 0.1;
    }
}
