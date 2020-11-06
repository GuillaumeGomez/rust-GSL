//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{multifit, MatrixF64, MultifitLinearWorkspace, Rng, RngType, VectorF64};

const N: usize = 1000; // number of observations
const P: usize = 2; // number of model parameters
const NPOINTS: usize = 200; // number of points on L-curve and GCV curve

fn main() {
    let mut r = Rng::new(RngType::default()).expect("Rng::new failed");
    let mut x = MatrixF64::new(N, P).expect("MatrixF64::new");
    let mut y = VectorF64::new(N).expect("MatrixF64::new");

    for i in 0..N {
        // generate first random variable u
        let ui = 5. * r.gaussian(1.);

        // set v = u + noise
        let vi = ui + r.gaussian(0.001);

        // set y = u + v + noise
        let yi = ui + vi + r.gaussian(1.);

        // since u =~ v, the matrix X is ill-conditioned
        x.set(i, 0, ui);
        x.set(i, 1, vi);

        // rhs vector
        y.set(i, yi);
    }

    let mut w = MultifitLinearWorkspace::new(N, P).expect("MultifitLinearWorkspace::new failed");
    // OLS solution
    let mut c = VectorF64::new(P).expect("VectorF64::new");
    // regularized solution (L-curve)
    let mut c_lcurve = VectorF64::new(P).expect("VectorF64::new");
    // regularized solution (GCV)
    let mut c_gcv = VectorF64::new(P).expect("VectorF64::new");

    let mut reg_param = VectorF64::new(NPOINTS).expect("VectorF64::new");
    // residual norms
    let mut rho = VectorF64::new(NPOINTS).expect("VectorF64::new");
    // solution norms
    let mut eta = VectorF64::new(NPOINTS).expect("VectorF64::new");
    // GCV function values
    let mut g = VectorF64::new(NPOINTS).expect("VectorF64::new");

    // compute SVD of X
    w.linear_svd(&mut x);

    // Get reciprocal condition number of X
    let rcond = w.linear_rcond();
    eprintln!("matrix condition number = {}", 1. / rcond);
    eprintln!("");

    // unregularized (standard) least squares fit, lambda = 0
    let (rnorm, snorm, _) = w.linear_solve(0., &x, &y, &mut c);
    let chisq = rnorm.powi(2);

    eprintln!("\n=== Unregularized fit ===");
    eprintln!("best fit: y = {} u + {} v", c.get(0), c.get(1));
    eprintln!("residual norm = {}", rnorm);
    eprintln!("solution norm = {}", snorm);
    eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

    // calculate L-curve and find its corner
    w.linear_lcurve(&y, &mut reg_param, &mut rho, &mut eta);
    let (reg_idx, _) = multifit::linear_lcorner(&rho, &eta);

    // store optimal regularization parameter
    let lambda_l = reg_param.get(reg_idx);

    // regularize with lambda_l
    let (rnorm, snorm, _) = w.linear_solve(lambda_l, &x, &y, &mut c);
    let chisq = rnorm.powi(2) + (lambda_l * snorm).powi(2);

    eprintln!("\n=== Regularized fit (L-curve) ===");
    eprintln!("optimal lambda: {}", lambda_l);
    eprintln!(
        "best fit: y = {} u + {} v",
        c_lcurve.get(0),
        c_lcurve.get(1)
    );
    eprintln!("residual norm = {}", rnorm);
    eprintln!("solution norm = {}", snorm);
    eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

    // calculate GCV curve and find its minimum
    let (lambda_gcv, g_gcv, _) = w.linear_gcv(&y, &mut reg_param, &mut g);

    // regularize with lambda_gcv
    let (rnorm, snorm, _) = w.linear_solve(lambda_gcv, &x, &y, &mut c_gcv);
    let chisq = rnorm.powi(2) + (lambda_gcv * snorm).powi(2);

    eprintln!("\n=== Regularized fit (GCV) ===\n");
    eprintln!("optimal lambda: {}", lambda_gcv);
    eprintln!("best fit: y = {} u + {} v", c_gcv.get(0), c_gcv.get(1));
    eprintln!("residual norm = {}", rnorm);
    eprintln!("solution norm = {}", snorm);
    eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

    // output L-curve and GCV curve
    for i in 0..NPOINTS {
        println!(
            "{:.6} {:.6} {:.6} {:.6}",
            reg_param.get(i),
            rho.get(i),
            eta.get(i),
            g.get(i)
        );
    }

    // output L-curve corner point
    println!("\n\n{:.6} {:.6}", rho.get(reg_idx), eta.get(reg_idx));

    // output GCV curve corner minimum
    println!("\n\n{:.6} {:.6}", lambda_gcv, g_gcv);
}
