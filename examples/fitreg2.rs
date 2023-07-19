//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{multifit, MatrixF64, MultifitLinearWorkspace, VectorF64};

const N: usize = 10; // number of observations
const P: usize = 8; // number of model parameters
const NPOINTS: usize = 200; // number of points on L-curve and GCV curve

fn hibert_matrix() -> Option<MatrixF64> {
    let mut x = MatrixF64::new(N, P)?;
    let n = x.size1();
    let m = x.size2();

    for i in 0..n {
        for j in 0..m {
            x.set(i, j, 1. / ((i + j) as f64 + 1.));
        }
    }
    Some(x)
}

fn main() {
    let mut y = VectorF64::new(N).expect("MatrixF64::new failed");

    // construct Hilbert matrix and rhs vector
    let mut x = hibert_matrix().expect("hibert_matrix failed");

    let mut val = 1.;
    for i in 0..N {
        y.set(i, val);
        val *= -1.;
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
    w.linear_svd(&mut x).unwrap();

    let rcond = w.linear_rcond();
    eprintln!("matrix condition number = {}", 1. / rcond);
    eprintln!("");

    // unregularized (standard) least squares fit, lambda = 0
    let (rnorm, snorm) = w.linear_solve(0., &x, &y, &mut c).unwrap();
    let chisq = rnorm.powi(2);

    eprintln!("\n=== Unregularized fit ===");
    eprintln!("residual norm = {}", rnorm);
    eprintln!("solution norm = {}", snorm);
    eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

    // calculate L-curve and find its corner
    w.linear_lcurve(&y, &mut reg_param, &mut rho, &mut eta)
        .unwrap();
    let reg_idx = multifit::linear_lcorner(&rho, &eta).unwrap();

    // store optimal regularization parameter
    let lambda_l = reg_param.get(reg_idx);

    // regularize with lambda_l
    let (rnorm, snorm) = w.linear_solve(lambda_l, &x, &y, &mut c_lcurve).unwrap();
    let chisq = rnorm.powi(2) + (lambda_l * snorm).powi(2);

    eprintln!("\n=== Regularized fit (L-curve) ===");
    eprintln!("optimal lambda: {}", lambda_l);
    eprintln!("residual norm = {}", rnorm);
    eprintln!("solution norm = {}", snorm);
    eprintln!("chisq/dof = {}", chisq / (N - P) as f64);

    // calculate GCV curve and find its minimum
    let (lambda_gcv, g_gcv) = w.linear_gcv(&y, &mut reg_param, &mut g).unwrap();

    // regularize with lambda_gcv
    let (rnorm, snorm) = w.linear_solve(lambda_gcv, &x, &y, &mut c_gcv).unwrap();
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
