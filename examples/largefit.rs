//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{
    blas, error, MatrixF64, MultilargeLinearType, MultilargeLinearWorkspace, Rng, RngType,
    VectorF64,
};

// number of observations
const N: usize = 50000;
// polynomial order + 1
const P: usize = 16;
// regularization parameter
const LAMBDA: f64 = 0.;
// number of blocks to accumulate
const NBLOCK: usize = 5;
// number of rows per block
const NROWS: usize = N / NBLOCK;
const NLCURVE: usize = 200;
const DT: f64 = 1. / (N as f64 - 1.);

fn func(t: f64) -> f64 {
    let x = (10. * t).sin();
    (x * x * x).exp()
}

fn build_row(t: f64, row: &mut VectorF64) {
    let p = row.len();
    let mut xj = 1.;

    for j in 0..p {
        row.set(j, xj);
        xj *= t;
    }
}

fn solve_system(print_data: bool, t: MultilargeLinearType, c: &mut VectorF64) {
    let mut w =
        MultilargeLinearWorkspace::new(t, P).expect("MultilargeLinearWorkspace::new failed");
    let mut x = MatrixF64::new(NROWS, P).expect("MatrixF64::new failed");
    let mut y = VectorF64::new(NROWS).expect("VectorF64::new failed");
    let mut r = Rng::new(RngType::default()).expect("Rng::new failed");

    let mut reg_param = VectorF64::new(NLCURVE).expect("VectorF64::new failed");
    let mut rho = VectorF64::new(NLCURVE).expect("VectorF64::new failed");
    let mut eta = VectorF64::new(NLCURVE).expect("VectorF64::new failed");

    let mut rowidx = 0;
    let mut t = 0.;

    while rowidx < N {
        // number of rows left to accumulate
        let nleft = N - rowidx;
        // number of rows in this block
        let nr = if NROWS > nleft { nleft } else { NROWS };

        let mut xv = x.submatrix(0, 0, nr, P);
        let mut yv = y.subvector(0, nr);

        // build (X,y) block with 'nr' rows
        for i in 0..nr {
            xv.matrix_mut(|mat| {
                mat.expect("Failed to get matrix").row(i, |row| {
                    let mut row = row.expect("Failed to get row...");
                    let fi = func(t);
                    // noise
                    let ei = r.gaussian(0.1 * fi);
                    let yi = fi + ei;

                    // construct this row of LS matrix
                    row.vector_mut(|vector| {
                        build_row(t, vector.expect("Failed to get vector"));
                    });

                    // set right hand side value with added noise
                    yv.vector_mut(|vector| {
                        vector.expect("Failed to get vector").set(i, yi);
                    });

                    if print_data && i % 100 == 0 {
                        println!("{} {}", t, yi);
                    }
                });

                t += DT;
            });
        }

        // accumulate (X,y) block into LS system
        xv.matrix_mut(|matrix| {
            yv.vector_mut(|vector| {
                w.accumulate(
                    matrix.expect("Failed to get matrix"),
                    vector.expect("Failed to get vector"),
                );
            });
        });

        rowidx += nr;
    }

    if print_data {
        println!("");
        println!("");
    }

    // compute L-curve
    w.lcurve(&mut reg_param, &mut rho, &mut eta);

    // solve large LS system and store solution in c
    let (_, rnorm, snorm) = w.solve(LAMBDA, c);

    // compute reciprocal condition number
    let (_, rcond) = w.rcond();

    eprintln!("=== Method {} ===\n", w.name().expect("Failed to get name"));
    eprintln!("condition number = {}", 1. / rcond);
    eprintln!("residual norm    = {}", rnorm);
    eprintln!("solution norm    = {}", snorm);

    // output L-curve
    for i in 0..NLCURVE {
        println!(
            "{:.12} {:.12} {:.12}",
            reg_param.get(i),
            rho.get(i),
            eta.get(i)
        );
    }
}

fn main() {
    let mut c_tsqr = VectorF64::new(P).expect("VectorF64::new failed");
    let mut c_normal = VectorF64::new(P).expect("VectorF64::new failed");

    // turn off error handler so normal equations method won't abort
    error::set_error_handler_off();

    // solve system with TSQR method
    solve_system(true, MultilargeLinearType::tsqr(), &mut c_tsqr);

    println!("");
    println!("");

    // solve system with Normal equations method
    solve_system(false, MultilargeLinearType::normal(), &mut c_normal);

    // output solution
    let mut v = VectorF64::new(P).expect("VectorF64::new failed");
    let mut t = 0.;

    while t <= 1. {
        let f_exact = func(t);
        build_row(t, &mut v);

        let (_, f_tsqr) = blas::level1::ddot(&v, &c_tsqr);
        let (_, f_normal) = blas::level1::ddot(&v, &c_normal);

        println!("{} {:.6} {:.6} {:.6}", t, f_exact, f_tsqr, f_normal);

        t += 0.01;
    }
}
