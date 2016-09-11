//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following example program fits a weighted exponential model with background to experimental data, Y = A \exp(-\lambda t) + b. The
// first part of the program sets up the functions expb_f and expb_df to calculate the model and its Jacobian. The appropriate fitting
// function is given by,
//
// f_i = ((A \exp(-\lambda t_i) + b) - y_i)/\sigma_i
// where we have chosen t_i = i. The Jacobian matrix J is the derivative of these functions with respect to the three parameters (A,
// \lambda, b). It is given by,
//
// J_{ij} = d f_i / d x_j
// where x_0 = A, x_1 = \lambda and x_2 = b.
//
// expfit.c -- model functions for exponential + background */

/*
The iteration terminates when the change in x is smaller than 0.0001, as both an absolute and relative change. Here are the results of
running the program:

iter: 0 x=1.00000000 0.00000000 0.00000000 |f(x)|=117.349
status=success
iter: 1 x=1.64659312 0.01814772 0.64659312 |f(x)|=76.4578
status=success
iter: 2 x=2.85876037 0.08092095 1.44796363 |f(x)|=37.6838
status=success
iter: 3 x=4.94899512 0.11942928 1.09457665 |f(x)|=9.58079
status=success
iter: 4 x=5.02175572 0.10287787 1.03388354 |f(x)|=5.63049
status=success
iter: 5 x=5.04520433 0.10405523 1.01941607 |f(x)|=5.44398
status=success
iter: 6 x=5.04535782 0.10404906 1.01924871 |f(x)|=5.44397
chisq/dof = 0.800996
A      = 5.04536 +/- 0.06028
lambda = 0.10405 +/- 0.00316
b      = 1.01925 +/- 0.03782
status = success
The approximate values of the parameters are found correctly, and the chi-squared value indicates a good fit (the chi-squared per degree
of freedom is approximately 1). In this case the errors on the parameters can be estimated from the square roots of the diagonal elements
of the covariance matrix.

If the chi-squared value shows a poor fit (i.e. chi^2/dof >> 1) then the error estimates obtained from the covariance matrix will be too
small. In the example program the error estimates are multiplied by \sqrt{\chi^2/dof} in this case, a common way of increasing the errors
for a poor fit. Note that a poor fit will result from the use an inappropriate model, and the scaled error estimates may then be outside
the range of validity for Gaussian errors.
*/

#![allow(non_snake_case)]

extern crate rgsl;

struct Data {
    n: usize,
    y: Vec<f64>,
    sigma: Vec<f64>
}

fn expb_f(x: &rgsl::VectorF64, data: &mut Data, f: &mut rgsl::VectorF64) -> rgsl::Value {
    let A = x.get(0);
    let lambda = x.get(1);
    let b = x.get(2);

    for i in 0..data.n {
        /* Model Yi = A * exp(-lambda * i) + b */
        let t = i as f64;
        let Yi = A * (-lambda * t).exp() + b;

        f.set(i, (Yi - data.y[i as usize]) / data.sigma[i as usize]);
    }

    rgsl::Value::Success
}

fn expb_df(x: &rgsl::VectorF64, data: &mut Data, J: &mut rgsl::MatrixF64) -> rgsl::Value {
    let A = x.get(0);
    let lambda = x.get(1);

    for i in 0..data.n {
        /* Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = (Yi - yi)/sigma[i],      */
        /*       Yi = A * exp(-lambda * i) + b  */
        /* and the xj are the parameters (A,lambda,b) */
        let t = i as f64;
        let s = data.sigma[i as usize];
        let e = (-lambda * t).exp();

        J.set(i, 0, e / s);
        J.set(i, 1, -t * A * e / s);
        J.set(i, 2, 1f64 / s);
    }
    rgsl::Value::Success
}

fn expb_fdf(x: &rgsl::VectorF64, data: &mut Data, f: &mut rgsl::VectorF64, J: &mut rgsl::MatrixF64) -> rgsl::Value {
    expb_f(x, data, f);
    expb_df(x, data, J);

    rgsl::Value::Success
}

// The main part of the program sets up a Levenberg-Marquardt solver and some simulated random data. The data uses the known parameters
// (1.0,5.0,0.1) combined with Gaussian noise (standard deviation = 0.1) over a range of 40 timesteps. The initial guess for the
// parameters is chosen as (0.0, 1.0, 0.0).

static N : usize = 40usize;

#[allow(unused_assignments)]
fn main() {
    let mut status = rgsl::Value::Success;
    let n = N;
    let p = 3;

    let mut covar = rgsl::MatrixF64::new(p, p).unwrap();
    let mut d = Data {
        n: n,
        y: ::std::iter::repeat(0f64).take(N).collect(),
        sigma: ::std::iter::repeat(0f64).take(N).collect()
    };
    let mut x_init : [f64; 3] = [1f64, 0f64, 0f64];
    let mut x = rgsl::VectorView::from_array(&mut x_init);

    rgsl::RngType::env_setup();
    let t : rgsl::RngType = rgsl::rng::default();
    let r = rgsl::Rng::new(&t).unwrap();
    let T = rgsl::MultiFitFdfSolverType::lmsder();

    let mut f = rgsl::MultiFitFunctionFdf {
        f: expb_f,
        df: Some(expb_df),
        fdf: Some(expb_fdf),
        n: n,
        p: p,
        params: &mut d
    };

    /* This is the data to be fitted */
    let mut iter = 0;

    for i in 0..n {
        let t = i as f64;

        f.params.y[i] = 1f64 + 5f64 * (-0.1f64 * t).exp() + rgsl::randist::gaussian::gaussian(&r, 0.1f64);
        f.params.sigma[i] = 0.1f64;
        println!("data: {:2} {:.5} {:.5}", i, f.params.y[i], f.params.sigma[i]);
    }

    let mut s = rgsl::MultiFitFdfSolver::new(&T, n, p).unwrap();

    s.set(&mut f, &x.vector());

    loop {
        iter += 1;
        status = s.iterate();

        println!("status = {}", rgsl::error::str_error(status));

        if status != rgsl::Value::Success {
            break;
        }

        status = rgsl::multifit::test_delta(&s.dx, &s.x, 1e-4, 1e-4);
        if status != rgsl::Value::Continue || iter >= 500 {
            break;
        }
    }

    rgsl::multifit::covar(&s.j, 0f64, &mut covar);

    {
        let chi = rgsl::blas::level1::dnrm2(&s.f);
        let dof = n as f64 - p as f64;
        let c = 1f64.max(chi / dof.sqrt());

        println!("chisq/dof = {}",  chi.powf(2f64) / dof);

        println!("A      = {:.5} +/- {:.5}", s.x.get(0), c * covar.get(0, 0).sqrt());
        println!("lambda = {:.5} +/- {:.5}", s.x.get(1), c * covar.get(1, 1).sqrt());
        println!("b      = {:.5} +/- {:.5}", s.x.get(2), c * covar.get(2, 2).sqrt());
    }

    println!("status = {}", rgsl::error::str_error(status));
}
