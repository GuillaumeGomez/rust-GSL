extern crate rgsl;

use rgsl::gsl_multifit_nlinear_basic;
use rgsl::gsl_multifit_nlinear_basic_df;

use rgsl::types::rng::Rng;
use rgsl::types::rng::RngType;


fn expb_f(params: Vec<f64>, t: f64, _args: Vec<f64>) -> f64 {

    let a = params.get(0).unwrap();
    let lambda = params.get(1).unwrap();
    let b = params.get(2).unwrap();

    a * f64::exp(-lambda * t) + b
}

fn expb_df_a(params: Vec<f64>, t: f64, _args: Vec<f64>) -> f64 {
    let lambda = params.get(1).unwrap();
    f64::exp(-lambda * t)
}

fn expb_df_lambda(params: Vec<f64>, t: f64, _args: Vec<f64>) -> f64 {
    let a = params.get(0).unwrap();
    let lambda = params.get(1).unwrap();
    -t * a * f64::exp(-lambda * t)
}

fn expb_df_b(_params: Vec<f64>, _t: f64, _args: Vec<f64>) -> f64 {
    1.0
}

fn main() {

    let params_out = vec![1.0, 1.0, 0.0];
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    let args = vec![];
    let expb_dfs = vec![expb_df_a, expb_df_lambda, expb_df_b];

    let mut rng = Rng::new(RngType::default()).unwrap();

    for i in 0..100 {

        let rand_flt = rng.uniform();

        let ti = (i as f64) * 3.0 / (100.0 - 1.0);
        let yi = 1.0 + 5.0 * f64::exp(-1.5 * ti);
        let si = 0.1 * yi;
        let dy = si * rand_flt;

        ts.push(ti);
        ys.push(yi + dy);
    }

    let (params, parerr, status) = match unsafe {
        gsl_multifit_nlinear_basic(expb_f, params_out.clone(), &ts, &ys, &args, 100)
    } {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error encountered during solve!");
            return;
        }
    };

    println!("{:?}", params);
    println!("{:?}", parerr);
    println!("{}", status);

    let (params, parerr, status) = match unsafe {
        gsl_multifit_nlinear_basic_df(expb_f, &expb_dfs, params_out.clone(), &ts, &ys, &args, 100)
    } {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Error encountered during solve!");
            return;
        }
    };

    println!("{:?}", params);
    println!("{:?}", parerr);
    println!("{}", status);
}
