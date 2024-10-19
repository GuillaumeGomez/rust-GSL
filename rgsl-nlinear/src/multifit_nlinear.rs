#[allow(improper_ctypes)]
#[link(name = "rgslmfnlin")]
extern "C" {
    fn run_gsl_multifit_nlinear(
        func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
        params: *const f64,
        covars: *const f64,
        params_len: usize,
        ts: *const f64,
        ys: *const f64,
        vars_len: usize,
        args: *const f64,
        args_len: usize,
        max_iters: u64,
        status: *const i8
    );
}

#[allow(improper_ctypes)]
#[link(name = "rgslmfnlin")]
extern "C" {
    fn run_gsl_multifit_nlinear_df(
        func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
        func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
        params: *const f64,
        parerr: *const f64,
        params_len: usize,
        ts: *const f64,
        ys: *const f64,
        vars_len: usize,
        args: *const f64,
        args_len: usize,
        max_iters: u64,
        status: *const i8
    );
}

pub unsafe fn gsl_multifit_nlinear_basic(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    params_in: Vec<f64>,
    ts: Vec<f64>,
    ys: Vec<f64>,
    args: Vec<f64>,
    max_iters: u64
) -> (Vec<f64>, Vec<f64>, i8) {

    if ts.len() != ys.len() {
        eprintln!("Time length does not match Ys length!");
        return (vec![], vec![], -1);
    }

    let mut params: Vec<f64> = params_in.clone();
    let mut parerr: Vec<f64> = Vec::with_capacity(params_in.len());
    let mut status: i8 = 0;

    parerr.resize(params_in.len(), 0.0);

    unsafe {
        run_gsl_multifit_nlinear(
            func_f,
            params.as_mut_ptr(),
            parerr.as_mut_ptr(),
            params.len(),
            ts.as_ptr(),
            ys.as_ptr(),
            ts.len(),
            args.as_ptr(),
            args.len(),
            max_iters,
            &mut status
        );
    }

    (params, parerr, status)
}

pub unsafe fn gsl_multifit_nlinear_basic_df(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
    params_in: Vec<f64>,
    ts: Vec<f64>,
    ys: Vec<f64>,
    args: Vec<f64>,
    max_iters: u64
) -> (Vec<f64>, Vec<f64>, i8) {

    if ts.len() != ys.len() {
        eprintln!("Time length does not match Ys length!");
        return (vec![], vec![], -1);
    }

    let mut params: Vec<f64> = params_in.clone();
    let mut parerr: Vec<f64> = Vec::with_capacity(params_in.len());
    let mut status: i8 = 0;

    parerr.resize(params_in.len(), 0.0);

    unsafe {
        run_gsl_multifit_nlinear_df(
            func_f,
            func_dfs,
            params.as_mut_ptr(),
            parerr.as_mut_ptr(),
            params.len(),
            ts.as_ptr(),
            ys.as_ptr(),
            ts.len(),
            args.as_ptr(),
            args.len(),
            max_iters,
            &mut status
        );
    }

    (params, parerr, status)
}
