use gsl_sys::gsl_vector;
use gsl_sys::gsl_vector_get;


#[no_mangle]
pub extern "C" fn rust_callback_f(
    func_f: fn(Vec<f64>, f64, Vec<f64>) -> f64,
    params: *const gsl_vector,
    params_len: usize,
    t: f64,
    args: *const gsl_vector,
    args_len: usize
) -> f64 {

    let mut params_vector = Vec::new();
    let mut args_vector = Vec::new();

    for i in 0..params_len {
        unsafe {
            params_vector.push(gsl_vector_get(params, i));
        }
    }

    for i in 0..args_len {
        unsafe {
            args_vector.push(gsl_vector_get(args, i));
        }
    }

    func_f(params_vector, t, args_vector)
}

#[no_mangle]
pub extern "C" fn rust_callback_dfs(
    func_dfs: &Vec<fn(Vec<f64>, f64, Vec<f64>) -> f64>,
    params: *const gsl_vector,
    params_len: usize,
    func_i: usize,
    t: f64,
    args: *const gsl_vector,
    args_len: usize
) -> f64 {

    let mut params_vector = Vec::new();
    let mut args_vector = Vec::new();
    let func_df = func_dfs[func_i];

    for i in 0..params_len {
        unsafe {
            params_vector.push(gsl_vector_get(params, i));
        }
    }

    for i in 0..args_len {
        unsafe {
            args_vector.push(gsl_vector_get(args, i));
        }
    }

    func_df(params_vector, t, args_vector)
}
