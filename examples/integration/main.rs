//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use std::intrinsics::sinf64;

struct FParams {
    // Amplitude
    a: f64,
    // Phase
    phi: f64
}

fn f(x: f64, p: &mut FParams) -> f64 {
    unsafe { sinf64(p.a * x + p.phi) }
}

fn main() {
    let mut params = FParams {a: 1f64, phi: 0f64};
    let mut result = 0f64;
    let mut error = 0f64;
    let mut n_eval = 0u64;
    
    let xlow = 0f64;
    let xhigh = 10f64;
    let eps_abs = 1e-4f64;
    let eps_rel = 1e-4f64;
  
    match rgsl::integration::qng(f, &mut params, xlow, xhigh, eps_abs, eps_rel, &mut result, &mut error, &mut n_eval) {
        rgsl::enums::Success => {
            println!("Result {} +/- {} from {} evaluations", result, error, n_eval);
        }
        e => {
            println!("There was a problem with integration: {}", e);
        }
    }
}