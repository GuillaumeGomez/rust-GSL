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

    println!("=== integration::qng ===");
    match rgsl::integration::qng(f, &mut params, xlow, xhigh, eps_abs, eps_rel, &mut result, &mut error, &mut n_eval) {
        rgsl::enums::Success => {
            println!("Result {} +/- {} from {} evaluations", result, error, n_eval);
        }
        e => {
            println!("There was a problem with integration: {}", e);
        }
    };

    println!("\n=== IntegrationWorkspace.qag ===");
    let iw = rgsl::IntegrationWorkspace::new(5).unwrap();

    match iw.qag(f, &mut params, xlow, xhigh, eps_abs, eps_rel, 1, rgsl::enums::Gauss15, &mut result, &mut error) {
        rgsl::enums::Success => {
            println!("Result {} +/- {}", result, error);
        }
        e => {
            println!("There was a problem with integration: {}", e);
        }
    };

    println!("\n=== IntegrationWorkspace.qagi ===");
    match  iw.qagi(f, &mut params, 1.0e-7f64, 0f64, iw.limit(), &mut result, &mut error) {
        rgsl::enums::Success => {
            println!("Result {} +/- {}", result, error);
        }
        e => {
            println!("There was a problem with integration: {}", e);
        }
    }
}