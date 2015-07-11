//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following code estimates the derivative of the function f(x) = x^{3/2} at x=2 and at x=0. The function f(x) is undefined for x<0 so
// the derivative at x=0 is computed using gsl_deriv_forward.

extern crate rgsl;

#[allow(unused_variables)]
fn f(x: f64, param: &mut i32) -> f64 {
    x.powf(1.5f64)
}

fn main() {
    let mut result = 0f64;
    let mut abs_err = 0f64;

    println!("f(x) = x^(3/2)");
    rgsl::numerical_differentiation::deriv_central(f, &mut 0i32, 2f64, 1e-8f64, &mut result, &mut abs_err);
    println!("x = 2.0");
    println!("f'(x) = {} +/- {}", result, abs_err);
    println!("exact = {}\n", 1.5 * 2f64.sqrt());

    rgsl::numerical_differentiation::deriv_central(f, &mut 0i32, 0f64, 1e-8f64, &mut result, &mut abs_err);
    println!("x = 0.0");
    println!("f'(x) = {} +/- {}", result, abs_err);
    println!("exact = {}", 0f64);
}