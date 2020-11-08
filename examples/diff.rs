//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::numerical_differentiation;

fn main() {
    println!("f(x) = x^(3/2)");
    let (_, result, abserr) = numerical_differentiation::deriv_central(|x| x.powf(1.5), 2., 1e-8);
    println!("x = 2.0");
    println!("f'(x) = {:.10} +/- {:.10}", result, abserr);
    println!("exact = {:.10}", 1.5 * 2f64.sqrt());
    println!("");

    let (_, result, abserr) = numerical_differentiation::deriv_forward(|x| x.powf(1.5), 0., 1e-8);
    println!("x = 0.0");
    println!("f'(x) = {:.10} +/- {:.10}", result, abserr);
    println!("exact = {:.10}", 0.0);
}
