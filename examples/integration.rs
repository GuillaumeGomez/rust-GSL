//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use std::f64::consts;

fn main() {
    // This is the example from https://www.gnu.org/software/gsl/doc/html/integration.html#adaptive-integration-example
    let mut w = rgsl::IntegrationWorkspace::new(1000).expect("IntegrationWorkspace::new failed");

    let alpha: f64 = 1.0;
    let expected: f64 = -4.0;

    let (result, error) = w
        .qags(|x| (alpha * x).ln() / x.sqrt(), 0., 1., 0., 1e-7, 1000)
        .unwrap();

    println!("== Adaptive integration ==");
    println!("result          = {}", result);
    println!("exact result    = {}", expected);
    println!("estimated error = {}", error);
    println!("actual error    = {}", result - expected);
    println!("intervals       = {}", w.size());

    // This is the example from https://www.gnu.org/software/gsl/doc/html/integration.html#fixed-point-quadrature-example
    let n = 6;
    let m = 10;
    let w = rgsl::IntegrationFixedWorkspace::new(
        rgsl::IntegrationFixedType::hermite(),
        n,
        0.,
        1.,
        0.,
        0.,
    )
    .expect("IntegrationFixedWorkspace::new failed");

    println!("");

    let result = w.fixed(|x| x.powf(m as _) + 1.).unwrap();

    let expected = if m % 2 == 0 {
        consts::PI.sqrt() + rgsl::gamma_beta::gamma::gamma(0.5 * (1. + m as f64))
    } else {
        consts::PI.sqrt()
    };

    println!("== Fixed-point quadrature ==");
    println!("m               = {}", m);
    println!("intervals       = {}", w.n());
    println!("result          = {}", result);
    println!("exact result    = {}", expected);
    println!("actual error    = {}", result - expected);
}
