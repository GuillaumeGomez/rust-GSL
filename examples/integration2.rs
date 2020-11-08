//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{gamma_beta, IntegrationFixedType, IntegrationFixedWorkspace};
use std::f64::consts::PI;

const M: usize = 10;
const N: usize = 6;

fn main() {
    let t = IntegrationFixedType::hermite();
    let w = IntegrationFixedWorkspace::new(t, N, 0., 1., 0., 0.)
        .expect("IntegrationFixedWorkspace::new failed");

    let (_, result) = w.fixed(|x| x.powi(M as _) + 1.);

    let expected = PI.sqrt() + gamma_beta::gamma::gamma(0.5 * (1. + M as f64));
    println!("m             = {}", M);
    println!("intervals     = {}", w.n());
    println!("result        = {:.18}", result);
    println!("expected      = {:.18}", expected);
    println!("actual error  = {:.18}", result - expected);
}
