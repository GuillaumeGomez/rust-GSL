//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//
use rgsl::integration::{qags, IntegrationWorkspace};
use std::f64::consts;

fn main() -> Result<(), rgsl::Error> {
    // This is the example from https://www.gnu.org/software/gsl/doc/html/integration.html#adaptive-integration-example
    let alpha: f64 = 1.0;
    let expected: f64 = -4.0;

    let mut w = IntegrationWorkspace::new(1000).expect("Workspace");
    let (result, error) = qags(|x| (alpha * x).ln() / x.sqrt(), 0., 1.)
        .epsabs(0.)
        .epsrel(1e-7)
        .workspace(&mut w)
        .val_err()?;

    println!("== Adaptive integration ==");
    println!("result          = {}", result);
    println!("exact result    = {}", expected);
    println!("estimated error = {}", error);
    println!("actual error    = {}", result - expected);
    println!("intervals       = {}", w.size());

    // This is the example from https://www.gnu.org/software/gsl/doc/html/integration.html#fixed-point-quadrature-example
    let n = 6;
    let m = 10;
    let w = rgsl::integration::IntegrationFixedWorkspace::new(
        rgsl::integration::IntegrationFixedType::hermite(),
        n,
        0.,
        1.,
        0.,
        0.,
    )
    .expect("IntegrationFixedWorkspace::new failed");

    println!();

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

    Ok(())
}
