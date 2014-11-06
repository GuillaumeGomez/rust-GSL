//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following program uses the Brent algorithm to find the minimum of the function f(x) = \cos(x) + 1, which occurs at x = \pi. The 
// starting interval is (0,6), with an initial guess for the minimum of 2.

extern crate rgsl;

use std::intrinsics::{cosf64};
use std::f64::consts::PI;

#[allow(unused_variables)]
fn fn1(x: f64, params: &mut f64) -> f64 {
    unsafe { cosf64(x) + 1f64 }
}

#[allow(unused_assignments)]
fn main() {
    let mut iter = 0u;
    let max_iter = 100u;
    let mut m = 2f64;
    let m_expected = PI;
    let mut a = 0f64;
    let mut b = 6f64;
    let mut status = rgsl::Value::Success;

    let t : rgsl::MinimizerType<f64> = rgsl::MinimizerType::brent();
    let mut s : rgsl::Minimizer<f64> = rgsl::Minimizer::new(&t).unwrap();

    s.set(fn1, &mut 0f64, m, a, b);

    println!("using {} method\n", s.name());

    println!("{:5} [{:9}, {:9}] {:9} {:10} {:9}", "iter", "lower", "upper", "min", "err", "err(est)");

    println!("{:5} [{:.7}, {:.7}] {:.7} {:+.7} {:.7}", iter, a, b, m, m - m_expected, b - a);

    loop {
        iter += 1;
        s.iterate();

        m = s.x_minimum();
        a = s.x_lower();
        b = s.x_upper();

        status = rgsl::minimizer::test_interval(a, b, 0.001f64, 0f64);

        if status == rgsl::Value::Success {
            println!("Converged:");
        }

        println!("{:5} [{:.7}, {:.7}] {:.7} {:+.7} {:.7}", iter, a, b, m, m - m_expected, b - a);
        if status != rgsl::Value::Continue || iter >= max_iter {
            break;
        }
    }
}