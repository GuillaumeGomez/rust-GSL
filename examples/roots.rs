//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

fn func_solver() {
    let r_expected = 5.0_f64.sqrt();
    let x_lo = 0.0_f64;
    let x_hi = 5.0_f64;

    let t = rgsl::RootFSolverType::brent();
    let mut solver = rgsl::RootFSolver::new(&t).unwrap();

    let a = 1.0;
    let b = 0.0;
    let c = -5.0;

    solver.set(|x| (a * x + b) + c, x_lo, x_hi);

    println!("using {} method", solver.name());
    println!("iter | [lower  ,  upper] |  root  |   err   | err(est)");

    let mut stop = false;
    for i in 0..100 {
        match solver.iterate() {
            rgsl::Value::Success | rgsl::Value::Continue => {}
            _ => panic!(),
        }

        let r = solver.root();
        let x_lo = solver.x_lower();
        let x_hi = solver.x_upper();
        if let rgsl::Value::Success = rgsl::roots::test_interval(x_lo, x_hi, 0.0, 0.001) {
            println!("Converged:");
            stop = true;
        }
        println!(
            "{:4} | [{x_lo:.4} , {x_hi:.4}] | {r:.4} | {e:+.4} | {err:.4}",
            i + 1,
            x_lo = x_lo,
            x_hi = x_hi,
            r = r,
            e = r - r_expected,
            err = x_hi - x_lo
        );
        if stop {
            break;
        }
    }
}

fn derivative_solver() {
    let r_expected = 5.0_f64.sqrt();
    let mut root = 5.0_f64;

    let t = rgsl::RootFdfSolverType::newton();
    let mut solver = rgsl::RootFdfSolver::new(&t).unwrap();

    let a = 1.0;
    let b = 0.0;
    let c = -5.0;

    solver.set(
        |x| (a * x + b) + c,
        |x| 2.0 * a * x + b,
        |x, y, dy| {
            *y = (a * x + b) + c;
            *dy = 2.0 * a * x + b;
        },
        root,
    );

    println!("using {} method", solver.name());
    println!("iter |  root  |   err   | err(est)");

    let mut stop = false;
    let mut x0;
    for i in 0..100 {
        match solver.iterate() {
            rgsl::Value::Success | rgsl::Value::Continue => {}
            _ => panic!(),
        }

        x0 = root;
        root = solver.root();
        if let rgsl::Value::Success = rgsl::roots::test_delta(root, x0, 0.0, 1e-3) {
            println!("Converged:");
            stop = true;
        }
        println!(
            "{:4} | {r:.4} | {e:+.4} | {err:.4}",
            i + 1,
            r = root,
            e = root - r_expected,
            err = root - x0
        );
        if stop {
            break;
        }
    }
}

fn main() {
    func_solver();
    println!();
    derivative_solver();
}
