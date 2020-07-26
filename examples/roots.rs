//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#![allow(non_snake_case)]

extern crate libc;
extern crate rgsl;

use libc::c_void;

#[repr(C)]
struct QuadraticParams {
    a: f64,
    b: f64,
    c: f64,
}

pub extern "C" fn quadratic(x: f64, params: *mut c_void) -> f64 {
    let params = unsafe { &mut *(params as *mut QuadraticParams) };
    (params.a * x + params.b) * x + params.c
}

fn func_solver() {
    let r_expected = 5.0_f64.sqrt();
    let x_lo = 0.0_f64;
    let x_hi = 5.0_f64;

    let t = rgsl::RootFSolverType::brent();
    let mut solver = rgsl::RootFSolver::new(&t).unwrap();
    let mut params = QuadraticParams {
        a: 1.0,
        b: 0.0,
        c: -5.0,
    };
    let mut f = rgsl::RootFunction {
        function: Some(quadratic),
        params: &mut params as *mut _ as *mut c_void,
    };

    solver.set(&mut f, x_lo, x_hi);

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

pub extern "C" fn quadratic_deriv(x: f64, params: *mut c_void) -> f64 {
    let params = unsafe { &mut *(params as *mut QuadraticParams) };
    2.0 * params.a * x + params.b
}

pub extern "C" fn quadratic_fdf(x: f64, params: *mut c_void, y: &mut f64, dy: &mut f64) {
    let params = unsafe { &mut *(params as *mut QuadraticParams) };
    *y = (params.a * x + params.b) * x + params.c;
    *dy = 2.0 * params.a * x + params.b;
}

fn derivative_solver() {
    let r_expected = 5.0_f64.sqrt();
    let mut x = 5.0_f64;
    let mut x0;

    let t = rgsl::RootFdfSolverType::newton();
    let mut solver = rgsl::RootFdfSolver::new(&t).unwrap();
    let mut params = QuadraticParams {
        a: 1.0,
        b: 0.0,
        c: -5.0,
    };
    let mut fdf = rgsl::RootFunctionFdf {
        f: Some(quadratic),
        df: Some(quadratic_deriv),
        fdf: Some(quadratic_fdf),
        params: &mut params as *mut _ as *mut c_void,
    };

    solver.set(&mut fdf, x);

    println!("using {} method", solver.name());
    println!("iter |  root  |   err   | err(est)");

    let mut stop = false;
    for i in 0..100 {
        match solver.iterate() {
            rgsl::Value::Success | rgsl::Value::Continue => {}
            _ => panic!(),
        }

        x0 = x;
        x = solver.root();
        if let rgsl::Value::Success = rgsl::roots::test_delta(x, x0, 0.0, 1e-3) {
            println!("Converged:");
            stop = true;
        }
        println!(
            "{:4} | {r:.4} | {e:+.4} | {err:.4}",
            i + 1,
            r = x,
            e = x - r_expected,
            err = x - x0
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
