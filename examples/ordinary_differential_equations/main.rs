//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

fn main() {
    let mu = 10.;

    let mut func_eval = 0;
    let mut jac_eval = 0;

    {
        let mut func = |_, y: &[f64], f: &mut [f64]| {
            func_eval += 1;
            f[0] = y[1];
            f[1] = -y[0] - mu * y[1] * (y[0] * y[0] - 1.);
            Ok(())
        };

        let mut jac = |_, y: &[f64], dfdy: &mut [f64], dfdt: &mut [f64]| {
            jac_eval += 1;
            dfdy[0] = 0.;
            dfdy[1] = 1.;
            dfdy[2] = -2. * mu * y[0] * y[1] - 1.;
            dfdy[3] = -mu * (y[0] * y[0] - 1.);

            dfdt[0] = 0.;
            dfdt[1] = 0.;
            Ok(())
        };

        let mut sys = rgsl::ODEiv2System::with_jacobian(2, &mut func, &mut jac);

        let mut d = rgsl::ODEiv2Driver::alloc_y_new(&mut sys, &rgsl::ODEiv2StepType::rk8pd(), 1e-6, 1e-6, 0.).unwrap();

        let mut t = 0.;
        let t1 = 100.;
        let mut y = [1., 0.];

        for i in 1..101 {
            let ti = i as f64 * t1 / 100.;

            match d.apply(&mut t, ti, &mut y) {
                rgsl::Value::Success => {}
                e => {
                    println!("error, return value={:?}", e);
                    break
                }
            }

            println!("{:.5} {:.5} {:.5}", t, y[0], y[1]);
        }
    }

    println!("\nfunc evaluated {} times, jac evaluated {} times", func_eval, jac_eval);
}
