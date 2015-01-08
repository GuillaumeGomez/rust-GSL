//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;
extern crate libc;
extern crate c_vec;

use libc::c_void;
use c_vec::CVec;

#[allow(unused_variables)]
fn func(t: f64, t_y: *const f64, t_f: *mut f64, params: *mut c_void) -> rgsl::Value {
    unsafe {
        let mu : &mut f64 = ::std::mem::transmute(params);
        let mut f = CVec::new(t_f, 2);
        let y = Vec::from_raw_buf(t_y as *mut f64, 2);
        
        f.as_mut_slice()[0] = y.as_slice()[1];
        f.as_mut_slice()[1] = -y.as_slice()[0] - *mu * y.as_slice()[1] * (y.as_slice()[0] *y.as_slice()[0] - 1f64);
        rgsl::Value::Success
    }
}

#[allow(unused_variables)]
fn jac(t: f64, t_y: *const f64, t_dfdy: *mut f64, t_dfdt: *mut f64, params: *mut c_void) -> rgsl::Value {
    unsafe {
        let mu : &mut f64 = ::std::mem::transmute(params);
        let mut dfdy = CVec::new(t_dfdy, 4);
        let mut dfdy_mat = rgsl::MatrixView::from_array(dfdy.as_mut_slice(), 2, 2);
        let m = dfdy_mat.matrix();
        let mut dfdt = CVec::new(t_dfdt, 2);
        let y = Vec::from_raw_buf(t_y as *mut f64, 2);

        m.set(0, 0, 0f64);
        m.set(0, 1, 1f64);
        m.set(1, 0, -2f64 * *mu * y.as_slice()[0] * y.as_slice()[1] - 1f64);
        m.set(1, 1, -*mu * (y.as_slice()[0] * y.as_slice()[0] - 1f64));
        dfdt.as_mut_slice()[0] = 0f64;
        dfdt.as_mut_slice()[1] = 0f64;
        rgsl::Value::Success
    }
}

fn main() {
    let mut mu = 10f64;
    let sys : rgsl::ODEiv2System = rgsl::ODEiv2System {
        function: func,
        jacobian: jac,
        dimension: 2,
        params: unsafe { ::std::mem::transmute(&mut mu) }
    };

    let d = rgsl::ODEiv2Driver::alloc_y_new(&sys, &rgsl::ODEiv2StepType::rk8pd(), 1e-6f64, 1e-6f64, 0f64).unwrap();
    let mut t = 0f64;
    let t1 = 100f64;
    let mut y : [f64; 2] = [1f64, 0f64];

    for i in range(1u, 101u) {
        let ti = i as f64 * t1 / 100f64;

        match d.apply(&mut t, ti, &mut y) {
            rgsl::Value::Success => {}
            e => {
                println!("error, return value={}", e);
                break
            }
        }

        println!("{:.5} {:.5} {:.5}", t, y[0], y[1]);
    }
}