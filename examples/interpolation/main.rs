//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

fn main() {
    let mut x : [f64; 10] = [0f64; 10];
    let mut y : [f64; 10] = [0f64; 10];

    println!("#m=0,S=2");

    for i in 0usize..10usize {
        x[i] = i as f64 + 0.5f64 * (i as f64).sin();
        y[i] = i as f64 + ((i * i) as f64).cos();
        println!("{} {}", x[i], y[i]);
    }

    println!("#m=1,S=0");

    {
        let mut acc = rgsl::InterpAccel::new();
        let spline = rgsl::Spline::new(&rgsl::InterpType::cspline(), 10).unwrap();

        spline.init(&x, &y);

        let mut xi = x[0];
        while xi < x[9] {
            let yi = spline.eval(xi, &mut acc);

            println!("{} {}", xi, yi);
            xi += 0.01f64;
        }
    }
}