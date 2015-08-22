//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The example program below uses the Monte Carlo routines to estimate the value of the following 3-dimensional integral from the theory of
// random walks,
// 
// I = \int_{-pi}^{+pi} {dk_x/(2 pi)} 
//     \int_{-pi}^{+pi} {dk_y/(2 pi)} 
//     \int_{-pi}^{+pi} {dk_z/(2 pi)} 
//      1 / (1 - cos(k_x)cos(k_y)cos(k_z)).
// The analytic value of this integral can be shown to be I = \Gamma(1/4)^4/(4 \pi^3) = 1.393203929685676859.... The integral gives the mean
// time spent at the origin by a random walk on a body-centered cubic lattice in three dimensions.
// 
// For simplicity we will compute the integral over the region (0,0,0) to (\pi,\pi,\pi) and multiply by 8 to obtain the full result. The
// integral is slowly varying in the middle of the region but has integrable singularities at the corners (0,0,0), (0,\pi,\pi), (\pi,0,\pi)
// and (\pi,\pi,0). The Monte Carlo routines only select points which are strictly within the integration region and so no special measures
// are needed to avoid these singularities.

extern crate rgsl;
extern crate num;

use std::f64::consts::PI;
use num::Float;

/* Computation of the integral,

      I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))

   over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
   is Gamma(1/4)^4/(4 pi^3).  This example is taken from
   C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
   Volume 1", Section 1.1, p21, which cites the original
   paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
   1800 (1977) */

/* For simplicity we compute the integral over the region 
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */

const EXACT : f64 = 1.3932039296856768591842462603255f64;

#[allow(unused_variables)]
fn g(k: &mut [f64], params: &mut f64) -> f64 {
    let a = 1f64 / (PI * PI * PI);
    
    a / (1.0 - k[0].cos() * k[1].cos() * k[2].cos())
}

fn display_results(title: &str, result: f64, error: f64) {
    println!("{} ==================", title);
    println!("result = {:.6}", result);
    println!("sigma  = {:.6}", error);
    println!("exact  = {:.6}", EXACT);
    println!("error  = {:.6} = {:.2} sigma", result - EXACT, (result - EXACT).abs() / error);
}

#[allow(unused_assignments)]
fn main() {
    let mut res = 0f64;
    let mut err = 0f64;

    let xl : [f64; 3] = [0f64; 3];
    let xu : [f64; 3] = [PI, PI, PI];

    let calls = 500000;

    rgsl::RngType::env_setup();
    let t : rgsl::RngType = rgsl::rng::default();
    let r = rgsl::Rng::new(&t).unwrap();

    {
        let s = rgsl::PlainMonteCarlo::new(3).unwrap();
        
        s.integrate(g, &mut 0f64, &xl, &xu, calls, &r, &mut res, &mut err);
        display_results("plain", res, err);
    }

    {
        let s = rgsl::MiserMonteCarlo::new(3).unwrap();
        
        s.integrate(g, &mut 0f64, &xl, &xu, calls, &r, &mut res, &mut err);
        display_results("miser", res, err);
    }

    {
        let s = rgsl::VegasMonteCarlo::new(3).unwrap();

        s.integrate(g, &mut 0f64, &xl, &xu, 10000, &r, &mut res, &mut err);
        display_results("vegas warm-up", res, err);

        println!("converging...");

        loop {
            s.integrate(g, &mut 0f64, &xl, &xu, calls / 5, &r, &mut res, &mut err);
            println!("result = {:.6} sigma = {:.6} chisq/dof = {:.1}", res, err, s.chisq());
            if (s.chisq() - 1f64).abs() <= 0.5f64 {
                break;
            }
        }

        display_results("vegas final", res, err);
    }
}