//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use std::default::Default;
use rgsl::{Pow/*, Elementary*/, Trigonometric};

fn main() {
    println!("=== VectorFloat tests ===");
    let mut tmp_vec = rgsl::VectorF32::from_slice(&[1f32, 0f32, 3f32, 2f32]).unwrap();
    let mut tmp_vec2 = rgsl::VectorF32::from_slice(&[14f32, 6f32, -3f32, 1.2f32]).unwrap();
    println!("min value : {}\nmin value index : {}", tmp_vec.min(), tmp_vec.min_index());
    println!("max value : {}\nmax value index : {}", tmp_vec.max(), tmp_vec.max_index());
    println!("{:?}", tmp_vec);
    println!("sswap : {:?}", rgsl::blas::level1::sswap(&mut tmp_vec, &mut tmp_vec2));

    println!("\n=== MatrixFloat tests ===");
    let mut tmp_mat = rgsl::MatrixF32::new(2, 3).unwrap();
    tmp_mat.set(1, 2, 42f32);
    tmp_mat.set(0, 0, 1f32);
    match tmp_mat.min_index() {
        (i, j) => {println!("min value : {}\nmin value index : {}-{}", tmp_mat.min(), i, j);}
    };
    match tmp_mat.max_index() {
        (i, j) => {println!("max value : {}\nmax value index : {}-{}", tmp_mat.max(), i, j);}
    };
    println!("{:?}", tmp_mat);

    println!("\n=== Complex tests ===");
    let mut tmp_complex : rgsl::ComplexF64 = Default::default();
    tmp_complex.data[0] = -1f64;
    println!("abs : {}", tmp_complex.abs());
    println!("add_real : {:?}", tmp_complex.add_real(14f64));
    tmp_complex.data[0] = 3f64;
    println!("sqrt : {:?}", tmp_complex.sqrt());
    tmp_complex.data[1] = 14f64;
    println!("sqrt : {:?}", tmp_complex.sqrt());

    println!("\n=== Modules tests ===");
    println!("Simple Airy::Ai test : {}", rgsl::airy::Ai(0.5f64, rgsl::Mode::PrecDouble));
    println!("Simple Bessel::I0 test : {}", rgsl::bessel::I0(0.5f64));
    println!("Simple Legendre::conical test : {}", rgsl::legendre::conical::half(0.37f64, 1.2f64));
    println!("Simple Legendre::conical test : {}", rgsl::legendre::conical::half(0.37f64, 1.2f64));
    println!("Simple BLAS::level1::snrm2 test : {}", rgsl::blas::level1::snrm2(&tmp_vec));
    println!("Simple Logarithm::log test : {}", rgsl::logarithm::log(0.1f64));
    tmp_mat = rgsl::MatrixF32::new(4, 4).unwrap();
    tmp_mat.set(1, 2, 42f32);
    tmp_mat.set(0, 0, 1f32);
    tmp_mat.set(1, 2, 3f32);
    tmp_mat.set(1, 3, 0.5f32);
    println!("\n=> Simple BLAS level2 test before :\n{:?}", tmp_mat);
    rgsl::blas::level2::sger(1.7f32, &tmp_vec, &rgsl::VectorF32::from_slice(&[0.4f32, 14f32, 3f32, 2f32]).unwrap(), &mut tmp_mat);
    println!("=> Simple BLAS level2 test after :\n{:?}", tmp_mat);
    println!("\nSimple CBLAS level1 test : {}", rgsl::cblas::level1::sdsdot(1i32, 0.6f32, &[1.1f32], 1i32, &[2.07f32], 1i32));
    println!("Simple Elementary test (acosh(1.0)) : {}", 1f64.acosh());
    println!("Simple Elementary test (asinh(1.0)) : {}", 1f64.asinh());
    println!("Simple Elementary test (atanh(1.0)) : {}", 1f64.atanh());
    println!("Simple Trigonometric test sin(1.0) : {}", 1f64.sin());
    println!("Simple Trigonometric test cos(1.0) : {}", 1f64.cos());
    println!("Simple Trigonometric test sf_hypot(1.0) : {}", 1f64.sf_hypot());

    println!("\n=== Fit tests ===");
    let x = [1970f64, 1980f64, 1990f64, 2000f64];
    let y = [12f64, 11f64, 14f64, 13f64];
    let w = [0.1f64, 0.2f64, 0.3f64, 0.4f64];
    let mut c0 = 0f64;
    let mut c1 = 0f64;
    let mut cov00 = 0f64;
    let mut cov01 = 0f64;
    let mut cov11 = 0f64;
    let mut chisq = 0f64;

    rgsl::fit::wlinear(&x, 1, &w, 1, &y, 1, x.len(), &mut c0, &mut c1, &mut cov00, &mut cov01, &mut cov11, &mut chisq);
    println!("=> wlinear test :");
    println!("best fit: Y = {} + {} X", c0, c1);
    println!("covariance matrix:");
    println!("[{}, {}\n {}, {}]", cov00, cov01, cov01, cov11);
    println!("chisq = {}\n", chisq);

    let dx = [2.3f64, 2.4f64, 2.1f64, 2.6f64, 3f64];
    let dfx = [10f64, 12f64, 15f64, 8f64, 16f64];
    let mut sumsq = 0f64;

    rgsl::fit::mul(&dx, 1, &dfx, 1, dx.len(), &mut c1, &mut cov11, &mut sumsq);
    println!("=> mul test :");
    for i in 0..dx.len() {
        println!("dfx[{}]/dx[{}] = {} / {} = {}", i, i, dfx[i], dx[i], dfx[i] / dx[i]);
    }
    println!("---------------------------");
    println!("Y = c1 X ==> c1 = {}", c1);
    println!("cov11 = {}", cov11);
    println!("sumsq = {}", sumsq);

    println!("\n=== Pow tests ===");
    println!("pow_int(2, 3) : {}", 2f64.pow_int(3));
    println!("pow3(2) : {}", 2f64.pow3());
    println!("pow9(2) : {}", 2f64.pow9());
}