//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use std::default::Default;

fn main() {
    println!("=== VectorFloat tests ===");
    let tmp_vec = rgsl::Gsl::VectorFloat::from_slice([1f32, 0f32, 3f32, 2f32]).unwrap();
    println!("min value : {}\nmin value index : {}", tmp_vec.min(), tmp_vec.min_index());
    println!("max value : {}\nmax value index : {}", tmp_vec.max(), tmp_vec.max_index());
    println!("{}", tmp_vec);

    println!("\n=== MatrixFloat tests ===");
    let mut tmp_mat = rgsl::Gsl::MatrixFloat::new(2u64, 3u64).unwrap();
    tmp_mat.set(1, 2, 42f32);
    tmp_mat.set(0, 0, 1f32);
    match tmp_mat.min_index() {
        (i, j) => {println!("min value : {}\nmin value index : {}-{}", tmp_mat.min(), i, j);}
    };
    match tmp_mat.max_index() {
        (i, j) => {println!("max value : {}\nmax value index : {}-{}", tmp_mat.max(), i, j);}
    };
    println!("{}", tmp_mat);

    println!("\n=== Complex tests ===");
    let mut tmp_complex : rgsl::Gsl::Complex = Default::default();
    tmp_complex.data[0] = -1f64;
    println!("abs : {}", tmp_complex.abs());
    println!("add_real : {}", tmp_complex.add_real(14f64));
    tmp_complex.data[0] = 3f64;
    println!("sqrt : {}", tmp_complex.sqrt());
    tmp_complex.data[1] = 14f64;
    println!("sqrt : {}", tmp_complex.sqrt());

    println!("\n=== Modules tests ===");
    println!("Simple Airy test : {}", rgsl::Airy::Ai(0.5f64, rgsl::mode::PrecDouble));
    println!("Simple Bessel test : {}", rgsl::Bessel::I0(0.5f64));
    println!("Simple Canonical test : {}", rgsl::Canonical::half(0.37f64, 1.2f64));
    println!("Simple BLAS level1 test : {}", rgsl::Blas::Level1::snrm2(&tmp_vec));
    tmp_mat = rgsl::Gsl::MatrixFloat::new(4u64, 4u64).unwrap();
    tmp_mat.set(1, 2, 42f32);
    tmp_mat.set(0, 0, 1f32);
    tmp_mat.set(1, 2, 3f32);
    tmp_mat.set(1, 3, 0.5f32);
    println!("\n=> Simple BLAS level2 test before :\n{}", tmp_mat);
    rgsl::Blas::Level2::sger(1.7f32, &tmp_vec, &rgsl::Gsl::VectorFloat::from_slice([0.4f32, 14f32, 3f32, 2f32]).unwrap(), &mut tmp_mat);
    println!("=> Simple BLAS level2 test after :\n{}", tmp_mat);
    println!("\nSimple CBLAS level1 test : {}", rgsl::Cblas::Level1::sdsdot(1i32, 0.6f32, [1.1f32], 1i32, [2.07f32], 1i32));
    println!("Simple Elementary test (acosh(1.0)) : {}", rgsl::Elementary::acosh(1f64));
    println!("Simple Elementary test (asinh(1.0)) : {}", rgsl::Elementary::asinh(1f64));
    println!("Simple Elementary test (atanh(1.0)) : {}", rgsl::Elementary::atanh(1f64));
}