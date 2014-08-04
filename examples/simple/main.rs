//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

fn main() {
    println!("Simple Airy test : {}", rgsl::Airy::Ai(0.5f64, rgsl::mode::PrecDouble));
    println!("Simple Bessel test : {}", rgsl::Bessel::I0(0.5f64));
    println!("Simple Canonical test : {}", rgsl::Canonical::half(0.37f64, 1.2f64));
    println!("Simple CBLAS level1 test : {}", rgsl::Cblas::Level1::sdsdot(1i32, 0.6f32, [1.1f32], 1i32, [2.07f32], 1i32));
    println!("Simple Elementary test (acosh(1.0)) : {}", rgsl::Elementary::acosh(1f64));
    println!("Simple Elementary test (asinh(1.0)) : {}", rgsl::Elementary::asinh(1f64));
    println!("Simple Elementary test (atanh(1.0)) : {}", rgsl::Elementary::atanh(1f64));

    println!("=== VectorFloat tests ===");
    let tmp_vec = rgsl::Gsl::VectorFloat::from_slice([1f32, 0f32, 3f32, 2f32]).unwrap();
    println!("min value : {}\nmin value index : {}\n", tmp_vec.min(), tmp_vec.min_index());
    println!("max value : {}\nmax value index : {}\n", tmp_vec.max(), tmp_vec.max_index());
}