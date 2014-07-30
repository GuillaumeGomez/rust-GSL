/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

extern crate rgsl;

fn main() {
	println!("Simple Airy test : {}", rgsl::Airy::Ai(0.5f64, rgsl::enums::Gsl::PrecDouble));
	println!("Simple Bessel test : {}", rgsl::Bessel::I0(0.5f64));
	println!("Simple Canonical test : {}", rgsl::Canonical::half(0.37f64, 1.2f64));
	println!("Simple CBLAS level1 test : {}", rgsl::Cblas::Level1::sdsdot(2i32, 0.6f32, 1.1f32, 2i32, 2.07f32, 2i32));
}