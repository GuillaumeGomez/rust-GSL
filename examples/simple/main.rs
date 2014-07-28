/*
 * A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
 */

extern crate rgsl;

fn main() {
	println!("Simple test : {}", rgsl::airy::Airy::Ai(0.5f64, rgsl::types::Gsl::PrecDouble));
}