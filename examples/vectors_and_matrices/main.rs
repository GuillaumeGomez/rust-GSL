//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// This program shows how to use matrices and vectors.

#![allow(unstable)]

extern crate rgsl;

fn main() {
    println!("=== VECTOR PART ===");
    let v = rgsl::VectorF64::new(3).unwrap();
  
    for i in range(0u64, 3u64) {
        v.set(i, 1.23f64 + i as f64);
    }

    for i in range(0u64, 3u64) {
        println!("v_{} = {}", i, v.get(i));
    }

    // or you can do :
    // println!("{}", c);

    println!("=== MATRIX PART ===");
    let m = rgsl::MatrixF64::new(10, 3).unwrap();
  
    for i in range(0, 10) {
        for j in range(0, 3) {
            m.set(i, j, 0.23f64 + 100f64 * i as f64 + j as f64);
        }
    }

    for i in range(0, 10) {
        for j in range(0, 3) {
            println!("m({},{}) = {}", i, j, m.get(i, j));
        }
    }
    // or you can do :
    // println!("{}", m);
}