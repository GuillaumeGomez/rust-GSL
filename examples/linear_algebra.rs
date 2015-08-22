//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

fn main() {
    let mut a_data : [f64; 16] =
        [0.18, 0.60, 0.57, 0.96,
         0.41, 0.24, 0.99, 0.58,
         0.14, 0.30, 0.97, 0.66,
         0.51, 0.13, 0.19, 0.85];
    let mut b_data : [f64; 4] = [1.0, 2.0, 3.0, 4.0];

    let mut m = rgsl::MatrixView::from_array(&mut a_data, 4, 4);
    let mut b = rgsl::VectorView::from_array(&mut b_data);
    let x = rgsl::VectorF64::new(4).unwrap();
    let mut s = 0i32;
    let p = rgsl::Permutation::new(4).unwrap();

    rgsl::linear_algebra::LU_decomp(&m.matrix(), &p, &mut s);
    rgsl::linear_algebra::LU_solve(&m.matrix(), &p, &b.vector(), &x);

    println!("x = \n{:?}", x);
}