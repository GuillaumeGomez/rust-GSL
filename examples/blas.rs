//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{blas, CblasTranspose, MatrixF64View};

fn main() {
    let a = &mut [0.11, 0.12, 0.13, 0.21, 0.22, 0.23];
    let b = &mut [1011., 1012., 1021., 1022., 1031., 1032.];
    let c = &mut [0., 0., 0., 0.];

    let mut view_a = MatrixF64View::from_array(a, 2, 3);
    let mut view_b = MatrixF64View::from_array(b, 3, 2);
    let mut view_c = MatrixF64View::from_array(c, 2, 2);

    view_a.matrix_mut(|mat_a| {
        view_b.matrix_mut(|mat_b| {
            view_c.matrix_mut(|mat_c| {
                blas::level3::dgemm(
                    CblasTranspose::NoTranspose,
                    CblasTranspose::NoTranspose,
                    1.,
                    mat_a.expect("Failed to get matrix"),
                    mat_b.expect("Failed to get matrix"),
                    0.,
                    mat_c.expect("Failed to get matrix"),
                );
            });
        });
    });

    println!("[ {}, {}", c[0], c[1]);
    println!("  {}, {} ]", c[2], c[3]);
}
