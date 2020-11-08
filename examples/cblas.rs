//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{cblas, CblasOrder, CblasTranspose};

fn main() {
    let lda = 3;
    let a = &[0.11, 0.12, 0.13, 0.21, 0.22, 0.23];

    let ldb = 2;
    let b = &[1011., 1012., 1021., 1022., 1031., 1032.];

    let ldc = 2;
    let c = &mut [0., 0., 0., 0.];

    // Compute C = A B
    cblas::level3::sgemm(
        CblasOrder::RowMajor,
        CblasTranspose::NoTranspose,
        CblasTranspose::NoTranspose,
        2,
        2,
        3,
        1.,
        a,
        lda,
        b,
        ldb,
        0.,
        c,
        ldc,
    );

    println!("[ {}, {}", c[0], c[1]);
    println!("  {}, {} ]", c[2], c[3]);
}
