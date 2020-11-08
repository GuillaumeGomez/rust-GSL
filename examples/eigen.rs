//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{eigen, EigenSort, EigenSymmetricVWorkspace, MatrixF64, MatrixF64View, VectorF64};

fn main() {
    let data = &mut [
        1.,
        1. / 2.,
        1. / 3.,
        1. / 4.,
        1. / 2.,
        1. / 3.,
        1. / 4.,
        1. / 5.,
        1. / 3.,
        1. / 4.,
        1. / 5.,
        1. / 6.,
        1. / 4.,
        1. / 5.,
        1. / 6.,
        1. / 7.,
    ];
    let mut m = MatrixF64View::from_array(data, 4, 4);
    let mut eval = VectorF64::new(4).expect("VectorF64::new");
    let mut evec = MatrixF64::new(4, 4).expect("MatrixF64::new failed...");
    let mut w = EigenSymmetricVWorkspace::new(4).expect("EigenSymmetricVWorkspace::new failed...");

    m.matrix_mut(|m| {
        w.symmv(m.expect("Failed to get matrix"), &mut eval, &mut evec);
    });

    eigen::symmv_sort(&mut eval, &mut evec, EigenSort::AbsAsc);

    for i in 0..4 {
        let eval_i = eval.get(i);
        evec.column(i, |evec_i| {
            println!("eigenvalue = {}", eval_i);
            evec_i.expect("Failed to get get column").vector(|v| {
                let v = v.expect("Failed to get vector from column");
                println!("eigenvector = {:?}", v);
            });
        });
    }
}
