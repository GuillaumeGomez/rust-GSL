//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{
    eigen, EigenNonSymmetricVWorkspace, EigenSort, MatrixComplexF64, MatrixF64View,
    VectorComplexF64,
};

fn main() {
    let data = &mut [
        -1., 1., -1., 1., -8., 4., -2., 1., 27., 9., 3., 1., 64., 16., 4., 1.,
    ];
    let mut m = MatrixF64View::from_array(data, 4, 4);
    let mut eval = VectorComplexF64::new(4).expect("VectorF64::new");
    let mut evec = MatrixComplexF64::new(4, 4).expect("MatrixF64::new failed...");
    let mut w =
        EigenNonSymmetricVWorkspace::new(4).expect("EigenNonSymmetricVWorkspace::new failed...");

    m.matrix_mut(|m| {
        w.nonsymmv(m.expect("Failed to get matrix"), &mut eval, &mut evec);
    });

    eigen::nonsymmv_sort(&mut eval, &mut evec, EigenSort::AbsDesc);

    for i in 0..4 {
        let eval_i = eval.get(i);
        evec.column(i, |evec_i| {
            println!("eigenvalue = {} + {}", eval_i.real(), eval_i.imaginary());
            evec_i.expect("Failed to get get column").vector(|v| {
                let v = v.expect("Failed to get vector from column");
                println!("eigenvector = ");
                for j in 0..4 {
                    let z = v.get(j);
                    println!("{} + {}", z.real(), z.imaginary());
                }
            });
        });
    }
}
