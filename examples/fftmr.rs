//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use num_complex::c64;
use rgsl::{FftComplexF64WaveTable, FftComplexF64Workspace};

const N: usize = 128;

fn main() {
    let data = &mut [c64(0., 0.); N];
    let wavetable = FftComplexF64WaveTable::new(N).expect("FftComplexF64WaveTable::new failed");
    let mut workspace = FftComplexF64Workspace::new(N).expect("FftComplexF64Workspace::new failed");

    data[0].re = 1.;

    for i in 1..=10 {
        data[i].re = 1.;
        data[N - i].re = 1.;
    }

    for i in 0..N {
        println!("{}: {}", i, data[i]);
    }
    println!();

    for i in 0..wavetable.nf() {
        println!("# factor {}: {}", i, wavetable.factor()[i]);
    }

    workspace.forward(data, &wavetable).unwrap();

    for i in 0..N {
        println!("{}: {}", i, data[i]);
    }
}
