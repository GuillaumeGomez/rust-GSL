//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::{FftComplexF64WaveTable, FftComplexF64Workspace};

macro_rules! real {
    ($z:ident, $i:expr) => {
        $z[2 * ($i)]
    };
}
macro_rules! imag {
    ($z:ident, $i:expr) => {
        $z[2 * ($i) + 1]
    };
}

const N: usize = 630;

fn main() {
    let data = &mut [0.; 2 * N];
    let wavetable = FftComplexF64WaveTable::new(N).expect("FftComplexF64WaveTable::new failed");
    let mut workspace = FftComplexF64Workspace::new(N).expect("FftComplexF64Workspace::new failed");

    data[0] = 1.;

    for i in 1..=10 {
        real!(data, i) = 1.;
        real!(data, N - i) = 1.;
    }

    for i in 0..N {
        println!("{}: {} {}", i, real!(data, i), imag!(data, i));
    }
    println!("");

    for i in 0..wavetable.nf() {
        println!("# factor {}: {}", i, wavetable.factor()[i]);
    }

    workspace.forward(data, 1, N, &wavetable).unwrap();

    for i in 0..N {
        println!("{}: {} {}", i, real!(data, i), imag!(data, i));
    }
}
