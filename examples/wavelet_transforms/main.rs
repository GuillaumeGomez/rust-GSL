//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The example program below prints all multisets elements containing the values {0,1,2,3} ordered by size. Multiset elements of the same
// size are ordered lexicographically.

#![feature(core)]
#![feature(io)]
#![feature(os)]
#![feature(collections)]
#![feature(path)]

extern crate rgsl;

use std::old_io::{File, Open, Read};
use rgsl::{wavelet_transforms, sort};
use std::os;
use std::num::Float;

pub const N : usize = 256;
pub const NC : usize = 20;

#[allow(unused_must_use)]
fn main() {
    let mut data : [f64; 256] = [0f64; 256];
    let mut abscoeff : [f64; 256] = [0f64; 256];
    let mut p : [u64; 256] = [0u64; 256];
    let args = os::args();
    let tmp = args.tail();

    if tmp.len() < 1 {
        panic!("USAGE: ./wavelet_transforms [file]");
    }

    {
        let p = Path::new(tmp[0].as_slice());
        let mut f = match File::open_mode(&p, Open, Read) {
            Ok(f) => f,
            Err(e) => panic!("file error: {}", e),
        };
        for i in range(0us, N) {
            match f.read_be_f64() {
                Ok(v) => {
                    data[i] = v;
                }
                Err(e) => panic!("Read error : {}", e),
            }
        }
    }

    let w = rgsl::Wavelet::new(&rgsl::WaveletType::daubechies(), 4u64).unwrap();
    let work = rgsl::WaveletWorkspace::new(N as u64).unwrap();

    wavelet_transforms::one_dimension::transform_forward(&w, &mut data, 1, N as u64, &work);

     for i in range(0us, N) {
        abscoeff[i] = data[i].abs();
    }

    sort::vectors::sort_index(&mut p, &abscoeff, 1, N as u64);

    let mut i = 0us;
    while i + NC < N {
        data[p[i] as usize] = 0f64;
        i += 1;
    }

    wavelet_transforms::one_dimension::transform_inverse(&w, &mut data, 1, N as u64, &work);

    for it in range(0us, N) {
        println!("{}", data[it]);
    }
}