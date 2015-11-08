//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The example program below prints all multisets elements containing the values {0,1,2,3} ordered by size. Multiset elements of the same
// size are ordered lexicographically.

extern crate rgsl;
extern crate num;

use std::fs::File;
use std::io::Read;
use rgsl::{wavelet_transforms, sort};
use num::Float;

pub const N : usize = 256;
pub const NC : usize = 20;

#[allow(unused_must_use)]
fn main() {
    let mut data : [f64; 256] = [0f64; 256];
    let mut abscoeff : [f64; 256] = [0f64; 256];
    let mut p : [usize; 256] = [0usize; 256];
    let mut args = Vec::new();

    for entry in std::env::args() {
        args.push(entry);
    }
    let tmp = args[1..].to_vec();

    if tmp.len() < 1 {
        panic!("USAGE: ./wavelet_transforms [file]");
    }

    {
        let mut f = match File::open(&tmp[0]) {
            Ok(f) => f,
            Err(e) => panic!("file error: {}", e),
        };
        for i in 0usize..N {
            // read 8 bytes and parse them as a f64
            let mut b = [0u8;8];
            match f.read(&mut b) {
                Ok(8) => {
                    data[i] = unsafe { ::std::mem::transmute(b) };
                }
                Ok(_) => { /* premature EOF, abort */ break },
                Err(e) => panic!("Read error : {}", e),
            }
        }
    }

    let w = rgsl::Wavelet::new(&rgsl::WaveletType::daubechies(), 4).unwrap();
    let work = rgsl::WaveletWorkspace::new(N).unwrap();

    wavelet_transforms::one_dimension::transform_forward(&w, &mut data, 1, N, &work);

     for i in 0usize..N {
        abscoeff[i] = data[i].abs();
    }

    sort::vectors::sort_index(&mut p, &abscoeff, 1, N);

    let mut i = 0usize;
    while i + NC < N {
        data[p[i] as usize] = 0f64;
        i += 1;
    }

    wavelet_transforms::one_dimension::transform_inverse(&w, &mut data, 1, N, &work);

    for it in 0usize..N {
        println!("{}", data[it]);
    }
}