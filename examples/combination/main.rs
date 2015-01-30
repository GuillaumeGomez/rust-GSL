//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The example program below prints all subsets of the set {0,1,2,3} ordered by size. Subsets of the same size are ordered lexicographically.

#![feature(core)]

extern crate rgsl;

use rgsl::Combination;

fn main() {
    println!("All subsets of [0,1,2,3] by size:");
    for i in range(0, 5) {
        let c = Combination::new_init_first(4, i).unwrap();
        let mut tmp = true;

        println!("size {}:", i);
        while tmp {
            println!("{:?}", c);
            tmp = c.next() == rgsl::Value::Success;
        }
    }
}