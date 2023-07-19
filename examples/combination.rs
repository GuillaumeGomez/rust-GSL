//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

use rgsl::Combination;

fn main() {
    println!("All subsets of {{0,1,2,3}} by size:");

    for i in 0..=4 {
        let mut c =
            Combination::new_with_init(4, i).expect("Combination::new_init_first failed...");
        loop {
            println!("{:?}", c);
            if c.next().is_err() {
                break;
            }
        }
    }
}
