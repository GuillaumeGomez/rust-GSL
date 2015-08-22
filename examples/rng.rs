//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/// Running the program without any environment variables uses the initial defaults, an mt19937 generator with a seed of 0,
/// 
/// ```Shell
/// > ./a.out
/// generator type: mt19937
/// seed = 0
/// first value = 4293858116
/// ```
/// By setting the two variables on the command line we can change the default generator and the seed,
/// 
/// ```Shell
/// > GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./a.out 
/// GSL_RNG_TYPE=taus
/// GSL_RNG_SEED=123
/// generator type: taus
/// seed = 123
/// first value = 2720986350
/// ```

extern crate rgsl;

use rgsl::{RngType};

fn main() {
    rgsl::RngType::env_setup();
    let t : RngType = rgsl::rng::default();
    let r = rgsl::Rng::new(&t).unwrap();

    println!("=== DEFAULT ===");
    println!("generator type: {}", r.get_name());
    println!("seed = {}", rgsl::Rng::default_seed());
    println!("first value = {}", r.get());

    println!("\n=== Available generators ===");
    let v = rgsl::RngType::types_setup();

    for tmp in v.iter() {
        println!("generator type: {}", tmp.name());
    }
}