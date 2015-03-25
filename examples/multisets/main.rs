//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The example program below prints all multisets elements containing the values {0,1,2,3} ordered by size. Multiset elements of the same
// size are ordered lexicographically.

#![feature(core)]
#![feature(old_io)]

extern crate rgsl;

#[allow(unused_must_use)]
fn main() {
    let stdio = &mut ::std::old_io::stdio::stdout().into_inner();

    println!("All multisets of {{0,1,2,3}} by size:");
    for i in 0u64..5u64 {
        let c = rgsl::MultiSet::new_init(4, i).unwrap();
      
        loop {
            write!(stdio, "{{");
            c.print(stdio);
            write!(stdio, " }}\n");
            if c.next() != rgsl::Value::Success {
                break;
            }
        }
    }
}