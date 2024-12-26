//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

fn main() {
    let data: [f64; 5] = [17.2, 18.1, 16.5, 18.3, 12.6];

    let mean = rgsl::statistics::mean(&data);
    let variance = rgsl::statistics::variance(&data);
    let largest = rgsl::statistics::max(&data);
    let smallest = rgsl::statistics::min(&data);

    println!(
        "The dataset is {}, {}, {}, {}, {}",
        data[0], data[1], data[2], data[3], data[4]
    );

    println!("The sample mean is {}", mean);
    println!("The estimated variance is {}", variance);
    println!("The largest value is {}", largest);
    println!("The smallest value is {}", smallest);
}
