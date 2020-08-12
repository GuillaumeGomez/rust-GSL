//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following example programs demonstrate the use of ntuples in managing a large dataset. The first part creates a set of 10,000
// simulated “events”, each with 3 associated values (x,y,z). These are generated from a Gaussian distribution with unit variance, for
// demonstration purposes, and written to the ntuple file test.dat
//
// The second part analyses the ntuple data in the file test.dat. The analysis procedure is to compute the squared-magnitude of each
// event, E^2=x^2+y^2+z^2, and select only those which exceed a lower limit of 1.5. The selected events are then histogrammed using their E^2 values.

extern crate rgsl;

struct Data {
    x: f64,
    y: f64,
    z: f64,
}

fn sel_func(data: &mut Data, scale: &mut f64) -> bool {
    let x = data.x;
    let y = data.y;
    let z = data.z;

    let e2 = x * x + y * y + z * z;

    e2 > *scale
}

#[allow(unused_variables)]
fn val_func(data: &mut Data, params: &mut i32) -> f64 {
    let x = data.x;
    let y = data.y;
    let z = data.z;

    x * x + y * y + z * z
}

fn first_part(r: &mut rgsl::Rng) {
    let mut ntuple_row = Data {
        x: 0f64,
        y: 0f64,
        z: 0f64,
    };
    let ntuple = rgsl::NTuples::create("test.dat", &mut ntuple_row).unwrap();

    for _ in 0usize..10000usize {
        ntuple_row.x = rgsl::randist::gaussian::ugaussian(r);
        ntuple_row.y = rgsl::randist::gaussian::ugaussian(r);
        ntuple_row.z = rgsl::randist::gaussian::ugaussian(r);

        ntuple.write();
    }
}

fn second_part() {
    let mut ntuple_row = Data {
        x: 0f64,
        y: 0f64,
        z: 0f64,
    };

    let ntuple = rgsl::NTuples::open("test.dat", &mut ntuple_row).unwrap();
    let mut lower = 1.5f64;

    let mut h = rgsl::Histogram::new(100).unwrap();
    h.set_ranges_uniform(0f64, 10f64);

    ntuple.project(&mut h, val_func, &mut 0i32, sel_func, &mut lower);
    //gsl_histogram_fprintf(stdout, h, "%f", "%f");
    h.print(&mut ::std::io::stdout()).expect("Failed to print");
}

fn main() {
    rgsl::RngType::env_setup();
    let t: rgsl::RngType = rgsl::rng::default();
    let mut r = rgsl::Rng::new(&t).unwrap();

    first_part(&mut r);
    second_part();
}
