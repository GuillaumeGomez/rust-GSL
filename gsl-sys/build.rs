extern crate pkg_config;

fn main() {
    if std::process::Command::new("pkg-config").output().is_err() {
        println!("cargo:rustc-link-lib=gsl");
        // println!("cargo:rustc-link-lib=gslcblas");
        return;
    }

    pkg_config::probe_library("gsl").expect("GSL library not found");
    /*pkg_config::probe_library("gslcblas")
    .expect("GSL cblas library not found");*/
}
