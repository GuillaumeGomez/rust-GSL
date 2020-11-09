extern crate pkg_config;

fn main() {
    if std::process::Command::new("pkg-config").output().is_err() {
        println!("cargo:rustc-link-lib=gsl");
        println!("cargo:rustc-link-lib=gslcblas");
        return;
    }

    if pkg_config::probe_library("gsl").is_err() {
        println!("cargo:rustc-link-lib=gsl");
    }
    if pkg_config::probe_library("gslcblas").is_err() {
        println!("cargo:rustc-link-lib=gslcblas");
    }
}
