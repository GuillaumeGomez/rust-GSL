extern crate pkg_config;

use std::path::Path;


fn main() {

    cc::Build::new()
        .file("src/multifit_nlinear.c")
        .compile("rgslmfnlin");

    if std::process::Command::new("pkg-config").output().is_err() {
        println!("cargo:rustc-link-lib=gsl");
        println!("cargo:rustc-link-lib=gslcblas");
        println!("cargo:rustc-link-lib=rgslmfnlin");
        return;
    }

    if pkg_config::probe_library("gsl").is_err() {
        println!("cargo:rustc-link-lib=gsl");
    }
    if pkg_config::probe_library("gslcblas").is_err() {
        println!("cargo:rustc-link-lib=gslcblas");
    }
    if pkg_config::probe_library("rgslmfnlin").is_err() {
        println!("cargo:rustc-link-lib=rgslmfnlin");
    }

    println!("cargo::rerun-if-changed=src/multifit_nlinear.c");
}
