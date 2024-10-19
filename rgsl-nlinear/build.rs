extern crate pkg_config;

fn main() {
    cc::Build::new()
        .file("src/multifit_nlinear.c")
        .include("./rgsl-cblib/target/debug")
        .flag("-lrgsl_cblib")
        .compile("rgslmfnlin");

    println!("cargo:rustc-link-search=./target/debug");
    println!("cargo:rustc-link-search=./rgsl-cblib/target/debug");

    if std::process::Command::new("pkg-config").output().is_err() {
        println!("cargo:rustc-link-lib=gsl");
        println!("cargo:rustc-link-lib=gslcblas");
        println!("cargo:rustc-link-lib=rgsl_cblib");
        println!("cargo:rustc-link-lib=rgslmfnlin");
        return;
    }

    if pkg_config::probe_library("gsl").is_err() {
        println!("cargo:rustc-link-lib=gsl");
    }
    if pkg_config::probe_library("gslcblas").is_err() {
        println!("cargo:rustc-link-lib=gslcblas");
    }
    if pkg_config::probe_library("rgsl_cblib").is_err() {
        println!("cargo:rustc-link-lib=rgsl_cblib");
    }
    if pkg_config::probe_library("rgslmfnlin").is_err() {
        println!("cargo:rustc-link-lib=rgslmfnlin");
    }

    println!("cargo::rerun-if-changed=src/multifit_nlinear.c");
}
