extern crate pkg_config;

fn main() {
    if std::process::Command::new("pkg-config").output().is_err() {
        println!("cargo:rustc-link-lib=gsl");
        println!("cargo:rustc-link-lib=gslcblas");
        return;
    }

    pkg_config::probe_library("gsl")
        .map(|lib| {
            if lib.version.starts_with("2.") {
                println!(r#"cargo:rustc-cfg=feature="v2""#);
            }
        })
        .expect("GSL library not found");
}
