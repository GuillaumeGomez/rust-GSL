fn main() {
    if std::process::Command::new("pkg-config").output().is_err() {
        println!("cargo:rustc-link-lib=gsl");
        println!("cargo:rustc-link-lib=gslcblas");
        return;
    }

    // pkg_config::probe_library("gsl")
    //     // .map(|lib| {
    //     //     if lib.version.starts_with("2.") && ::std::env::var_os("IGNORE_VERSION").is_none() {
    //     //         println!(r#"cargo:rustc-cfg=feature="v2""#);
    //     //     }
    //     // })
    //     .expect("GSL library not found");
    // pkg_config::probe_library("gslcblas")
    //     .expect("GSL cblas library not found");
}
