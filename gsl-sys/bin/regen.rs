extern crate regex;

use std::collections::HashSet;
use std::env;
use std::ffi::OsStr;
use std::fs::{read_dir, write, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::process::Command;

use regex::Regex;

const HEADER_FILE: &str = "wrapper.h";
const GSL_REPOSITORY: &str = "git://git.savannah.gnu.org/gsl.git";
const STRUCTS_TO_IGNORE: &[&str] = &["pub struct _IO_", "pub struct _G_"];
const TYPES_TO_IGNORE: &[&str] = &["pub type _IO_lock_t "];
const STATICS_TO_IGNORE: &[&str] = &[
    "pub static mut stdin:",
    "pub static mut stderr:",
    "pub static mut stdout:",
    "pub static mut sys_nerr:",
    "pub static mut sys_errlist:",
];
const VERSIONS: &[(&str, &str)] = &[
    // 2.1
    ("pub fn gsl_multifit_linear_rcond(", "v2_1"),
    // --> little exception here...
    ("pub fn gsl_multilarge_linear_matrix_ptr(", "v2_7"),
    // --> little exception here...
    ("pub fn gsl_multilarge_linear_rhs_ptr(", "v2_7"),
    ("pub struct gsl_multilarge_*;", "v2_1"),
    // --> little exception here...
    ("pub fn gsl_multilarge_linear_lcurve(", "v2_2"),
    ("pub fn gsl_multilarge_*(", "v2_1"),
    ("pub static mut gsl_multilarge_*:", "v2_1"),
    // 2.2
    ("pub fn gsl_linalg_tri_*(", "v2_2"),
    ("pub fn gsl_linalg_COD_*(", "v2_2"),
    ("pub fn gsl_rstat_quantile_reset(", "v2_2"),
    ("pub fn gsl_ran_multivariate_gaussian*(", "v2_2"),
    ("pub fn gsl_linalg_cholesky_rcond(", "v2_2"),
    ("pub fn gsl_linalg_QRPT_rcond(", "v2_2"),
    ("pub fn gsl_linalg_QRPT_lssolve(", "v2_2"),
    ("pub fn gsl_linalg_QRPT_lssolve2(", "v2_2"),
    ("pub fn gsl_permute_matrix(", "v2_2"),
    ("pub fn gsl_linalg_mcholesky_*(", "v2_2"),
    ("pub fn gsl_linalg_pcholesky_*(", "v2_2"),
    ("pub fn gsl_rstat_rms(", "v2_2"),
    // 2.3
    ("pub fn gsl_multifit_linear_tsvd(", "v2_3"),
    ("pub fn gsl_multifit_wlinear_tsvd(", "v2_3"),
    ("pub fn gsl_multifit_linear_rank(", "v2_3"),
    // 2.5
    ("pub fn gsl_ran_wishart*(", "v2_5"),
    ("pub const gsl_filter_*:", "v2_5"),
    ("pub type gsl_filter_* = ", "v2_5"),
    ("pub struct gsl_filter_*;", "v2_5"),
    ("pub fn gsl_filter_*(", "v2_5"),
    ("pub const gsl_movstat_*:", "v2_5"),
    ("pub static mut gsl_movstat_*:", "v2_5"),
    ("pub type gsl_movstat_* = ", "v2_5"),
    ("pub struct gsl_movstat_*;", "v2_5"),
    ("pub fn gsl_movstat_*(", "v2_5"),
    ("pub fn gsl_stats_median(", "v2_5"),
    ("pub fn gsl_stats_select(", "v2_5"),
    ("pub fn gsl_stats_mad(", "v2_5"),
    ("pub fn gsl_stats_mad0(", "v2_5"),
    ("pub fn gsl_stats_Sn_from_sorted_data(", "v2_5"),
    ("pub fn gsl_stats_Qn_from_sorted_data(", "v2_5"),
    ("pub fn gsl_stats_gastwirth_from_sorted_data(", "v2_5"),
    ("pub fn gsl_stats_trmean_from_sorted_data(", "v2_5"),
    ("pub struct gsl_integration_romberg_workspace;", "v2_5"),
    ("pub fn gsl_integration_romberg*(", "v2_5"),
    // 2.6
    ("pub fn gsl_vector_axpby(", "v2_6"),
    ("pub fn gsl_linalg_ldlt_*(", "v2_6"),
    ("pub static mut gsl_bst_*: ", "v2_6"),
    ("pub type gsl_bst_* = ", "v2_6"),
    ("pub struct gsl_bst_*;", "v2_6"),
    ("pub fn gsl_bst_*(", "v2_6"),
    ("pub fn gsl_spmatrix_scale_columns(", "v2_6"),
    ("pub fn gsl_spmatrix_scale_rows(", "v2_6"),
    ("pub fn gsl_spmatrix_add_to_dense(", "v2_6"),
    ("pub fn gsl_spmatrix_min_index(", "v2_6"),
    // --> little exception here...
    ("pub fn gsl_linalg_cholesky_band_solvem(", "v2_7"),
    // --> little exception here...
    ("pub fn gsl_linalg_cholesky_band_svxm(", "v2_7"),
    ("pub fn gsl_linalg_cholesky_band_*(", "v2_6"),
    ("pub fn gsl_linalg_LQ_lssolve(", "v2_6"),
    // 2.7
    ("pub fn gsl_linalg_LU_band_*(", "v2_7"),
    ("pub fn gsl_matrix_norm1(", "v2_7"),
    ("pub fn gsl_spmatrix_norm1(", "v2_7"),
    ("pub fn gsl_matrix_complex_conjtrans_memcpy(", "v2_7"),
    ("pub fn gsl_linalg_QL_*(", "v2_7"),
    ("pub fn gsl_linalg_complex_QR_*(", "v2_7"),
    ("pub fn gsl_vector_sum(", "v2_7"),
    ("pub fn gsl_matrix_scale_rows(", "v2_7"),
    ("pub fn gsl_matrix_scale_columns(", "v2_7"),
    ("pub fn gsl_spmatrix_dense_sub(", "v2_7"),
];

fn get_all_headers(folder: &Path, extra: &mut Vec<String>, headers: &mut Vec<String>) {
    println!("=> Entering `{:?}`", folder);
    for entry in read_dir(folder).expect("Failed to read gsl directory...") {
        let entry = entry.expect("failed to get entry...");
        let entry = entry.path();
        if entry.is_dir() {
            extra.push(
                entry
                    .file_name()
                    .expect("failed to get file name")
                    .to_str()
                    .expect("failed to convert to str")
                    .to_owned(),
            );
            get_all_headers(&entry, extra, headers);
            extra.pop();
        } else if entry.is_file() && *entry.extension().as_ref().unwrap_or(&OsStr::new("")) == "h" {
            let file_name = entry
                .file_name()
                .expect("failed to get file name")
                .to_str()
                .expect("failed to convert to str");
            headers.push(format!("#include \"{}/{}\"", extra.join("/"), file_name));
            println!("--> Added `{}` to the list!", file_name);
        }
    }
    println!("<= Leaving `{:?}`", folder);
}

fn create_header_file(folder: &Path) {
    println!("=> Creating header file...");
    let mut headers = Vec::new();
    let mut extra = vec!["gsl".to_owned()];

    get_all_headers(&folder.join("gsl"), &mut extra, &mut headers);
    write(
        HEADER_FILE,
        format!(
            "#ifndef __WRAPPER__\n#define __WRAPPER__\n{}\n#endif\n",
            headers.join("\n")
        ),
    )
    .expect("failed to write content to wrapper header file...");
    println!("<= Done");
}

fn first_pass(mut content: Vec<&str>) -> Vec<String> {
    println!("=> Running first pass...");
    let mut consts = HashSet::new();
    let mut pos = 0;

    while pos < content.len() {
        if content[pos].starts_with("pub const _") {
            content.remove(pos);
            continue;
        } else if content[pos].starts_with("pub const ") {
            if !consts.insert(content[pos].split(":").next().unwrap()) {
                content.remove(pos);
                continue;
            }
        } else if STRUCTS_TO_IGNORE
            .iter()
            .any(|s| content[pos].starts_with(s))
        {
            while pos > 1 && content[pos - 1].starts_with("#[") {
                pos -= 1;
            }
            while pos < content.len() && !content[pos].starts_with("}") {
                content.remove(pos);
            }
            if pos < content.len() {
                content.remove(pos);
            }
            continue;
        } else if content[pos].starts_with("pub type FILE = ") {
            content[pos] = "pub type FILE = libc::FILE;";
        } else if content[pos].starts_with("pub struct ") {
            pos += 1;
            while pos < content.len() && !content[pos].starts_with("}") {
                content.remove(pos);
            }
            if pos < content.len() {
                content.remove(pos);
            }
        } else if content[pos].starts_with("pub union ") {
            while pos > 0 && content[pos - 1].starts_with("#[") {
                content.remove(pos - 1);
                pos -= 1;
            }
            content.remove(pos);
            while pos < content.len() && !content[pos].starts_with("}") {
                content.remove(pos);
            }
            if pos < content.len() {
                content.remove(pos);
            }
        } else if TYPES_TO_IGNORE.iter().any(|s| content[pos].starts_with(s)) {
            content.remove(pos);
            continue;
        } else {
            let should_remove = if let Some(fn_name) = content[pos]
                .trim_start()
                .split("(")
                .next()
                .unwrap()
                .split("pub fn ")
                .skip(1)
                .next()
            {
                !fn_name.starts_with("gsl_") && !fn_name.starts_with("cblas_")
            } else {
                false
            };
            if should_remove {
                while !content[pos].starts_with("extern \"C\" {") {
                    if pos > 0 {
                        pos -= 1;
                    } else {
                        break;
                    }
                }
                while !content[pos].starts_with("}") && pos < content.len() {
                    content.remove(pos);
                }
                if pos < content.len() {
                    content.remove(pos);
                }
                continue;
            }
        }
        pos += 1;
    }
    content
        .into_iter()
        .map(|x| x.to_owned())
        .collect::<Vec<_>>()
}

fn add_features(content: &mut Vec<String>) {
    println!("=> Running add_features...");
    let mut versions = VERSIONS
        .iter()
        .map(|(matcher, version)| {
            (
                Regex::new(
                    &matcher
                        .replace("(", "\\(")
                        .replace("{", "\\{")
                        .replace("*", ".*"),
                )
                .expect("failed to create regex"),
                version,
                0,
                matcher,
            )
        })
        .collect::<Vec<_>>();
    for line in content.iter_mut() {
        if !line.trim_start().starts_with("pub ") {
            continue;
        }
        if let Some(pos) = versions
            .iter()
            .position(|(matcher, _, _, _)| matcher.is_match(line))
        {
            let mut tmp_str = String::new();
            for c in line.chars() {
                if c != ' ' {
                    break;
                }
                tmp_str.push(' ');
            }
            println!("==> Adding feature for {:?}", versions[pos].3);
            line.insert_str(
                0,
                &format!("{}#[cfg(feature = \"{}\")]\n", tmp_str, versions[pos].1),
            );
            versions[pos].2 += 1;
        }
    }
    if versions.iter().any(|(_, _, count, _)| *count == 0) {
        println!("Some version features were not added:");
        for (_, _, _, v) in versions.iter().filter(|(_, _, count, _)| *count == 0) {
            println!("==> {:?}", v);
        }
    }
}

fn clean_structs(content: &mut Vec<String>) {
    for line in content.iter_mut() {
        if !line.starts_with("pub struct ") {
            continue;
        }
        // remove " {" at the end of the struct.
        line.pop();
        line.pop();
        line.push(';');
    }
}

fn run_bindgen(folder: &Path, commit_hash: String) {
    println!("=> Running bindgen...");
    let bindings = bindgen::Builder::default()
        .header(HEADER_FILE)
        .layout_tests(false)
        .size_t_is_usize(true)
        .clang_args(&[format!("-I{}", folder.display())])
        .whitelist_function("(gsl|cblas)_.*")
        .whitelist_type("(gsl|cblas)_.*")
        .whitelist_var("(GSL|CBLAS|gsl|cblas)_.*")
        .generate()
        .expect("Unable to generate bindings");

    println!("<= Done");

    let content = bindings.to_string();
    let mut content = first_pass(content.lines().collect::<Vec<_>>());
    clean_structs(&mut content);
    add_features(&mut content);

    let out = "../src/auto.rs";
    println!("=> Writing content into `{}`...", out);

    let mut f = OpenOptions::new()
        .truncate(true)
        .create(true)
        .write(true)
        .open(out)
        .expect("Failed to open binding file...");
    writeln!(
        f,
        "// Generated on commit {} from {}",
        commit_hash, GSL_REPOSITORY
    )
    .unwrap();
    writeln!(f, "// DO NOT EDIT THIS FILE!!!",).unwrap();
    writeln!(f, "").unwrap();
    write!(
        f,
        "{}",
        content
            .iter()
            .map(|x| {
                if x.starts_with("pub const GSL_ODEIV_HADJ_INC: u32 =") {
                    format!("{}\n", x.replace(": u32 =", ": i32 ="))
                } else if x.starts_with("pub const GSL_ODEIV_HADJ_NIL: u32 =") {
                    format!("{}\n", x.replace(": u32 =", ": i32 ="))
                } else {
                    format!("{}\n", x)
                }
            })
            .collect::<String>()
    )
    .unwrap();

    println!("<= Done");
}

fn ready_gsl_lib(folder: &Path) {
    if Command::new("git")
        .arg("clone")
        .arg(GSL_REPOSITORY)
        .arg("--depth")
        .arg("1")
        .arg(
            folder
                .join("gsl")
                .to_str()
                .expect("failed to convert path to str"),
        )
        .status()
        .is_err()
    {
        panic!("Failed to clone gsl repository...");
    }
    if Command::new("bash")
        .arg("-c")
        .arg(&format!("cd {}/gsl && ./autogen.sh", folder.display()))
        .status()
        .is_err()
    {
        panic!("Failed to run autogen.sh");
    }
    if Command::new("bash")
        .arg("-c")
        .arg(&format!("cd {}/gsl && ./configure", folder.display()))
        .status()
        .is_err()
    {
        panic!("Failed to run configure");
    }
    if Command::new("bash")
        .arg("-c")
        .arg(&format!("cd {}/gsl && make", folder.display()))
        .status()
        .is_err()
    {
        panic!("Failed to run make");
    }
}

fn get_current_commit_hash(folder: &Path) -> String {
    let commit_hash = Command::new("bash")
        .arg("-c")
        .arg(&format!(
            "cd {} && git rev-parse --short HEAD",
            folder.display()
        ))
        .output()
        .expect("Failed to retrieve current gsl commit hash");
    if !commit_hash.status.success() {
        panic!("Commit hash retrieval failed....");
    }
    String::from_utf8(commit_hash.stdout)
        .expect("Invalid commit hash received...")
        .trim()
        .to_owned()
}

fn run_everything(folder: &Path, ready_gsl: bool) {
    if ready_gsl {
        ready_gsl_lib(folder);
    }
    create_header_file(folder);
    run_bindgen(folder, get_current_commit_hash(folder));
}

fn main() {
    if env::args().skip(1).count() != 0 {
        let dir = env::args().skip(1).next().unwrap();
        println!(
            "Using `{}` path as gsl directory. No initialization will be performed on it",
            dir
        );

        run_everything(&Path::new(&dir), false);
        return;
    }

    let dir = tempfile::tempdir().expect("failed to create temporary directory");
    println!("Created temporary directory: {:?}", dir.path());

    run_everything(&dir.path(), true);
}
