# pylint: disable=C0103,C0114,C0115,C0116,C0301,R0913

from os import walk
from os.path import join
import sys


MACROS_TO_IGNORE = ["ffi_wrap", "wrap_callback", "ffi_wrapper"]
FUNC_NAME_TO_IGNORE = [
    "new",
    "new_with_init",
    "from_slice",
    "nm_simplex",
    "nm_simplex2",
    "nm_simplex2_rand",
]


def read_file(path):
    with open(path, 'r') as fd:
        return fd.read()


def read_dirs(dir_path, excluded, errors, totals, functions_to_call):
    for root, _, files in walk(dir_path):
        for file in files:
            file_path = join(root, file)
            if file_path in excluded:
                continue
            content = read_file(file_path)
            for func in functions_to_call:
                func(file_path, content, errors, totals)


def add_error(file_path, line_nb, errors, err):
    errors.append("[{}:{}] => {}".format(file_path, line_nb, err))


def to_rust_name(n):
    if n.startswith("_"):
        n = n[1:]
    # Simple/stupid rule...
    if n.startswith("is"):
        n = "is_" + n[2:]
    return {
        "is_nonneg": "is_non_neg",
        "swap_rowcol": "swap_row_col",
    }.get(n, n)


def validate_name(rust_name, pending_func_name):
    if pending_func_name == "drop":
        return True
    if rust_name not in ["alloc", "calloc"]:
        for start in ["from_", "new"]:
            if pending_func_name.startswith(start):
                return True
    if (pending_func_name in "copy" or pending_func_name.startswith("copy_")) and rust_name.endswith("memcpy"):
        return True
    return rust_name.startswith("gsl_") and pending_func_name in rust_name


def init_check_macro_data(data, totals, ignored):
    data["pending_func_name"] = None
    data["sys_names"] = []
    data["pending_doc_aliases"] = []
    data["indent"] = 0
    data["func_line"] = 0
    added = 0
    if data.get("is_test", False):
        totals["test_count"] += 1
        ignored = True
    data["is_test"] = False
    if not ignored and data.get("is_in_macro", False):
        totals["in_macro_count"] += 1
        added += 1
    data["is_in_macro"] = False
    if not ignored and data.get("is_ffi_wrap", False):
        totals["ffi_wrap_count"] += 1
        added += 1
    data["is_ffi_wrap"] = False
    if not ignored and data.get("is_in_func", False):
        totals["in_func_count"] += 1
        added += 1
    data["is_in_func"] = False
    if not ignored and data.get("ignore_next", False):
        totals["ignored"] += 1
    if not ignored and added == 0:
        totals["rust_func_count"] += 1
    data["ignore_next"] = False


def check_macro_names(file_path, errors, data):
    if len(data["sys_names"]) < 1:
        # Nothing to check if there is no sys call.
        return
    # For now, we only check the last sys name.
    sys_name = data["sys_names"][-1]
    if len(data["pending_doc_aliases"]) == 0:
        add_error(
            file_path,
            data["func_line"],
            errors,
            "Missing `#[doc(alias = {})]`".format(sys_name))
    elif sys_name not in data["pending_doc_aliases"]:
        if len(data["pending_doc_aliases"]) == 1:
            add_error(
                file_path,
                data["func_line"],
                errors,
                "Mismatch between doc alias and sys call: `{}` != `{}`".format(
                    data["pending_doc_aliases"][0], sys_name))
        else:
            add_error(
                file_path,
                data["func_line"],
                errors,
                "Mismatch between doc alias and sys call: no doc alias named `{}` in {}".format(
                    sys_name, data["pending_doc_aliases"]))
    if data["ignore_next"]:
        return
    rust_name = to_rust_name(sys_name.split(" ")[-1])
    if rust_name != data["pending_func_name"]:
        if (data["pending_func_name"] not in FUNC_NAME_TO_IGNORE
                and not validate_name(rust_name, data["pending_func_name"])):
            add_error(
                file_path,
                data["func_line"],
                errors,
                "Mismatch between function name and sys name: should be `{}` (but currently is `{}`)".format(
                    rust_name, data["pending_func_name"]))


def check_counts(errors, data):
    if data["in_macro_count"] < 1:
        errors.append("No checks run on macros...")
    else:
        print("Checked {} items in macros".format(data["in_macro_count"]))
    if data["in_func_count"] < 1:
        errors.append("No checks run on functions (outside of macros!)...")
    else:
        print("Checked {} items not in macros".format(data["in_func_count"]))
    if data["ffi_wrap_count"] < 1:
        errors.append("No checks run on `ffi_wrap` macro...")
    else:
        print("Checked {} `ffi_wrap` items".format(data["ffi_wrap_count"]))
    print("Rust functions: {}".format(data["rust_func_count"]))
    print("Rust test functions: {}".format(data["test_count"]))
    print("Ignored items: {}".format(data["ignored"]))


# This test is used to ensure that the methods generated through macros are correct (doc alias and
# FFI call conform to what they should be).
# pylint: disable=R0912
def check_macros(file_path, content, errors, totals):
    is_in_macro = False
    is_in_comment = False
    is_in_decl = False
    data = {}
    init_check_macro_data(data, totals, True)
    for line_nb, line in enumerate(content.split('\n')):
        if len(line) == 0:
            continue
        stripped_line = line.strip()
        # comment part
        if stripped_line.startswith("/*"):
            is_in_comment = True
        if is_in_comment:
            if line.endswith("*/"):
                is_in_comment = False
            continue
        if stripped_line.startswith("//"):
            if stripped_line.split("//")[1].strip() == "checker:ignore":
                data["ignore_next"] = True
            continue
        # macro parsing part
        if not is_in_macro:
            if stripped_line.startswith("macro_rules! "):
                macro_name = stripped_line.split("macro_rules! ")[1].split("{")[0]
                if not macro_name in MACROS_TO_IGNORE:
                    is_in_macro = True
                continue
        if not is_in_decl:
            if (stripped_line.startswith("pub trait ") or stripped_line.startswith("trait ")
                    or stripped_line.startswith("pub struct ") or stripped_line.startswith("struct ")):
                is_in_decl = True
                continue
        if stripped_line == "}":
            if line == "}":
                is_in_macro = False
                is_in_decl = False
            current_indent = len(line) - len(line.lstrip())
            if current_indent == data["indent"]:
                ignored = True
                if data["pending_func_name"] is not None and not is_in_decl:
                    check_macro_names(file_path, errors, data)
                    ignored = False
                init_check_macro_data(data, totals, ignored)
        if is_in_decl:
            continue
        elif " ffi_wrap!(" in line:
            sys_name = stripped_line.split("ffi_wrap!(")[1].split(")")[0].strip()
            data["sys_names"].append(sys_name)
            data["is_ffi_wrap"] = True
        else:
            if stripped_line.startswith("#[doc(alias"):
                alias_name = stripped_line.split("(alias")[1].split("=")[1].split(")]")[0].strip()
                data["pending_doc_aliases"].append(alias_name.replace("\"", ""))
            if stripped_line == "#[test]":
                data["is_test"] = True
            elif data["pending_func_name"] is None and (stripped_line.startswith("pub fn ") or stripped_line.startswith("fn ")):
                data["pending_func_name"] = line.split("(")[0].split("<")[0].split("fn ")[1].strip()
                data["indent"] = len(line) - len(line.lstrip())
                data["func_line"] = line_nb + 1
            elif data["pending_func_name"] is not None:
                if " sys::[<" in stripped_line:
                    # A sys call in a macro.
                    sys_name = stripped_line.split(" sys::[<")[1].split(">](")[0].strip()
                    data["sys_names"].append(sys_name)
                    data["is_in_macro"] = True
                elif "sys::gsl_" in stripped_line:
                    # A sys call NOT in a macro.
                    part = "gsl_" + stripped_line.split("sys::gsl_")[1]
                    if "(" in part:
                        sys_name = part.split("(")[0].split("{")[0].strip()
                        data["sys_names"].append(sys_name)
                        data["is_in_func"] = True

def check_file_header(file_path, content, errors, _totals):
    if content.split('\n')[:3] != [
        "//",
        "// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)",
        "//"
    ]:
        add_error(file_path, 1, errors, "Invalid header (take a look at `lib.rs`)")


def main():
    errors = []
    totals = {
        "in_macro_count": 0,
        "ffi_wrap_count": 0,
        "in_func_count": 0,
        "rust_func_count": 0,
        "test_count": 0,
        "ignored": 0,
    }
    excluded = ["src/_docs/header.html", "src/integration.rs"]
    read_dirs("src", excluded, errors, totals, [check_file_header, check_macros])
    check_counts(errors, totals)
    if len(errors) > 0:
        for err in errors:
            print("=> {}".format(err))
        print()
        print("{} error{} occurred".format(len(errors), "s" if len(errors) > 1 else "" ))
    return 0 if len(errors) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
