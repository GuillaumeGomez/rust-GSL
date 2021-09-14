# pylint: disable=C0103,C0114,C0115,C0116,C0301,R0913

from os import walk
from os.path import join
import sys


MACROS_TO_IGNORE = ["ffi_wrap", "wrap_callback", "ffi_wrapper"]
FUNC_NAME_TO_IGNORE = ["new", "new_with_init", "from_slice"]


def read_file(path):
    with open(path, 'r') as fd:
        return fd.read()


def read_dirs(dir_path, errors, functions_to_call):
    for root, _, files in walk(dir_path):
        for file in files:
            file_path = join(root, file)
            content = read_file(file_path)
            for func in functions_to_call:
                func(file_path, content, errors)


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
    if rust_name not in ["alloc", "calloc"]:
        for start in ["from_", "new"]:
            if pending_func_name.startswith(start):
                return True
    if (pending_func_name == "copy" or pending_func_name.startswith("copy_")) and rust_name.endswith("memcpy"):
        return True
    return rust_name.startswith("gsl_") and pending_func_name in rust_name


def init_check_macro_data(data):
    data["pending_func_name"] = None
    data["sys_names"] = []
    data["pending_doc_aliases"] = []
    data["indent"] = 0
    data["ignore_next"] = False
    data["func_line"] = 0


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


# This test is used to ensure that the methods generated through macros are correct (doc alias and
# FFI call conform to what they should be).
# pylint: disable=R0912
def check_macros(file_path, content, errors):
    is_in_macro = False
    is_in_comment = False
    data = {}
    init_check_macro_data(data)
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
        if line == "}":
            is_in_macro = False
        elif stripped_line == "}":
            current_indent = len(line) - len(line.lstrip())
            if current_indent == data["indent"]:
                if data["pending_func_name"] is not None:
                    check_macro_names(file_path, errors, data)
                init_check_macro_data(data)
        else:
            if stripped_line.startswith("#[doc(alias"):
                alias_name = stripped_line.split("(alias")[1].split("=")[1].split(")]")[0].strip()
                data["pending_doc_aliases"].append(alias_name.replace("\"", ""))
            elif stripped_line.startswith("pub fn "):
                data["pending_func_name"] = stripped_line.split("(")[0].split("<")[0].split("pub fn ")[1].strip()
                data["indent"] = len(line) - len(line.lstrip())
                data["func_line"] = line_nb + 1
            elif data["pending_func_name"] is not None:
                if " sys::[<" in stripped_line:
                    # A sys call in a macro.
                    sys_name = stripped_line.split(" sys::[<")[1].split(">](")[0].strip()
                    data["sys_names"].append(sys_name)
                elif "sys::gsl_" in stripped_line:
                    # A sys call NOT in a macro.
                    part = "gsl_" + stripped_line.split("sys::gsl_")[1]
                    if "(" in part:
                        sys_name = part.split("(")[0].split("{")[0].strip()
                        data["sys_names"].append(sys_name)


def main():
    errors = []
    read_dirs("src", errors, [check_macros])
    if len(errors) > 0:
        for err in errors:
            print("=> {}".format(err))
        print()
        print("{} error{} occurred".format(len(errors), "s" if len(errors) > 1 else "" ))
    return 0 if len(errors) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
