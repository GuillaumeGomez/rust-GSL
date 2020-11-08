# gsl-sys

This is the FFI counter-part of the Rust GSL crate. It is better to use the GSL crate
directly rather than this one (unless something is missing the Rust binding!).

## Update FFI

Most of the FFI is generated using `bindgen`. Therefore, it'll be used to update. You'll also
need the `gsl` repository locally because we need its headers.

Basically, an update would look like this:

```bash
cd bin && cargo run
```

However, if you want to change the output a bit, it's **strongly** recommended that you
instead clone and setup the GSL repository on your computer directly:

```bash
git clone git://git.savannah.gnu.org/gsl.git
cd gsl
./autogen.sh
./configure
make
# The headers should now all be in the gsl subfolder!
```

Then run the FFI generation like this:

```bash
cd bin && cargo run -- [PATH TO GSL]
```

It'll prevent all the clone and rebuild every time.
