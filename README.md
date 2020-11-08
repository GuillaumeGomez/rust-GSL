# rust-GSL

A __Rust__ binding for the [GSL library][GSL library] (the GNU Scientific Library).

## Installation

This binding requires the [GSL library] library (version >= 2) to be installed:

### Linux

```bash
# on debian based systems:
sudo apt-get install libgsl0-dev
```

### OSX

```bash
brew install gsl
```

### Windows

Instructions are available there: <https://www.gnu.org/software/gsl/extras/native_win_builds.html>.

## Usage

This crate works with Cargo and is on [crates.io]. Just add the following to your `Cargo.toml` file:

```toml
[dependencies]
GSL = "2.0"
```

You can see examples in the `examples` folder.

## Building

To build `rgsl`, just run `cargo build`. However, if you want to use a specific version, you'll
need to use the `cargo` features. For example:

```bash
cargo build --features v2_1
```

If a project depends on this version, don't forget to add in your `Cargo.toml`:

```toml
[dependencies.GSL]
version = "2"
features = ["v2_1"]
```

## Documentation

You can access the __rgsl__ documentation locally, just build it:

```Shell
> cargo doc --open
```

Then open this file with an internet browser: `file:///{rgsl_location}/target/doc/rgsl/index.html`

You can also access the latest build of the documentation via the internet [here](https://docs.rs/crate/GSL/).

## License

__rust-GSL__ is a wrapper for __GSL__, therefore inherits the [GPL license](http://www.gnu.org/copyleft/gpl.html).

[crates.io]: https://crates.io/crates/GSL
[GSL library]: http://www.gnu.org/software/gsl/
