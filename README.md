rust-GSL [![Build Status](https://api.travis-ci.org/GuillaumeGomez/rust-GSL.png?branch=master)](https://travis-ci.org/GuillaumeGomez/rust-GSL)
========

A __Rust__ binding for the [GSL library] (the GNU Scientific Library).

## Installation

This binding requires the [GSL library] library to be installed:

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

Instructions from https://www.gnu.org/software/gsl/extras/native_win_builds.html:

```
Building GSL on Windows Using Native Tools

Several packages are available for building GSL on Windows using Microsoft Visual Studio with the native Microsoft and Intel compilers.

1. Using Visual Studio Build Files

These files are designed for building the repository version of GSL using Visual Studio 2015. No additional build tools are required.

     https://github.com/BrianGladman/gsl

2. Using Perl Generated Visual Studio Project Files

This approach, which offers project files for Visual Studio 2010 and 2012, has been tested with GSL versions from 1.8 to 2.3. This can also be used to obtain build files for later version of Visual Studio since these can load files intended for use in the earlier versions of Visual Studio. This approach requires that Perl is installed.

     https://github.com/aweatherguy/GSL-on-Windows

3. Using Cmake

Cmake can produce build files for most versions of Visual Studio and for MSBUILD. A Cmake installation is required

     https://github.com/ampl/gsl/commit/0cf47581583de7ed4ee1cf0a095066af0d2f97d1

This approach will be especially useful for Visual Studio 2017, which has built-in support for Cmake build files.

4. Installation using NuGet

For those using NuGet there are win32 and x64 installation packages at:

     https://www.nuget.org/packages/gsl-msvc14-x86/

     https://www.nuget.org/packages/gsl-msvc14-x64/
```

This crate works with Cargo and is on [crates.io]. Just add the
following to your `Cargo.toml` file:

```toml
[dependencies]
GSL = "1.1"
```

Add the following line to your source code:

```rust
extern crate rgsl;
```

## Building

Just run `cargo build`. However, if you have issues with GSL version 2.x, try building it with:

```bash
cargo build --features v2
```

If a project depends on this version, don't forget to add in your `Cargo.toml`:

```toml
[features]
v2 = ["GSL/v2"]
```

If you always want to run the v2 features you can, instead, use the following in your `Cargo.toml`: 

```toml
[features]
default = ["GSL/v2"]
```

## Documentation

You can access the __rgsl__ documentation locally, just build it:

```Shell
> cargo doc --open
```

Then open this file with an internet browser: `file:///{rgsl_location}/target/doc/rgsl/index.html`

You can also access the latest build of the documentation via the internet [here](https://docs.rs/crate/GSL/).

## Donations

If you appreciate my work and want to support me, you can do it here:

[![Become a patron](https://c5.patreon.com/external/logo/become_a_patron_button.png)](https://www.patreon.com/GuillaumeGomez)

## License

__rust-GSL__ is a wrapper for __GSL__, therefore inherits the [GPL license](http://www.gnu.org/copyleft/gpl.html).

[crates.io]: https://crates.io/crates/GSL
[GSL library]: http://www.gnu.org/software/gsl/
