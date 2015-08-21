rust-GSL [![Build Status](https://api.travis-ci.org/GuillaumeGomez/rust-GSL.png?branch=master)](https://travis-ci.org/GuillaumeGomez/rust-GSL) [![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/GuillaumeGomez/rust-GSL?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
========

A __Rust__ binding for the [GSL library] (the GNU Scientific Library).

##Installation

This binding requires the [GSL library] library to be installed.

This crate works with Cargo and is on [crates.io].  Just add the
following to your `Cargo.toml` file:

```toml
[dependencies]
GSL = "*"
```

Add the following line to your source code:

```rust
extern crate rgsl;
```

##Documentation

You can access the __rgsl__ documentation locally, just build it:

```Shell
> cargo doc --open
```

Then open this file with an internet browser:
`file:///{rgsl_location}/target/doc/rgsl/index.html`

You can also access the latest build of the documentation via the internet [here](http://rust-ci.org/GuillaumeGomez/rust-GSL/doc/rgsl/).

## License
__rust-GSL__ is a wrapper for __GSL__, therefore inherits the [GPL license](http://www.gnu.org/copyleft/gpl.html).

[crates.io]: https://crates.io/crates/GSL
[GSL library]: http://www.gnu.org/software/gsl/
