rust-GSL [![Build Status](https://api.travis-ci.org/GuillaumeGomez/rust-GSL.png?branch=master)](https://travis-ci.org/GuillaumeGomez/rust-GSL) [![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/GuillaumeGomez/rust-GSL?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
========

A __Rust__ binding for the [GSL library](http://www.gnu.org/software/gsl/) (the GNU Scientific Library).

##Installation

This binding requires the [GSL library](http://www.gnu.org/software/gsl/) library to be installed.

To build it, please use :

```Shell
> make
```

This command build __rgsl__, the examples and the documentation.

You can build them separatly too.

```Shell
> make rgsl
> make examples
> make doc
```

Since this project supports cargo, you can also build it like this :

```Shell
> cargo build
```

##Documentation

You can access the __rgsl__ documentation locally, just build it :

```Shell
> make doc
```

Then open this file with an internet browser :
file:///{rgsl_location}/doc/rgsl/index.html

You can also access the latest build of the documentation via the internet [here](http://rust-ci.org/GuillaumeGomez/rust-GSL/doc/rgsl/).

## License
__rust-GSL__ is a wrapper for __GSL__, therefore inherits the [GPL licencs](http://www.gnu.org/copyleft/gpl.html).
