# `rgsl` examples

This folder contains the `rgsl` examples. To run one, just use `cargo`:

```bash
$ cargo run --bin intro
```

And that's it!

Some examples might require a higher GSL version. `rgsl` supports versions through features.
So for example:

```bash
$ cargo run --bin largefit --features v2_2
```

# Original examples

I ported some of the original GSL examples. You can find the original ones from the
[git repository](git://git.savannah.gnu.org/gsl.git) in the folder `doc/examples`.
