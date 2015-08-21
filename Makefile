#
# A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
#

all: rgsl examples doc

rgsl:
	mkdir -p lib
	rustc --out-dir=lib src/rgsl.rs

examples: rgsl
	  mkdir -p bin
	  rustc -o bin/simple -L ./lib examples/simple.rs
	  rustc -o bin/rng -L ./lib examples/rng.rs
	  rustc -o bin/random_number_distribution -L ./lib examples/random_number_distribution.rs
	  rustc -o bin/permutation -L ./lib examples/permutation.rs
	  rustc -o bin/sorting -L ./lib examples/sorting.rs
	  rustc -o bin/chebyshev_approximation -L ./lib examples/chebyshev_approximation.rs
	  rustc -o bin/combination_example -L ./lib examples/combination.rs
	  rustc -o bin/roots_of_polynomial -L ./lib examples/roots_of_polynomial.rs
	  rustc -o bin/numerical_differentiation -L ./lib examples/numerical_differentiation.rs
	  rustc -o bin/radix -L ./lib examples/radix.rs
	  rustc -o bin/integration -L ./lib examples/integration.rs
	  rustc -o bin/interpolation -L ./lib examples/interpolation.rs
	  rustc -o bin/linear_algebra -L ./lib examples/linear_algebra.rs
	  rustc -o bin/minimization -L ./lib examples/minimization.rs
	  rustc -o bin/monte_carlo -L ./lib examples/monte_carlo.rs
	  rustc -o bin/n_tuples -L ./lib examples/n_tuples.rs
	  rustc -o bin/multisets -L ./lib examples/multisets.rs
	  rustc -o bin/ordinary_differential_equations -L ./lib examples/ordinary_differential_equations.rs
	  rustc -o bin/qrng -L ./lib examples/qrng.rs
	  rustc -o bin/statistics -L ./lib examples/statistics.rs
	  rustc -o bin/series_acceleration -L ./lib examples/series_acceleration.rs
	  rustc -o bin/wavelet_transforms -L ./lib examples/wavelet_transforms.rs
	  rustc -o bin/physical_constant -L ./lib examples/physical_constant.rs
	  rustc -o bin/vectors_and_matrices -L ./lib examples/vectors_and_matrices.rs
	  rustc -o bin/multifit_solver -L ./lib examples/multifit_solver.rs

doc:
	rustdoc -o doc src/rgsl.rs

clean:
	rm -rf lib
	rm -rf bin

re: clean all
