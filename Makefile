#
# A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
#

all: rgsl examples doc

rgsl:
	mkdir -p lib
	rustc --out-dir=lib src/rgsl.rs

examples: rgsl
	  mkdir -p bin
	  rustc -o bin/simple -L ./lib examples/simple/main.rs
	  rustc -o bin/rng -L ./lib examples/rng/main.rs
	  rustc -o bin/random_number_distribution -L ./lib examples/random_number_distribution/main.rs
	  rustc -o bin/permutation -L ./lib examples/permutation/main.rs
	  rustc -o bin/sorting -L ./lib examples/sorting/main.rs
	  rustc -o bin/chebyshev_approximation -L ./lib examples/chebyshev_approximation/main.rs
	  rustc -o bin/combination_example -L ./lib examples/combination/main.rs
	  rustc -o bin/roots_of_polynomial -L ./lib examples/roots_of_polynomial/main.rs
	  rustc -o bin/numerical_differentiation -L ./lib examples/numerical_differentiation/main.rs
	  rustc -o bin/radix -L ./lib examples/radix/main.rs
	  rustc -o bin/integration -L ./lib examples/integration/main.rs
	  rustc -o bin/interpolation -L ./lib examples/interpolation/main.rs
	  rustc -o bin/linear_algebra -L ./lib examples/linear_algebra/main.rs
	  rustc -o bin/minimization -L ./lib examples/minimization/main.rs
	  rustc -o bin/monte_carlo -L ./lib examples/monte_carlo/main.rs
	  rustc -o bin/n_tuples -L ./lib examples/n_tuples/main.rs
	  rustc -o bin/multisets -L ./lib examples/multisets/main.rs
	  rustc -o bin/ordinary_differential_equations -L ./lib examples/ordinary_differential_equations/main.rs
	  rustc -o bin/qrng -L ./lib examples/qrng/main.rs
	  rustc -o bin/statistics -L ./lib examples/statistics/main.rs
	  rustc -o bin/series_acceleration -L ./lib examples/series_acceleration/main.rs
	  rustc -o bin/wavelet_transforms -L ./lib examples/wavelet_transforms/main.rs
	  rustc -o bin/physical_constant -L ./lib examples/physical_constant/main.rs
	  rustc -o bin/vectors_and_matrices -L ./lib examples/vectors_and_matrices/main.rs
	  rustc -o bin/multifit_solver -L ./lib examples/multifit_solver/main.rs

doc:
	rustdoc -o doc src/rgsl.rs

clean:
	rm -rf lib
	rm -rf bin

re: clean all
