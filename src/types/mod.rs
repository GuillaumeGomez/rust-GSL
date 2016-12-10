//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

#[cfg(not(feature = "v2"))]
pub use self::basis_spline::{
    BSpLineWorkspace,
    BSpLineDerivWorkspace
};
#[cfg(feature = "v2")]
pub use self::basis_spline::BSpLineWorkspace;

pub use self::chebyshev::ChebSeries;
pub use self::combination::Combination;
pub use self::complex::{ComplexF32, ComplexF64};
pub use self::discrete_hankel::DiscreteHankel;
pub use self::eigen_symmetric_workspace::{EigenSymmetricWorkspace, EigenSymmetricVWorkspace, EigenHermitianWorkspace,
    EigenHermitianVWorkspace, EigenNonSymmWorkspace, EigenNonSymmVWorkspace, EigenGenSymmWorkspace, EigenGenSymmVWorkspace,
    EigenGenHermWorkspace, EigenGenHermVWorkspace, EigenGenWorkspace, EigenGenVWorkspace};
pub use self::fast_fourier_transforms::{FftComplexWaveTable, FftComplexWorkspace};
pub use self::histograms::{Histogram, HistogramPdf, Histogram2D, Histogram2DPdf};
pub use self::integration::{IntegrationWorkspace, IntegrationQawsTable, IntegrationQawoTable, CquadWorkspace, GLFixedTable};
pub use self::interpolation::{InterpAccel, Interp, InterpType, Spline};
pub use self::mathieu::MathieuWorkspace;
pub use self::matrix::{MatrixF32, MatrixF64, MatrixView};
pub use self::matrix_complex::{MatrixComplexF32, MatrixComplexF64};
pub use self::minimizer::{Minimizer, MinimizerType};
pub use self::monte_carlo::{PlainMonteCarlo, MiserMonteCarlo, MiserParams, VegasMonteCarlo};
pub use self::multifit_solver::{MultiFitFdfSolver, MultiFitFunction, MultiFitFdfSolverType, MultiFitFunctionFdf};
pub use self::multiset::MultiSet;
pub use self::n_tuples::NTuples;
pub use self::ordinary_differential_equations::{ODEiv2System, ODEiv2Step, ODEiv2StepType, ODEiv2Control, ODEiv2Evolve, ODEiv2Driver};
pub use self::permutation::Permutation;
pub use self::polynomial::PolyComplex;
pub use self::qrng::{QRng, QRngType};
pub use self::ran_discrete::RanDiscrete;
pub use self::result::{Result, ResultE10};
pub use self::rng::{Rng, RngType};
pub use self::series_acceleration::{LevinUWorkspace, LevinUTruncWorkspace};
pub use self::vector::{VectorF32, VectorF64, VectorView};
pub use self::vector_complex::{VectorComplexF32, VectorComplexF64};
pub use self::wavelet_transforms::{Wavelet, WaveletType, WaveletWorkspace};

pub mod basis_spline;
pub mod chebyshev;
pub mod combination;
pub mod complex;
pub mod discrete_hankel;
pub mod eigen_symmetric_workspace;
pub mod fast_fourier_transforms;
pub mod histograms;
pub mod integration;
pub mod interpolation;
pub mod mathieu;
pub mod matrix;
pub mod matrix_complex;
pub mod minimizer;
pub mod monte_carlo;
pub mod multifit_solver;
pub mod multiset;
pub mod n_tuples;
pub mod ordinary_differential_equations;
pub mod permutation;
pub mod polynomial;
pub mod qrng;
pub mod ran_discrete;
pub mod result;
pub mod series_acceleration;
pub mod rng;
pub mod vector;
pub mod vector_complex;
pub mod wavelet_transforms;
