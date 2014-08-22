//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub use self::basis_spline::{BSpLineWorkspace, BSpLineDerivWorkspace};
pub use self::complex::{ComplexF32, ComplexF64};
pub use self::matrix::{MatrixF32, MatrixF64};
pub use self::matrix_complex::{MatrixComplexF32, MatrixComplexF64};
pub use self::result::{Result, ResultE10};
pub use self::rng::{Rng, RngType};
pub use self::vector::{VectorF32, VectorF64};
pub use self::vector_complex::{VectorComplexF32, VectorComplexF64};
pub use self::mathieu::MathieuWorkspace;
pub use self::ran_discrete::RanDiscrete;
pub use self::permutation::Permutation;
pub use self::chebyshev::ChebSeries;
pub use self::combination::Combination;

pub mod basis_spline;
pub mod chebyshev;
pub mod combination;
pub mod complex;
pub mod mathieu;
pub mod matrix;
pub mod matrix_complex;
pub mod permutation;
pub mod ran_discrete;
pub mod result;
pub mod rng;
pub mod vector;
pub mod vector_complex;