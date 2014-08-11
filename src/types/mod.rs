//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub use self::complex::{ComplexF32, ComplexF64};
pub use self::matrix::{MatrixF32, MatrixF64};
pub use self::matrix_complex::{MatrixComplexF32, MatrixComplexF64};
pub use self::result::{Result, ResultE10};
pub use self::vector::{VectorF32, VectorF64};
pub use self::vector_complex::{VectorComplexF32, VectorComplexF64};

pub mod complex;
pub mod matrix;
pub mod matrix_complex;
pub mod result;
pub mod vector;
pub mod vector_complex;