//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub use self::mode::Mode;
pub use self::cblas_order::CblasOrder;
pub use self::cblas_side::CblasSide;
pub use self::cblas_transpose::CblasTranspose;
pub use self::cblas_uplo::CblasUplo;
pub use self::cblas_diag::CblasDiag;

pub mod mode {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum Mode {
        PrecDouble,
        PrecSingle,
        PrecApprox
    }
}

/// Indicates whether a matrix is in Row Major or Column Major order.
/// Row major order is the native order for C programs, while Column major order is native for Fortran.
pub mod cblas_order {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum CblasOrder {
        RowMajor = 101,
        ColMajoyr
    }
}

/// Used to indicate the order of a matrix-matrix multiply.
pub mod cblas_side {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum CblasSide {
        /// Means __A__ __B__
        Left = 141,
        /// Means __B__ __A__
        Right
    }
}

/// Used to represent transpose operations on a matrix.
pub mod cblas_transpose {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum CblasTranspose {
        /// Represents __X__
        NoTrans = 111,
        /// Represents __X^T__
        Trans,
        /// Represents __X^H__
        ConjTrans
    }
}

/// Used to indicate which part of a symmetric matrix to use.
pub mod cblas_uplo {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum CblasUplo {
        /// Means user the upper triagle of the matrix.
        Upper = 121,
        /// Means use the lower triange of the matrix.
        Lower
    }
}

pub mod cblas_diag {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum CblasDiag {
        NonUnit = 131,
        Unit
    }
}