//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod Gsl {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum Mode {
        PrecDouble,
        PrecSingle,
        PrecApprox
    }

    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    /// Indicates whether a matrix is in Row Major or Column Major order.
    /// Row major order is the native order for C programs, while Column major order is native for Fortran.
    pub enum CblasOrder {
        RowMajor = 101,
        ColMajoyr
    }

    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    /// Used to indicate the order of a matrix-matrix multiply.
    pub enum CblasSide {
        /// Means __A__ __B__
        Left = 141,
        /// Means __B__ __A__
        Right
    }

    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    /// Used to represent transpose operations on a matrix.
    pub enum CblasTranspose {
        /// Represents __X__
        NoTrans = 111,
        /// Represents __X^T__
        Trans,
        /// Represents __X^H__
        ConjTrans
    }

    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    /// Used to indicate which part of a symmetric matrix to use.
    pub enum CblasUplo {
        /// Means user the upper triagle of the matrix.
        Upper = 121,
        /// Means use the lower triange of the matrix.
        Lower
    }

    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum CblasDiag {
        NonUnit = 131,
        Unit
    }
}