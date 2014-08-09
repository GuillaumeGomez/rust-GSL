//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub use self::mode::Mode;
pub use self::cblas_order::CblasOrder;
pub use self::cblas_side::CblasSide;
pub use self::cblas_transpose::CblasTranspose;
pub use self::cblas_uplo::CblasUplo;
pub use self::cblas_diag::CblasDiag;
pub use gsl_value::GslValue;

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

/// Used with _e functions
pub mod gsl_value {
    #[deriving(PartialEq, PartialOrd, Show)]
    #[repr(C)]
    pub enum GslValue {
        Success = 0,
        Failure = -1,
        /// iteration has not converged
        Continue = -2,
        /// input domain error, e.g sqrt(-1)
        Dom = 1,
        /// output range error, e.g. exp(1e100)
        Range = 2,
        /// invalid pointer
        Fault = 3,
        /// invalid argument supplied by user
        Inval = 4,
        /// generic failure
        Failed = 5,
        /// factorization failed
        Factor = 6,
        /// sanity check failed - shouldn't happen
        Sanity = 7,
        /// malloc failed
        NoMem = 8,
        /// problem with user-supplied function
        BadFunc = 9,
        /// iterative process is out of control
        RunAway = 10,
        /// exceeded max number of iterations
        MaxIter = 11,
        /// tried to divide by zero
        ZeroDiv = 12,
        /// user specified an invalid tolerance
        BadTol = 13,
        /// failed to reach the specified tolerance
        Tol = 14,
        /// underflow
        UndrFlw = 15,
        /// overflow
        OvrFlw = 16,
        /// loss of accuracy
        Loss = 17,
        /// failed because of roundoff error
        Round = 18,
        /// matrix, vector lengths are not conformant
        BadLen = 19,
        /// matrix not square
        NotSqr = 20,
        /// apparent singularity detected
        Sing = 21,
        /// integral or series is divergent
        Diverge = 22,
        /// requested feature is not supported by the hardware
        Unsup = 23,
        /// requested feature not (yet) implemented
        Unimpl = 24,
        /// cache limit exceeded
        Cache = 25,
        /// table limit exceeded
        Table = 26,
        /// iteration is not making progress towards solution
        NoProg = 27,
        /// jacobian evaluations are not improving the solution
        NoProgJ = 28,
        /// cannot reach the specified tolerance in F
        TolF = 29,
        /// cannot reach the specified tolerance in X
        TolX = 30,
        /// cannot reach the specified tolerance in gradient
        TolG = 31,
        /// cannot reach the specified tolerance in gradient
        EOF = 32
    }
}