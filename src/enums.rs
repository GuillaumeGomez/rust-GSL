//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub type gsl_mode_t = u32;

/// The maximum x such that gamma(x) is not considered an overflow.
pub static SF_GAMMA_XMAX : f64 = 171.0;
/// The maximum n such that gsl_sf_fact(n) does not give an overflow.
pub static SF_FACT_NMAX : f64 = 170.0;
/// The maximum n such that gsl_sf_doublefact(n) does not give an overflow.
pub static SF_DOUBLEFACT_NMAX : f64 = 297.0;

#[deriving(PartialEq, PartialOrd, Show)]
#[repr(C)]
pub enum Mode {
    PrecDouble,
    PrecSingle,
    PrecApprox
}

#[deriving(PartialEq, PartialOrd, Show)]
#[repr(C)]
/// Used with _e functions
pub enum Value {
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