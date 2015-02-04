//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub mod mode {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    pub enum Mode {
        PrecDouble,
        PrecSingle,
        PrecApprox
    }
}

pub mod value {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
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
}

pub mod eigen_sort {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    pub enum EigenSort {
        /// ascending order in numerical value
        ValAsc,
        /// descending order in numerical value
        VasDesc,
        /// ascending order in magnitude
        AbsAsc,
        /// descending order in magnitude
        AbsDesc
    }
}

pub mod fft_direction {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    /// this gives the sign in the formula
    /// 
    /// h(f) = \sum x(t) exp(+/- 2 pi i f t)
    ///     
    /// where - is the forward transform direction and + the inverse direction
    pub enum FftDirection {
        Forward = -1,
        Backward = 1
    }
}

pub mod gauss_konrod_rule {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    /// The low-level integration rules in QUADPACK are identified by small integers (1-6). We'll use symbolic constants to refer to them.
    pub enum GaussKonrodRule {
        /// 15 point Gauss-Kronrod rule
        Gauss15 = 1,
        /// 21 point Gauss-Kronrod rule
        Gauss21 = 2,
        /// 31 point Gauss-Kronrod rule
        Gauss31 = 3,
        /// 41 point Gauss-Kronrod rule
        Gauss41 = 4,
        /// 51 point Gauss-Kronrod rule
        Gauss51 = 5,
        /// 61 point Gauss-Kronrod rule
        Gauss61 = 6,
    }
}

pub mod integration_qawo {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    /// Used by workspace for QAWO integrator
    pub enum IntegrationQawo {
        Cosine,
        Sine
    }
}

pub mod vegas_mode {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    /// Used by VegasMonteCarlo struct
    pub enum VegasMode {
        Importance = 1, 
        ImportanceOnly = 0, 
        Stratified = -1
    }
}

pub mod odeiv {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    /// Possible return values for an hadjust() evolution method for ordinary differential equations
    pub enum ODEiv {
        /// step was increased
        Inc = 1,
        /// step unchanged
        Nil = 0,
        /// step decreased
        Dec = -1
    }
}

pub mod wavelet_direction {
    #[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
    #[repr(C)]
    pub enum WaveletDirection {
        Forward = 1,
        Backward = -1
    }
}