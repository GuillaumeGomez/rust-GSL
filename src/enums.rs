//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::ffi::CStr;
use std::os::raw::{c_char, c_int};

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Mode {
    PrecDouble,
    PrecSingle,
    PrecApprox,
}

#[allow(clippy::from_over_into)]
impl Into<sys::gsl_mode_t> for Mode {
    fn into(self) -> sys::gsl_mode_t {
        match self {
            Mode::PrecDouble => sys::GSL_PREC_DOUBLE,
            Mode::PrecSingle => sys::GSL_PREC_SINGLE,
            Mode::PrecApprox => sys::GSL_PREC_APPROX,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_mode_t> for Mode {
    fn from(v: sys::gsl_mode_t) -> Mode {
        match v {
            sys::GSL_PREC_DOUBLE => Mode::PrecDouble,
            sys::GSL_PREC_SINGLE => Mode::PrecSingle,
            sys::GSL_PREC_APPROX => Mode::PrecApprox,
            _ => panic!("Unknown Mode value"),
        }
    }
}

#[deprecated(since = "8.0.0", note = "Use rgsl::sf::Error instead")]
pub type Value = Error;

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Error {
    Failure,
    /// iteration has not converged
    Continue,
    /// input domain error, e.g sqrt(-1)
    Domain,
    /// output range error, e.g. exp(1e100)
    Range,
    /// invalid pointer
    Fault,
    /// invalid argument supplied by user
    Invalid,
    /// generic failure
    Failed,
    /// factorization failed
    Factorization,
    /// sanity check failed - shouldn't happen
    Sanity,
    /// malloc failed
    NoMemory,
    /// problem with user-supplied function
    BadFunction,
    /// iterative process is out of control
    RunAway,
    /// exceeded max number of iterations
    MaxIteration,
    /// tried to divide by zero
    ZeroDiv,
    /// user specified an invalid tolerance
    BadTolerance,
    /// failed to reach the specified tolerance
    Tolerance,
    /// underflow
    UnderFlow,
    /// overflow
    OverFlow,
    /// loss of accuracy
    Loss,
    /// failed because of roundoff error
    Round,
    /// matrix, vector lengths are not conformant
    BadLength,
    /// matrix not square
    NotSquare,
    /// apparent singularity detected
    Singularity,
    /// integral or series is divergent
    Diverge,
    /// requested feature is not supported by the hardware
    Unsupported,
    /// requested feature not (yet) implemented
    Unimplemented,
    /// cache limit exceeded
    Cache,
    /// table limit exceeded
    Table,
    /// iteration is not making progress towards solution
    NoProgress,
    /// jacobian evaluations are not improving the solution
    NoProgressJacobian,
    /// cannot reach the specified tolerance in F
    ToleranceF,
    /// cannot reach the specified tolerance in X
    ToleranceX,
    /// cannot reach the specified tolerance in gradient
    ToleranceG,
    /// cannot reach the specified tolerance in gradient
    #[allow(clippy::upper_case_acronyms)]
    EOF,
    /// Unknown value.
    Unknown(i32),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use Error::*;
        let err = match self {
            Failure => "Failure",
            Continue => "The iteration has not converged yet",
            Domain => "Input domain error",
            Range => "Output range error",
            Fault => "Invalid pointer",
            Invalid => "Invalid argument supplied by user",
            Failed => "generic failure",
            Factorization => "Factorization failed",
            Sanity => "Sanity check failed - shouldn't happen",
            NoMemory => "Malloc failed",
            BadFunction => "Problem with user-supplied function",
            RunAway => "Iterative process is out of control",
            MaxIteration => "Exceeded max number of iterations",
            ZeroDiv => "Tried to divide by zero",
            BadTolerance => {
                "Specified tolerance is invalid or \
                            theoretically unattainable"
            }
            Tolerance => "Failed to reach the specified tolerance",
            UnderFlow => "Underflow",
            OverFlow => "Overflow",
            Loss => "Loss of accuracy",
            Round => "Roundoff error",
            BadLength => "Matrix/vector sizes are not conformant",
            NotSquare => "Matrix not square",
            Singularity => {
                "Singularity or extremely bad function \
                           behavior detected"
            }
            Diverge => "Integral or series is divergent",
            Unsupported => {
                "The required feature is not supported by \
                           this hardware platform"
            }
            Unimplemented => "The requested feature is not (yet) implemented",
            Cache => "Cache limit exceeded",
            Table => "Table limit exceeded",
            NoProgress => "Iteration is not making progress towards solution",
            NoProgressJacobian => {
                "Jacobian evaluations are not improving \
                                   the solution"
            }
            ToleranceF => "Cannot reach the specified tolerance in F",
            ToleranceX => "Cannot reach the specified tolerance in X",
            ToleranceG => "Cannot reach the specified tolerance in gradient",
            EOF => "End of file",
            Unknown(_) => "Unknown error",
        };
        write!(f, "{}", err)
    }
}

impl std::error::Error for Error {}

impl Error {
    pub(crate) fn handle<T>(v: c_int, x: T) -> Result<T, Error> {
        match v {
            sys::GSL_SUCCESS => Ok(x),
            sys::GSL_FAILURE => Err(Self::Failure),
            sys::GSL_CONTINUE => Err(Self::Continue),
            sys::GSL_EDOM => Err(Self::Domain),
            sys::GSL_ERANGE => Err(Self::Range),
            sys::GSL_EFAULT => Err(Self::Fault),
            sys::GSL_EINVAL => Err(Self::Invalid),
            sys::GSL_EFAILED => Err(Self::Failed),
            sys::GSL_EFACTOR => Err(Self::Factorization),
            sys::GSL_ESANITY => Err(Self::Sanity),
            sys::GSL_ENOMEM => Err(Self::NoMemory),
            sys::GSL_EBADFUNC => Err(Self::BadFunction),
            sys::GSL_ERUNAWAY => Err(Self::RunAway),
            sys::GSL_EMAXITER => Err(Self::MaxIteration),
            sys::GSL_EZERODIV => Err(Self::ZeroDiv),
            sys::GSL_EBADTOL => Err(Self::BadTolerance),
            sys::GSL_ETOL => Err(Self::Tolerance),
            sys::GSL_EUNDRFLW => Err(Self::UnderFlow),
            sys::GSL_EOVRFLW => Err(Self::OverFlow),
            sys::GSL_ELOSS => Err(Self::Loss),
            sys::GSL_EROUND => Err(Self::Round),
            sys::GSL_EBADLEN => Err(Self::BadLength),
            sys::GSL_ENOTSQR => Err(Self::NotSquare),
            sys::GSL_ESING => Err(Self::Singularity),
            sys::GSL_EDIVERGE => Err(Self::Diverge),
            sys::GSL_EUNSUP => Err(Self::Unsupported),
            sys::GSL_EUNIMPL => Err(Self::Unimplemented),
            sys::GSL_ECACHE => Err(Self::Cache),
            sys::GSL_ETABLE => Err(Self::Table),
            sys::GSL_ENOPROG => Err(Self::NoProgress),
            sys::GSL_ENOPROGJ => Err(Self::NoProgressJacobian),
            sys::GSL_ETOLF => Err(Self::ToleranceF),
            sys::GSL_ETOLX => Err(Self::ToleranceX),
            sys::GSL_ETOLG => Err(Self::ToleranceG),
            sys::GSL_EOF => Err(Self::EOF),
            x => Err(Self::Unknown(x)),
        }
    }

    pub(crate) fn to_c(x: Result<(), Error>) -> c_int {
        match x {
            Ok(()) => sys::GSL_SUCCESS,
            Err(Error::Failure) => sys::GSL_FAILURE,
            Err(Error::Continue) => sys::GSL_CONTINUE,
            Err(Error::Domain) => sys::GSL_EDOM,
            Err(Error::Range) => sys::GSL_ERANGE,
            Err(Error::Fault) => sys::GSL_EFAULT,
            Err(Error::Invalid) => sys::GSL_EINVAL,
            Err(Error::Failed) => sys::GSL_EFAILED,
            Err(Error::Factorization) => sys::GSL_EFACTOR,
            Err(Error::Sanity) => sys::GSL_ESANITY,
            Err(Error::NoMemory) => sys::GSL_ENOMEM,
            Err(Error::BadFunction) => sys::GSL_EBADFUNC,
            Err(Error::RunAway) => sys::GSL_ERUNAWAY,
            Err(Error::MaxIteration) => sys::GSL_EMAXITER,
            Err(Error::ZeroDiv) => sys::GSL_EZERODIV,
            Err(Error::BadTolerance) => sys::GSL_EBADTOL,
            Err(Error::Tolerance) => sys::GSL_ETOL,
            Err(Error::UnderFlow) => sys::GSL_EUNDRFLW,
            Err(Error::OverFlow) => sys::GSL_EOVRFLW,
            Err(Error::Loss) => sys::GSL_ELOSS,
            Err(Error::Round) => sys::GSL_EROUND,
            Err(Error::BadLength) => sys::GSL_EBADLEN,
            Err(Error::NotSquare) => sys::GSL_ENOTSQR,
            Err(Error::Singularity) => sys::GSL_ESING,
            Err(Error::Diverge) => sys::GSL_EDIVERGE,
            Err(Error::Unsupported) => sys::GSL_EUNSUP,
            Err(Error::Unimplemented) => sys::GSL_EUNIMPL,
            Err(Error::Cache) => sys::GSL_ECACHE,
            Err(Error::Table) => sys::GSL_ETABLE,
            Err(Error::NoProgress) => sys::GSL_ENOPROG,
            Err(Error::NoProgressJacobian) => sys::GSL_ENOPROGJ,
            Err(Error::ToleranceF) => sys::GSL_ETOLF,
            Err(Error::ToleranceX) => sys::GSL_ETOLX,
            Err(Error::ToleranceG) => sys::GSL_ETOLG,
            Err(Error::EOF) => sys::GSL_EOF,
            Err(Error::Unknown(x)) => x,
        }
    }
}

// FIXME: Can do better?
static mut CALLBACK: Option<fn(&str, &str, u32, Error)> = None;

/// `f` is the type of GSL error handler functions. An error handler will be passed four arguments
/// which specify the reason for the error (a string), the name of the source file in which it
/// occurred (also a string), the line number in that file (an integer) and the error number (an
/// integer). The source file and line number are set at compile time using the __FILE__ and
/// __LINE__ directives in the preprocessor. An error handler function returns type void. Error
/// handler functions should be defined like this,
///
/// This function sets a new error handler, new_handler, for the GSL library routines. The previous
/// handler is returned (so that you can restore it later). Note that the pointer to a user defined
/// error handler function is stored in a static variable, so there can be only one error handler
/// per program. This function should be not be used in multi-threaded programs except to set up a
/// program-wide error handler from a master thread. The following example shows how to set and
/// restore a new error handler,
///
/// ```
/// use rgsl::{Error, set_error_handler};
///
/// fn error_handling(error_str: &str, file: &str, line: u32, error_value: Error) {
///     println!("[{:?}] '{}:{}': {}", error_value, file, line, error_str);
/// }
///
/// /* save original handler, install new handler */
/// let old_handler = set_error_handler(Some(error_handling));
///
/// /* code uses new handler */
/// // ...
///
/// /* restore original handler */
/// set_error_handler(old_handler);
/// ```
///
/// To use the default behavior (abort on error) set the error handler
/// to `None`:
///
/// ```
/// # use rgsl::set_error_handler;
/// let old_handler = set_error_handler(None);
/// ```
#[doc(alias = "gsl_set_error_handler")]
#[allow(static_mut_refs)]
pub fn set_error_handler(
    f: Option<fn(&str, &str, u32, Error)>,
) -> Option<fn(&str, &str, u32, Error)> {
    unsafe {
        let out = CALLBACK.take();
        match f {
            Some(f) => {
                CALLBACK = Some(f);
                sys::gsl_set_error_handler(Some(inner_error_handler));
            }
            None => {
                sys::gsl_set_error_handler(None);
            }
        }
        out
    }
}

/// This function turns off the error handler by defining an error handler which does nothing. This
/// will cause the program to continue after any error, so the return values from any library
/// routines must be checked. This is the recommended behavior for production programs. The previous
/// handler is returned (so that you can restore it later).
#[doc(alias = "gsl_set_error_handler_off")]
#[allow(static_mut_refs)]
pub fn set_error_handler_off() -> Option<fn(&str, &str, u32, crate::Error)> {
    unsafe {
        sys::gsl_set_error_handler_off();
        CALLBACK.take()
    }
}

extern "C" fn inner_error_handler(
    reason: *const c_char,
    file: *const c_char,
    line: c_int,
    gsl_errno: c_int,
) {
    unsafe {
        if let Some(ref call) = CALLBACK {
            let s = CStr::from_ptr(reason);
            let f = CStr::from_ptr(file);
            if let Err(e) = Error::handle(gsl_errno, ()) {
                // Do nothing on success.
                call(
                    s.to_str().unwrap_or("Unknown"),
                    f.to_str().unwrap_or("Unknown"),
                    line as _,
                    e,
                );
            }
        }
    }
}

#[cfg(test)]
#[test]
fn test_error_handler() {
    use crate::{bessel, Error};

    set_error_handler_off();
    match bessel::K0_e(1e3) {
        Err(Error::UnderFlow) => println!("K0(1e3) underflowed"),
        _ => panic!("unexpected"),
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum EigenSort {
    /// ascending order in numerical value
    ValAsc,
    /// descending order in numerical value
    ValDesc,
    /// ascending order in magnitude
    AbsAsc,
    /// descending order in magnitude
    AbsDesc,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_eigen_sort_t> for EigenSort {
    fn into(self) -> sys::gsl_eigen_sort_t {
        match self {
            Self::ValAsc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_ASC,
            Self::ValDesc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_DESC,
            Self::AbsAsc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_ASC,
            Self::AbsDesc => sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_DESC,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_eigen_sort_t> for EigenSort {
    fn from(v: sys::gsl_eigen_sort_t) -> EigenSort {
        match v {
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_ASC => Self::ValAsc,
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_VAL_DESC => Self::ValDesc,
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_ASC => Self::AbsAsc,
            sys::gsl_eigen_sort_t_GSL_EIGEN_SORT_ABS_DESC => Self::AbsDesc,
            _ => panic!("Unknown EigenSort value"),
        }
    }
}

/// This gives the sign in the formula:
///
/// ```text
/// h(f) = \sum x(t) exp(+/- 2 pi i f t)
/// ```
///
/// where - is the forward transform direction and + the inverse direction
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum FftDirection {
    Forward,
    Backward,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_fft_direction> for FftDirection {
    fn into(self) -> sys::gsl_fft_direction {
        match self {
            Self::Forward => sys::gsl_fft_direction_gsl_fft_forward,
            Self::Backward => sys::gsl_fft_direction_gsl_fft_backward,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_fft_direction> for FftDirection {
    fn from(v: sys::gsl_fft_direction) -> FftDirection {
        match v {
            sys::gsl_fft_direction_gsl_fft_forward => Self::Forward,
            sys::gsl_fft_direction_gsl_fft_backward => Self::Backward,
            _ => panic!("Unknown FftDirection value"),
        }
    }
}

/// Used by [`VegasParams`][crate::VegasParams].
///
/// This determines whether vegas will use importance sampling or
/// stratified sampling, or whether it can pick on its own.  In low
/// dimensions vegas uses strict stratified sampling (more precisely,
/// stratified sampling is chosen if there are fewer than 2 bins per
/// box).
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum VegasMode {
    /// Importance sampling: allocate more sample points where the
    /// integrand is larger.
    Importance,
    /// Exclusively use importance sampling without any stratification.
    ImportanceOnly,
    /// Stratified sampling: divides the integration region into
    /// sub-regions and sample each sub-region separately.
    Stratified,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<std::os::raw::c_int> for VegasMode {
    fn into(self) -> std::os::raw::c_int {
        match self {
            Self::Importance => sys::GSL_VEGAS_MODE_IMPORTANCE,
            Self::ImportanceOnly => sys::GSL_VEGAS_MODE_IMPORTANCE_ONLY,
            Self::Stratified => sys::GSL_VEGAS_MODE_STRATIFIED,
        }
    }
}

#[doc(hidden)]
impl From<std::os::raw::c_int> for VegasMode {
    fn from(v: std::os::raw::c_int) -> VegasMode {
        match v {
            sys::GSL_VEGAS_MODE_IMPORTANCE => Self::Importance,
            sys::GSL_VEGAS_MODE_IMPORTANCE_ONLY => Self::ImportanceOnly,
            sys::GSL_VEGAS_MODE_STRATIFIED => Self::Stratified,
            _ => panic!("Unknown VegasMode value"),
        }
    }
}

/// Possible return values for an hadjust() evolution method for ordinary differential equations
#[allow(clippy::upper_case_acronyms)]
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum ODEiv {
    /// step was increased
    Inc,
    /// step unchanged
    Nil,
    /// step decreased
    Dec,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<c_int> for ODEiv {
    fn into(self) -> c_int {
        match self {
            Self::Inc => sys::GSL_ODEIV_HADJ_INC,
            Self::Nil => sys::GSL_ODEIV_HADJ_NIL,
            Self::Dec => sys::GSL_ODEIV_HADJ_DEC,
        }
    }
}

#[doc(hidden)]
impl From<c_int> for ODEiv {
    fn from(v: c_int) -> ODEiv {
        match v {
            sys::GSL_ODEIV_HADJ_INC => Self::Inc,
            sys::GSL_ODEIV_HADJ_NIL => Self::Nil,
            sys::GSL_ODEIV_HADJ_DEC => Self::Dec,
            _ => panic!("Unknown ODEiv value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum WaveletDirection {
    Forward,
    Backward,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_wavelet_direction> for WaveletDirection {
    fn into(self) -> sys::gsl_wavelet_direction {
        match self {
            Self::Forward => sys::gsl_wavelet_direction_gsl_wavelet_forward,
            Self::Backward => sys::gsl_wavelet_direction_gsl_wavelet_backward,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_wavelet_direction> for WaveletDirection {
    fn from(v: sys::gsl_wavelet_direction) -> WaveletDirection {
        match v {
            sys::gsl_wavelet_direction_gsl_wavelet_forward => Self::Forward,
            sys::gsl_wavelet_direction_gsl_wavelet_backward => Self::Backward,
            _ => panic!("Unknown WaveletDirection value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum SfLegendreNorm {
    Schmidt,
    SphericalHarmonic,
    Full,
    None,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_sf_legendre_t> for SfLegendreNorm {
    fn into(self) -> sys::gsl_sf_legendre_t {
        match self {
            SfLegendreNorm::Schmidt => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SCHMIDT,
            SfLegendreNorm::SphericalHarmonic => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SPHARM,
            SfLegendreNorm::Full => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_FULL,
            SfLegendreNorm::None => sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_NONE,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_sf_legendre_t> for SfLegendreNorm {
    fn from(v: sys::gsl_sf_legendre_t) -> SfLegendreNorm {
        match v {
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SCHMIDT => SfLegendreNorm::Schmidt,
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_SPHARM => SfLegendreNorm::SphericalHarmonic,
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_FULL => SfLegendreNorm::Full,
            sys::gsl_sf_legendre_t_GSL_SF_LEGENDRE_NONE => SfLegendreNorm::None,
            _ => panic!("Unknown SfLegendreNorm value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum CblasTranspose {
    NoTranspose,
    Transpose,
    ConjugateTranspose,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_TRANSPOSE> for CblasTranspose {
    fn into(self) -> sys::CBLAS_TRANSPOSE {
        match self {
            Self::NoTranspose => sys::CBLAS_TRANSPOSE_CblasNoTrans,
            Self::Transpose => sys::CBLAS_TRANSPOSE_CblasTrans,
            Self::ConjugateTranspose => sys::CBLAS_TRANSPOSE_CblasConjTrans,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_TRANSPOSE> for CblasTranspose {
    fn from(v: sys::CBLAS_TRANSPOSE) -> CblasTranspose {
        match v {
            sys::CBLAS_TRANSPOSE_CblasNoTrans => Self::NoTranspose,
            sys::CBLAS_TRANSPOSE_CblasTrans => Self::Transpose,
            sys::CBLAS_TRANSPOSE_CblasConjTrans => Self::ConjugateTranspose,
            _ => panic!("Unknown CblasTranspose value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum CblasUplo {
    Upper,
    Lower,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_UPLO> for CblasUplo {
    fn into(self) -> sys::CBLAS_UPLO {
        match self {
            Self::Upper => sys::CBLAS_UPLO_CblasUpper,
            Self::Lower => sys::CBLAS_UPLO_CblasLower,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_UPLO> for CblasUplo {
    fn from(v: sys::CBLAS_UPLO) -> CblasUplo {
        match v {
            sys::CBLAS_UPLO_CblasUpper => Self::Upper,
            sys::CBLAS_UPLO_CblasLower => Self::Lower,
            _ => panic!("Unknown CblasUplo value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum CblasOrder {
    RowMajor,
    ColumnMajor,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_ORDER> for CblasOrder {
    fn into(self) -> sys::CBLAS_ORDER {
        match self {
            Self::RowMajor => sys::CBLAS_ORDER_CblasRowMajor,
            Self::ColumnMajor => sys::CBLAS_ORDER_CblasColMajor,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_ORDER> for CblasOrder {
    fn from(v: sys::CBLAS_ORDER) -> CblasOrder {
        match v {
            sys::CBLAS_ORDER_CblasRowMajor => Self::RowMajor,
            sys::CBLAS_ORDER_CblasColMajor => Self::ColumnMajor,
            _ => panic!("Unknown CblasOrder value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum CblasSide {
    Left,
    Right,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_SIDE> for CblasSide {
    fn into(self) -> sys::CBLAS_SIDE {
        match self {
            Self::Left => sys::CBLAS_SIDE_CblasLeft,
            Self::Right => sys::CBLAS_SIDE_CblasRight,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_SIDE> for CblasSide {
    fn from(v: sys::CBLAS_SIDE) -> CblasSide {
        match v {
            sys::CBLAS_SIDE_CblasLeft => Self::Left,
            sys::CBLAS_SIDE_CblasRight => Self::Right,
            _ => panic!("Unknown CblasSide value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum CblasDiag {
    NonUnit,
    Unit,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::CBLAS_SIDE> for CblasDiag {
    fn into(self) -> sys::CBLAS_SIDE {
        match self {
            Self::NonUnit => sys::CBLAS_DIAG_CblasNonUnit,
            Self::Unit => sys::CBLAS_DIAG_CblasUnit,
        }
    }
}

#[doc(hidden)]
impl From<sys::CBLAS_SIDE> for CblasDiag {
    fn from(v: sys::CBLAS_SIDE) -> CblasDiag {
        match v {
            sys::CBLAS_DIAG_CblasNonUnit => Self::NonUnit,
            sys::CBLAS_DIAG_CblasUnit => Self::Unit,
            _ => panic!("Unknown CblasDiag value"),
        }
    }
}

#[cfg(feature = "v2_5")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_5")))]
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum FilterEnd {
    PadZero,
    PadValue,
    Truncate,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
#[cfg(feature = "v2_5")]
impl Into<sys::gsl_filter_end_t> for FilterEnd {
    fn into(self) -> sys::gsl_filter_end_t {
        match self {
            Self::PadZero => sys::gsl_filter_end_t_GSL_FILTER_END_PADZERO,
            Self::PadValue => sys::gsl_filter_end_t_GSL_FILTER_END_PADVALUE,
            Self::Truncate => sys::gsl_filter_end_t_GSL_FILTER_END_TRUNCATE,
        }
    }
}

#[doc(hidden)]
#[cfg(feature = "v2_5")]
impl From<sys::gsl_filter_end_t> for FilterEnd {
    fn from(v: sys::gsl_filter_end_t) -> FilterEnd {
        match v {
            sys::gsl_filter_end_t_GSL_FILTER_END_PADZERO => Self::PadZero,
            sys::gsl_filter_end_t_GSL_FILTER_END_PADVALUE => Self::PadValue,
            sys::gsl_filter_end_t_GSL_FILTER_END_TRUNCATE => Self::Truncate,
            _ => panic!("Unknown FilterEnd value"),
        }
    }
}

#[cfg(feature = "v2_5")]
#[cfg_attr(feature = "dox", doc(cfg(feature = "v2_5")))]
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum FilterScale {
    MedianAbsoluteDeviation,
    InterQuartileRange,
    SN,
    QN,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
#[cfg(feature = "v2_5")]
impl Into<sys::gsl_filter_scale_t> for FilterScale {
    fn into(self) -> sys::gsl_filter_scale_t {
        match self {
            Self::MedianAbsoluteDeviation => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_MAD,
            Self::InterQuartileRange => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_IQR,
            Self::SN => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_SN,
            Self::QN => sys::gsl_filter_scale_t_GSL_FILTER_SCALE_QN,
        }
    }
}

#[doc(hidden)]
#[cfg(feature = "v2_5")]
impl From<sys::gsl_filter_scale_t> for FilterScale {
    fn from(v: sys::gsl_filter_scale_t) -> FilterScale {
        match v {
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_MAD => Self::MedianAbsoluteDeviation,
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_IQR => Self::InterQuartileRange,
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_SN => Self::SN,
            sys::gsl_filter_scale_t_GSL_FILTER_SCALE_QN => Self::QN,
            _ => panic!("Unknown FilterScale value"),
        }
    }
}
