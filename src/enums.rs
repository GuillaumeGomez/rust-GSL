//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

use std::os::raw::c_int;

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

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum Value {
    Success,
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

impl Value {
    pub fn is_success(self) -> bool {
        self == Self::Success
    }
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<c_int> for Value {
    fn into(self) -> c_int {
        match self {
            Self::Success => sys::GSL_SUCCESS,
            Self::Failure => sys::GSL_FAILURE,
            Self::Continue => sys::GSL_CONTINUE,
            Self::Domain => sys::GSL_EDOM,
            Self::Range => sys::GSL_ERANGE,
            Self::Fault => sys::GSL_EFAULT,
            Self::Invalid => sys::GSL_EINVAL,
            Self::Failed => sys::GSL_EFAILED,
            Self::Factorization => sys::GSL_EFACTOR,
            Self::Sanity => sys::GSL_ESANITY,
            Self::NoMemory => sys::GSL_ENOMEM,
            Self::BadFunction => sys::GSL_EBADFUNC,
            Self::RunAway => sys::GSL_ERUNAWAY,
            Self::MaxIteration => sys::GSL_EMAXITER,
            Self::ZeroDiv => sys::GSL_EZERODIV,
            Self::BadTolerance => sys::GSL_EBADTOL,
            Self::Tolerance => sys::GSL_ETOL,
            Self::UnderFlow => sys::GSL_EUNDRFLW,
            Self::OverFlow => sys::GSL_EOVRFLW,
            Self::Loss => sys::GSL_ELOSS,
            Self::Round => sys::GSL_EROUND,
            Self::BadLength => sys::GSL_EBADLEN,
            Self::NotSquare => sys::GSL_ENOTSQR,
            Self::Singularity => sys::GSL_ESING,
            Self::Diverge => sys::GSL_EDIVERGE,
            Self::Unsupported => sys::GSL_EUNSUP,
            Self::Unimplemented => sys::GSL_EUNIMPL,
            Self::Cache => sys::GSL_ECACHE,
            Self::Table => sys::GSL_ETABLE,
            Self::NoProgress => sys::GSL_ENOPROG,
            Self::NoProgressJacobian => sys::GSL_ENOPROGJ,
            Self::ToleranceF => sys::GSL_ETOLF,
            Self::ToleranceX => sys::GSL_ETOLX,
            Self::ToleranceG => sys::GSL_ETOLG,
            Self::EOF => sys::GSL_EOF,
            Self::Unknown(x) => x,
        }
    }
}

#[doc(hidden)]
impl From<c_int> for Value {
    fn from(v: c_int) -> Value {
        match v {
            sys::GSL_SUCCESS => Self::Success,
            sys::GSL_FAILURE => Self::Failure,
            sys::GSL_CONTINUE => Self::Continue,
            sys::GSL_EDOM => Self::Domain,
            sys::GSL_ERANGE => Self::Range,
            sys::GSL_EFAULT => Self::Fault,
            sys::GSL_EINVAL => Self::Invalid,
            sys::GSL_EFAILED => Self::Failed,
            sys::GSL_EFACTOR => Self::Factorization,
            sys::GSL_ESANITY => Self::Sanity,
            sys::GSL_ENOMEM => Self::NoMemory,
            sys::GSL_EBADFUNC => Self::BadFunction,
            sys::GSL_ERUNAWAY => Self::RunAway,
            sys::GSL_EMAXITER => Self::MaxIteration,
            sys::GSL_EZERODIV => Self::ZeroDiv,
            sys::GSL_EBADTOL => Self::BadTolerance,
            sys::GSL_ETOL => Self::Tolerance,
            sys::GSL_EUNDRFLW => Self::UnderFlow,
            sys::GSL_EOVRFLW => Self::OverFlow,
            sys::GSL_ELOSS => Self::Loss,
            sys::GSL_EROUND => Self::Round,
            sys::GSL_EBADLEN => Self::BadLength,
            sys::GSL_ENOTSQR => Self::NotSquare,
            sys::GSL_ESING => Self::Singularity,
            sys::GSL_EDIVERGE => Self::Diverge,
            sys::GSL_EUNSUP => Self::Unsupported,
            sys::GSL_EUNIMPL => Self::Unimplemented,
            sys::GSL_ECACHE => Self::Cache,
            sys::GSL_ETABLE => Self::Table,
            sys::GSL_ENOPROG => Self::NoProgress,
            sys::GSL_ENOPROGJ => Self::NoProgressJacobian,
            sys::GSL_ETOLF => Self::ToleranceF,
            sys::GSL_ETOLX => Self::ToleranceX,
            sys::GSL_ETOLG => Self::ToleranceG,
            sys::GSL_EOF => Self::EOF,
            x => Self::Unknown(x),
        }
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

/// The low-level integration rules in QUADPACK are identified by small integers (1-6). We'll use
/// symbolic constants to refer to them.
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum GaussKronrodRule {
    /// 15 point Gauss-Kronrod rule
    Gauss15,
    /// 21 point Gauss-Kronrod rule
    Gauss21,
    /// 31 point Gauss-Kronrod rule
    Gauss31,
    /// 41 point Gauss-Kronrod rule
    Gauss41,
    /// 51 point Gauss-Kronrod rule
    Gauss51,
    /// 61 point Gauss-Kronrod rule
    Gauss61,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<::std::os::raw::c_int> for GaussKronrodRule {
    fn into(self) -> ::std::os::raw::c_int {
        let x = match self {
            Self::Gauss15 => sys::GSL_INTEG_GAUSS15,
            Self::Gauss21 => sys::GSL_INTEG_GAUSS21,
            Self::Gauss31 => sys::GSL_INTEG_GAUSS31,
            Self::Gauss41 => sys::GSL_INTEG_GAUSS41,
            Self::Gauss51 => sys::GSL_INTEG_GAUSS51,
            Self::Gauss61 => sys::GSL_INTEG_GAUSS61,
        };
        x as _
    }
}

#[doc(hidden)]
impl From<::std::os::raw::c_int> for GaussKronrodRule {
    fn from(v: ::std::os::raw::c_int) -> GaussKronrodRule {
        match v as _ {
            sys::GSL_INTEG_GAUSS15 => Self::Gauss15,
            sys::GSL_INTEG_GAUSS21 => Self::Gauss21,
            sys::GSL_INTEG_GAUSS31 => Self::Gauss31,
            sys::GSL_INTEG_GAUSS41 => Self::Gauss41,
            sys::GSL_INTEG_GAUSS51 => Self::Gauss51,
            sys::GSL_INTEG_GAUSS61 => Self::Gauss61,
            _ => panic!("Unknown GaussKronrodRule value"),
        }
    }
}

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
/// Used by workspace for QAWO integrator
pub enum IntegrationQawo {
    Cosine,
    Sine,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_integration_qawo_enum> for IntegrationQawo {
    fn into(self) -> sys::gsl_integration_qawo_enum {
        match self {
            Self::Cosine => sys::gsl_integration_qawo_enum_GSL_INTEG_COSINE,
            Self::Sine => sys::gsl_integration_qawo_enum_GSL_INTEG_SINE,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_integration_qawo_enum> for IntegrationQawo {
    fn from(v: sys::gsl_integration_qawo_enum) -> IntegrationQawo {
        match v {
            sys::gsl_integration_qawo_enum_GSL_INTEG_COSINE => Self::Cosine,
            sys::gsl_integration_qawo_enum_GSL_INTEG_SINE => Self::Sine,
            _ => panic!("Unknown IntegrationQawo value"),
        }
    }
}

/// Used by VegasMonteCarlo struct
///
/// The possible choices are GSL_VEGAS_MODE_IMPORTANCE, GSL_VEGAS_MODE_
/// STRATIFIED, GSL_VEGAS_MODE_IMPORTANCE_ONLY. This determines whether vegas
/// will use importance sampling or stratified sampling, or whether it can pick on
/// its own. In low dimensions vegas uses strict stratified sampling (more precisely,
/// stratified sampling is chosen if there are fewer than 2 bins per box).
#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
pub enum VegasMode {
    Importance,
    ImportanceOnly,
    Stratified,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<::std::os::raw::c_int> for VegasMode {
    fn into(self) -> ::std::os::raw::c_int {
        match self {
            Self::Importance => sys::GSL_VEGAS_MODE_IMPORTANCE,
            Self::ImportanceOnly => sys::GSL_VEGAS_MODE_IMPORTANCE_ONLY,
            Self::Stratified => sys::GSL_VEGAS_MODE_STRATIFIED,
        }
    }
}

#[doc(hidden)]
impl From<::std::os::raw::c_int> for VegasMode {
    fn from(v: ::std::os::raw::c_int) -> VegasMode {
        match v {
            sys::GSL_VEGAS_MODE_IMPORTANCE => Self::Importance,
            sys::GSL_VEGAS_MODE_IMPORTANCE_ONLY => Self::ImportanceOnly,
            sys::GSL_VEGAS_MODE_STRATIFIED => Self::Stratified,
            _ => panic!("Unknown VegasMode value"),
        }
    }
}

/// Possible return values for an hadjust() evolution method for ordinary differential equations
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
