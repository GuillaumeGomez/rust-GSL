//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

//! Numerical Integration
/*!
## Introduction

Each algorithm computes an approximation to a definite integral of the form,

$$I = ∫_a^b f(x) w(x) dx$$

where $w(x)$ is a weight function (for general integrands $w(x)=1$). The
user provides absolute and relative error bounds (`epsabs`, `epsrel`)
which specify the following accuracy requirement,

$$|\RESULT - I| ≤ \max(\epsabs, \epsrel · |I|)$$

where RESULT is the numerical approximation obtained by the
algorithm.  The algorithms attempt to estimate the absolute error
$\ABSERR = |\RESULT - I|$ in such a way that the following
inequality holds,

$$|\RESULT - I| ≤ \ABSERR ≤ \max(\epsabs, \epsrel · |I|)$$

In short, the routines return the first approximation which has an
absolute error smaller than epsabs or a relative error smaller than
`epsrel`.

Note that this is an either-or constraint, not simultaneous. To
compute to a specified absolute error, set epsrel to zero. To compute
to a specified relative error, set epsabs to zero. The routines will
fail to converge if the error bounds are too stringent, but always
return the best approximation obtained up to that stage.

The algorithms in QUADPACK use a naming convention based on the
following letters,

- Q — quadrature routine

- N — non-adaptive integrator
- A — adaptive integrator

- G — general integrand (user-defined)
- W — weight function with integrand

- S — singularities can be more readily integrated
- P — points of special difficulty can be supplied
- I — infinite range of integration
- O — oscillatory weight function, cos or sin
- F — Fourier integral
- C — Cauchy principal value

The algorithms are built on pairs of quadrature rules, a higher order
rule and a lower order rule. The higher order rule is used to compute
the best approximation to an integral over a small range. The
difference between the results of the higher order rule and the lower
order rule gives an estimate of the error in the approximation.

 * [Integrands without weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-without-weight-functions.html#Integrands-without-weight-functions)
 * [Integrands with weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-with-weight-functions.html#Integrands-with-weight-functions)
 * [Integrands with singular weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-with-singular-weight-functions.html#Integrands-with-singular-weight-functions)

## QNG non-adaptive Gauss-Kronrod integration

The QNG algorithm is a non-adaptive procedure which uses fixed
Gauss-Kronrod-Patterson abscissae to sample the integrand at a maximum
of 87 points. It is provided for fast integration of smooth functions.
- [`qng`]

## QAG adaptive integration

The QAG algorithm is a simple adaptive integration procedure. The
integration region is divided into subintervals, and on each iteration
the subinterval with the largest estimated error is bisected. This
reduces the overall error rapidly, as the subintervals become
concentrated around local difficulties in the integrand. These
subintervals are managed by an [`IntegrationWorkspace`], which
handles the memory for the subinterval ranges, results and error
estimates.
- [`qag`]

### Adaptive integration with singularities

The presence of an integrable singularity in the integration region
causes an adaptive routine to concentrate new subintervals around the
singularity.  As the subintervals decrease in size the successive
approximations to the integral converge in a limiting fashion.  This
approach to the limit can be accelerated using an extrapolation
procedure.  The QAGS algorithm combines adaptive bisection with the
Wynn epsilon-algorithm to speed up the integration of many types of
integrable singularities.
- [`qags`]
- [`qagp`]: adaptive integration with known singular points.

### Adaptive integration on infinite intervals

- [`qagi`]: adaptive integration on $(-∞, ∞)$,
- [`qagiu`]: adaptive integration on $(a, ∞)$,
- [`qagil`]: adaptive integration on $(-∞, b)$.

### Adaptive integration for Cauchy principal values

- [`qawc`]

## Adaptive integration for singular weight

- [`qaws`]: weights of the form $(x-a)^α (b-x)^β \ln^μ (x-a) \ln^ν (b-x)$.

## Adaptive integration for oscillatory functions

- [`qawo`]: weights of the form $\cos(ωx)$ or $\sin(ωx)$,
- [`qawf`]: adaptive integration for Fourier integrals.

## Doubly-adaptive integration

- [`cquad`]

## Gauss-Legendre integration

The fixed-order Gauss-Legendre integration routines are provided for
fast integration of smooth functions with known polynomial order. The
$n$-point Gauss-Legendre rule is exact for polynomials of order $2n - 1$
or less.  For example, these rules are useful when integrating basis
functions to form mass matrices for the Galerkin method.  Unlike other
numerical integration routines within the library, these routines do
not accept absolute or relative error bounds.

## Gauss-Legendre integration

- [`GLFixedTable::glfixed`]

## Fixed point quadratures

The routines in this section approximate an integral by the sum

$$∫_a^b f(x) w(x) dx ≈ \sum_{i=1}^n w_i f(x_i)$$

where $f$ is the function to be integrated and $w$ is a weighting
function.  The $n$ weights $w_i$ and nodes $x_i$ are carefully chosen
so that the result is exact when $f$ is a polynomial of degree $2n -
1$ or less.  Once the user chooses the order $n$ and weighting
function $w$, the weights $w_i$ and nodes $x_i$ can be precomputed and
used to efficiently evaluate integrals for any number of functions $f$.

This method works best when $f$ is well approximated by a polynomial
on the interval $(a,b)$, and so is not suitable for functions with
singularities.  Since the user specifies ahead of time how many
quadrature nodes will be used, these routines do not accept absolute
or relative error bounds.  The table below lists the weighting
functions currently supported (see also [`IntegrationFixedType`]).

| Name             | Interval | Weighting function $w$    | Constraints     |
|------------------|----------|---------------------------|-----------------|
| Legendre         | $(a,b)$  | 1                         | $b > a$         |
| Chebyshev Type 1 | $(a,b)$  | $1/ \sqrt{(b-x)(x-a)}$    | $b > a$         |
| Gegenbauer       | $(a,b)$  | $((b-x) (x-a))^α$         | $α > -1, b > a$ |
| Jacobi           | $(a,b)$  | $(b-x)^α (x-a)^β$         | $α,β>-1, b > a$ |
| Laguerre         | $(a,∞)$  | $(x-a)^α \exp(-b(x-a))$   | $α > -1, b > 0$ |
| Hermite          | $(-∞,∞)$ | $\|x-a\|^α \exp(-b(x-a)^2)$ | $α > -1, b > 0$ |
| Exponential      | $(a,b)$  | $\|x - (a + b)/2\|^α$     | $α > -1, b > a$ |
| Rational         | $(a,∞)$  | $(x - a)^α (x + b)^β$     | $α > -1, α + β + 2n < 0, a + b > 0$ |
| Chebyshev Type 2 | $(a,b)$  | $\sqrt{(b - x) (x - a)}$  | $b > a$         |

- [`IntegrationFixedWorkspace::fixed`]


## References and Further Reading

The following book is the definitive reference for QUADPACK, and was
written by the original authors. It provides descriptions of the
algorithms, program listings, test programs and examples. It also
includes useful advice on numerical integration and many references to
the numerical integration literature used in developing QUADPACK.

- R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber,
  D.K. Kahaner. QUADPACK A subroutine package for automatic integration
  Springer Verlag, 1983.

The CQUAD integration algorithm is described in the following paper:

- P. Gonnet, “Increasing the Reliability of Adaptive Quadrature Using
  Explicit Interpolants”, ACM Transactions on Mathematical Software,
  Volume 37 (2010), Issue 3, Article 26.
*/

use crate::ffi::FFI;
use crate::Error;

/// Add to `self` methods `epsabs` and `epsrel` assuming the struct
/// has fields of the same name.
macro_rules! integ_builder {
    ($(#[$doc: meta])* $name: ident <$($l: lifetime,)*>,
        $($field: ident $ty: ty ,)*
    ) => {
        integ_builder!(with epsrel; $(#[$doc])* $name <$($l,)*>,
            $($field $ty,)*);
    };
    // qawf does not support `epsrel`.
    (with $($epsrel: ident)?;
        $(#[$doc: meta])* $name: ident <$($l: lifetime,)*>,
        $($field: ident $ty: ty ,)*
    ) => {
        $(#[$doc])*
        #[must_use]
        pub struct $name<$($l,)* F> {
            f: F,
            epsabs: f64,
            $($epsrel: f64,)?
            $($field: $ty,)*
        }

        impl<$($l,)* F> $name<$($l,)* F> {
            fn new(f: F, $($field: $ty,)*) -> Self {
                Self { f, epsabs: 1e-7, $($epsrel: 1e-7,)? $($field,)* }
            }

            /// Set the absolute error bound.
            pub fn epsabs(mut self, eps: f64) -> Self {
                if eps < 0. {
                    panic!("{}::epsabs: eps = {} should be >= 0",
                        stringify!($name), eps);
                }
                self.epsabs = eps;
                self
            }

            $(
            /// Set the relative error bound.
            pub fn $epsrel(mut self, eps: f64) -> Self {
                if eps < 0. {
                    panic!("{}::epsrel: eps = {} should be >= 0",
                        stringify!($name), eps);
                }
                self.epsrel = eps;
                self
            }
            )?
        }
    };
}

enum BorrowedOrOwned<'a, T> {
    Borrowed(&'a mut T),
    Owned(T),
}

impl<T> BorrowedOrOwned<'_, T> {
    fn to_mut(&mut self) -> &mut T {
        match self {
            Self::Borrowed(t) => t,
            Self::Owned(ref mut t) => t,
        }
    }
}

/// Create a builder with a workspace and limit on the number of
/// subintervals.  Optional lifetimes may be needed if the constructor
/// depends on borrowed ressources.
macro_rules! integ_builder_workspace {
    ($(#[$doc: meta])* $name: ident <$($l: lifetime,)*>,
        $($field: ident $ty: ty ,)*
    ) => {
        integ_builder!{
            $(#[$doc])* $name<'w, $($l,)*>,
            workspace Option<BorrowedOrOwned<'w, IntegrationWorkspace>>,
            limit usize,
            $($field $ty,)*
        }
        integ_builder_workspace!(methods $name <$($l,)*>, $($field $ty,)*);
    };
    (no-epsrel; $(#[$doc: meta])* $name: ident <$($l: lifetime,)*>,
        $($field: ident $ty: ty ,)*
    ) => {
        integ_builder!{with; // No epsrel
            $(#[$doc])* $name<'w, $($l,)*>,
            workspace Option<BorrowedOrOwned<'w, IntegrationWorkspace>>,
            limit usize,
            $($field $ty,)*
        }
        integ_builder_workspace!(methods $name <$($l,)*>, $($field $ty,)*);
    };
    (methods $name: ident <$($l: lifetime,)*>,
        $($field: ident $ty: ty ,)*
    ) => {
        impl<$($l,)* F> $name<'static, $($l,)* F> {
            fn with_workspace(f: F, $($field: $ty,)*) -> Self {
                Self::new(f, None, 1000, $($field,)*)
            }
        }

        impl<'w, $($l,)* F> $name<'w, $($l,)* F> {
            /// Use the workspace `w` to compute the integral.
            pub fn workspace(
                self,
                w: &mut IntegrationWorkspace
            ) -> $name<'_, $($l,)* F> {
                $name {
                    workspace: Some(BorrowedOrOwned::Borrowed(w)),
                    .. self
                }
            }

            /// Set the maximum number of subintervals in the adaptive
            /// procedure.
            pub fn limit(mut self, limit: usize) -> Self {
                self.limit = limit;
                self
            }

            /// Add an owned workspace if none was present.  This is
            /// to call before launching the computation to make sure
            /// `limit` has its final value.
            fn add_workspace(&mut self) -> Result<(), Error> {
                if self.workspace.is_none() {
                    let w = IntegrationWorkspace::new(self.limit)
                    .ok_or(Error::NoMemory)?;
                    self.workspace = Some(BorrowedOrOwned::Owned(w));
                }
                Ok(())
            }
        }
    }
}

/// Get the workspace, allocating one if necessary.  This has to be a
/// macro to use the projection on the struct field — not mutably
/// borrowing the entire `$s`.
macro_rules! get_workspace {
    ($s: ident) => {{
        $s.add_workspace()?;
        $s.workspace.as_mut().unwrap().to_mut()
    }};
}

/// Non-adaptive Gauss-Kronrod integration.
///
/// This function applies the Gauss-Kronrod 10-point, 21-point,
/// 43-point and 87-point integration rules in succession until an
/// estimate of the integral of `f` over the interval $(a,b)$ is
/// achieved within the desired absolute and relative error limits,
/// `eps_abs` and `eps_rel`.  The function returns the final
/// approximation, result, an estimate of the absolute error, `abserr`
/// and the number of function evaluations used, `neval`.  The
/// Gauss-Kronrod rules are designed in such a way that each rule uses
/// all the results of its predecessors, in order to minimize the
/// total number of function evaluations.
#[doc(alias = "gsl_integration_qng")]
pub fn qng<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> Qng<F> {
    Qng::new(f, a, b)
}

integ_builder! {
    /// Builder to compute integrals using [`qng`].
    Qng<>,
    a f64, b f64,
}

impl<F: Fn(f64) -> f64> Qng<F> {
    /// Return `(result, abs_err, n_eval)`.
    pub fn val_err_n(&self) -> Result<(f64, f64, usize), Error> {
        let function = wrap_callback!(self.f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut n_eval = 0;

        let ret = unsafe {
            sys::gsl_integration_qng(
                &function,
                self.a,
                self.b,
                self.epsabs,
                self.epsrel,
                &mut result,
                &mut abs_err,
                &mut n_eval,
            )
        };
        Error::handle(ret, (result, abs_err, n_eval))
    }
}

/// Low-level integration rules in QUADPACK.
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
impl Into<std::os::raw::c_int> for GaussKronrodRule {
    fn into(self) -> std::os::raw::c_int {
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
impl From<std::os::raw::c_int> for GaussKronrodRule {
    fn from(v: std::os::raw::c_int) -> GaussKronrodRule {
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

ffi_wrapper!(
    IntegrationWorkspace,
    *mut sys::gsl_integration_workspace,
    gsl_integration_workspace_free,
    "Manage subintervals for adaptive integration.  It handles the memory for the subinterval ranges, results and error estimates."
);

impl IntegrationWorkspace {
    /// This function allocates a workspace sufficient to hold `n`
    /// double precision intervals, their integration results and
    /// error estimates. One workspace may be used multiple times as
    /// all necessary reinitialization is performed automatically by
    /// the integration routines.
    #[doc(alias = "gsl_integration_workspace_alloc")]
    pub fn new(n: usize) -> Option<IntegrationWorkspace> {
        let tmp = unsafe { sys::gsl_integration_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// Return the maximum double precision intervals that the
    /// workspace can handle.
    pub fn limit(&self) -> usize {
        unsafe { (*self.unwrap_shared()).limit }
    }
    pub fn size(&self) -> usize {
        unsafe { (*self.unwrap_shared()).size }
    }
    pub fn nrmax(&self) -> usize {
        unsafe { (*self.unwrap_shared()).nrmax }
    }
    pub fn i(&self) -> usize {
        unsafe { (*self.unwrap_shared()).i }
    }
    pub fn maximum_level(&self) -> usize {
        unsafe { (*self.unwrap_shared()).maximum_level }
    }
}

/// Simple adaptive integration procedure.
///
/// This function applies an integration rule adaptively until an
/// estimate of the integral of `f` over $(a,b)$ is achieved within
/// the desired absolute and relative error limits, `epsabs` and
/// `epsrel`.  The function returns the final approximation, `result`,
/// and an estimate of the absolute error, `abserr` (see
/// [`Qag::val_err`]).  The integration rule is determined by the
/// value of `key` corresponding to the 15, 21, 31, 41, 51 and 61
/// point Gauss-Kronrod rules (set with [`Qag::rule`], the default
/// being the 21 point rule).  The higher-order rules give better
/// accuracy for smooth functions, while lower-order rules save time
/// when the function contains local difficulties, such as
/// discontinuities.
///
/// On each iteration the adaptive integration strategy bisects the
/// interval with the largest error estimate.  The subintervals and
/// their results are stored in the memory provided by a workspace
/// (you can pass one with [`Qag::workspace`] or one will be allocated
/// automatically).  The maximum number of subintervals is set by
/// [`Qag::limit`], which may not exceed the allocated size of the workspace.
///
/// # Example
///
/// ```
/// use rgsl::integration::qag;
/// let (v, errabs) = qag(|x| x, 0., 1.).val_err()?;
/// assert_eq!(v, 0.5);
/// assert!(errabs < 1e-14);
/// # Ok::<(), rgsl::Error>(())
/// ```
#[doc(alias = "gsl_integration_qag")]
pub fn qag<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> Qag<'static, F> {
    Qag::with_workspace(f, a, b, GaussKronrodRule::Gauss21)
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qag`].
    Qag<>,
    a f64, b f64, key GaussKronrodRule,
}

impl<F: Fn(f64) -> f64> Qag<'_, F> {
    pub fn rule(mut self, key: GaussKronrodRule) -> Self {
        self.key = key;
        self
    }

    /// Return $(∫_a^b f, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let function = wrap_callback!(self.f, F);
        let ret = unsafe {
            sys::gsl_integration_qag(
                &function,
                self.a,
                self.b,
                self.epsabs,
                self.epsrel,
                self.limit,
                self.key.into(),
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

/// Adaptive integration with singularities.
///
/// This function applies the Gauss-Kronrod 21-point integration rule
/// adaptively until an estimate of the integral of `f` over the
/// interval $(a,b)$ is achieved within the desired absolute and
/// relative error limits (set with [`Qags::epsabs`] and
/// [`Qags::epsrel`]).  The results are extrapolated using the
/// epsilon-algorithm, which accelerates the convergence of the
/// integral in the presence of discontinuities and integrable
/// singularities.  [`Qags::val_err`] returns the final approximation
/// from the extrapolation and an estimate of the absolute error.  The
/// subintervals and their results are stored in the memory provided
/// by workspace. The maximum number of subintervals is set by
/// [`Qags::limit`], which may not exceed the allocated size of the
/// workspace.
#[doc(alias = "gsl_integration_qags")]
pub fn qags<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> Qags<'static, F> {
    Qags::with_workspace(f, a, b)
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qags`].
    Qags<>,
    a f64, b f64,
}

impl<F: Fn(f64) -> f64> Qags<'_, F> {
    /// Return $(∫_a^b f, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let function = wrap_callback!(self.f, F);
        let ret = unsafe {
            sys::gsl_integration_qags(
                &function,
                self.a,
                self.b,
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

/// Adaptive integration with known singular points.
///
/// This function applies the adaptive integration algorithm QAGS
/// taking account of the user-supplied locations of singular
/// points. The array pts of length npts should contain the endpoints
/// of the integration ranges defined by the integration region and
/// locations of the singularities.
///
/// For example, to integrate over the interval $(a,b)$ with break-points
/// at `x1`, `x2`, `x3` (where `a < x1 < x2 < x3 < b`) the following pts
/// array should be used
///
/// ```text
/// pts[0] = a
/// pts[1] = x1
/// pts[2] = x2
/// pts[3] = x3
/// pts[4] = b
/// ```
///
/// If you know the locations of the singular points in the
/// integration region then this routine will be faster than QAGS.
#[doc(alias = "gsl_integration_qagp")]
pub fn qagp<F: Fn(f64) -> f64>(f: F, pts: &mut [f64]) -> Qagp<'static, '_, F> {
    Qagp::with_workspace(f, pts)
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qagp`].
    Qagp<'pts,>,
    pts &'pts mut [f64],
}

impl<F: Fn(f64) -> f64> Qagp<'_, '_, F> {
    /// Return $(∫_a^b f, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let function = wrap_callback!(self.f, F);

        let ret = unsafe {
            sys::gsl_integration_qagp(
                &function,
                self.pts.as_mut_ptr(),
                self.pts.len() as _,
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

/// Adaptive integration on $(-∞, ∞)$.
///
/// This function computes the integral of the function `f` over the
/// infinite interval $(-∞,+∞)$.  The integral is mapped
/// onto the semi-open interval $(0,1]$ using the transformation $ x =
/// (1-t)/t$,
///
/// $$∫_{-∞}^{+∞} f(x) dx
/// = ∫_0^1 \bigl( f((1-t)/t) + f((-1+t)/t) \bigr)/t^2 dt.$$
///
/// It is then integrated using the QAGS algorithm. The normal
/// 21-point Gauss-Kronrod rule of QAGS is replaced by a 15-point
/// rule, because the transformation can generate an integrable
/// singularity at the origin.  In this case a lower-order rule is
/// more efficient.
#[doc(alias = "gsl_integration_qagi")]
pub fn qagi<F: Fn(f64) -> f64>(f: F) -> Qagi<'static, F> {
    Qagi::with_workspace(f)
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qagi`].
    Qagi<>,
}

impl<F: Fn(f64) -> f64> Qagi<'_, F> {
    /// Return $(∫_{-∞}^∞ f, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let mut function = wrap_callback!(self.f, F);

        let ret = unsafe {
            sys::gsl_integration_qagi(
                &mut function,
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

/// Adaptive integration on $(a, ∞)$.
///
/// This function computes the integral of the function f over the
/// semi-infinite interval $(a, +∞)$.  The integral is mapped
/// onto the semi-open interval $(0,1]$ using the transformation:
/// $x = a + (1-t)/t$,
///
/// $$∫_a^{+∞} f(x) dx = ∫_0^1 f(a + (1-t)/t)/t^2 dt$$
///
/// and then integrated using the QAGS algorithm.
#[doc(alias = "gsl_integration_qagiu")]
pub fn qagiu<F: Fn(f64) -> f64>(f: F, a: f64) -> Qagiu<'static, F> {
    Qagiu::with_workspace(f, a)
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qagiu`].
    Qagiu<>,
    a f64,
}

impl<F: Fn(f64) -> f64> Qagiu<'_, F> {
    /// Return $(∫_a^∞ f, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let mut function = wrap_callback!(self.f, F);

        let ret = unsafe {
            sys::gsl_integration_qagiu(
                &mut function,
                self.a,
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

/// Adaptive integration on $(-∞, b)$.
///
/// This function computes the integral of the function `f` over the
/// semi-infinite interval $(-∞,b)$.  The integral is mapped onto
/// the semi-open interval $(0,1]$ using the transformation:
/// $x = b - (1-t)/t$,
///
/// $$∫_{-∞}^{b} f(x) dx = ∫_0^1 f(b - (1-t)/t)/t^2 dt$$
///
/// and then integrated using the QAGS algorithm.
#[doc(alias = "gsl_integration_qagil")]
pub fn qagil<F: Fn(f64) -> f64>(f: F, b: f64) -> Qagil<'static, F> {
    Qagil::with_workspace(f, b)
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qagil`].
    Qagil<>,
    b f64,
}

impl<F: Fn(f64) -> f64> Qagil<'_, F> {
    /// Return $(∫_{-∞}^b f, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let mut function = wrap_callback!(self.f, F);

        let ret = unsafe {
            sys::gsl_integration_qagil(
                &mut function,
                self.b,
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

/// Adaptive integration for Cauchy principal values.
///
/// This function computes the Cauchy principal value of the integral
/// of f over the interval $(a,b)$, with a singularity at `c`,
///
/// $$I = ∫_a^b \frac{f(x)}{x - c} dx
///     = \lim _{ε→0} \biggl( ∫_a^{c-ε} \frac{f(x)}{x - c} dx
///                         + ∫ _{c+ε}^b \frac{f(x)}{x - c} dx \biggr).$$
///
/// The adaptive bisection algorithm of QAG is used, with
/// modifications to ensure that subdivisions do not occur at the
/// singular point $x = c$.
///
/// When a subinterval contains the point $x = c$ or is close to it
/// then a special 25-point modified Clenshaw-Curtis rule is used to
/// control the singularity. Further away from the singularity the
/// algorithm uses an ordinary 15-point Gauss-Kronrod integration
/// rule.
#[doc(alias = "gsl_integration_qawc")]
pub fn qawc<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, c: f64) -> Qawc<'static, F> {
    Qawc::with_workspace(f, a, b, c)
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qawc`].
    Qawc<>,
    a f64, b f64, c f64,
}

impl<F: Fn(f64) -> f64> Qawc<'_, F> {
    /// Return $(∫_a^b f, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let mut function = wrap_callback!(self.f, F);

        let ret = unsafe {
            sys::gsl_integration_qawc(
                &mut function,
                self.a,
                self.b,
                self.c,
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

ffi_wrapper!(
    QawsTable,
    *mut sys::gsl_integration_qaws_table,
    gsl_integration_qaws_table_free,
    "Table of Chebyshev moments for [`qaws`]."
);

impl QawsTable {
    /// This function allocates space for a gsl_integration_qaws_table
    /// struct describing a singular weight function W(x) with the
    /// parameters α = `alpha`, β = `beta`, μ = `mu` and ν = `nu`,
    ///
    /// $$w(x) = (x-a)^α (b-x)^β \ln^μ (x-a) \ln^ν (b-x)$$
    ///
    /// where $α > -1$, $β > -1$, and $μ ∈ \{0, 1\}$, $ν ∈ \{0, 1\}$.
    /// The weight function can take four different forms depending on
    /// the values of μ and ν,
    ///
    /// - $w(x) = (x-a)^α (b-x)^β$                         ($μ = 0, ν = 0$)
    /// - $w(x) = (x-a)^α (b-x)^β \ln(x-a)$                ($μ = 1, ν = 0$)
    /// - $w(x) = (x-a)^α (b-x)^β \ln(b-x)$                ($μ = 0, ν = 1$)
    /// - $w(x) = (x-a)^α (b-x)^β \ln(x-a) \ln(b-x)$       ($μ = 1, ν = 1$)
    ///
    /// The singular points $(a,b)$ do not have to be specified until
    /// the integral is computed, where they are the endpoints of the
    /// integration range.
    #[doc(alias = "gsl_integration_qaws_table_alloc")]
    pub fn new(alpha: f64, beta: f64, mu: i32, nu: i32) -> Option<QawsTable> {
        let tmp = unsafe { sys::gsl_integration_qaws_table_alloc(alpha, beta, mu, nu) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function modifies the parameters (`alpha`, `beta`, `mu`, `nu`).
    #[doc(alias = "gsl_integration_qaws_table_set")]
    pub fn set(&mut self, alpha: f64, beta: f64, mu: i32, nu: i32) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_integration_qaws_table_set(self.unwrap_unique(), alpha, beta, mu, nu)
        };
        Error::handle(ret, ())
    }
}

/// Adaptive integration with a singular weight.
///
/// This function computes the integral of the function `f` over the
/// interval $(a,b)$ with the singular weight function $(x-a)^α
/// (b-x)^β \ln^μ (x-a) \ln^ν (b-x)$.  The parameters of the weight
/// function ($α, β, μ, ν$) are set by the calling method
/// [`QawsWeight::table`] or [`QawsWeight::w`] on the return value of
/// this function.  The integral is,
///
/// $$I = ∫_a^b f(x) (x-a)^α (b-x)^β \ln^μ (x-a) \ln^ν (b-x) dx.$$
///
/// The adaptive bisection algorithm of QAG is used.  When a
/// subinterval contains one of the endpoints then a special 25-point
/// modified Clenshaw-Curtis rule is used to control the
/// singularities.  For subintervals which do not include the
/// endpoints an ordinary 15-point Gauss-Kronrod integration rule is used.
#[doc(alias = "gsl_integration_qaws")]
pub fn qaws<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> QawsWeight<F> {
    QawsWeight { f, a, b }
}

/// Set weight for a [`qaws`] integration rule.
pub struct QawsWeight<F> {
    f: F,
    a: f64,
    b: f64,
}

impl<F> QawsWeight<F> {
    pub fn table(self, t: &mut QawsTable) -> Qaws<'static, '_, F> {
        Qaws::with_workspace(self.f, self.a, self.b, BorrowedOrOwned::Borrowed(t))
    }

    pub fn w(self, alpha: f64, beta: f64, mu: i32, nu: i32) -> Qaws<'static, 'static, F> {
        let table = QawsTable::new(alpha, beta, mu, nu)
            // FIXME: delay until eval?  Can then return Error::NoMemory
            .expect("rgsl::integration::QawsWeight: cannot allocate");
        Qaws::with_workspace(self.f, self.a, self.b, BorrowedOrOwned::Owned(table))
    }
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qaws`].
    Qaws<'t,>,
    a f64, b f64, table BorrowedOrOwned<'t, QawsTable>,
}

impl<F: Fn(f64) -> f64> Qaws<'_, '_, F> {
    /// Return $(∫_a^b f(x) w(x) dx, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;
        let w = get_workspace!(self);
        let mut function = wrap_callback!(self.f, F);

        let ret = unsafe {
            sys::gsl_integration_qaws(
                &mut function,
                self.a,
                self.b,
                self.table.to_mut().unwrap_unique(),
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

ffi_wrapper!(
    QawoTable,
    *mut sys::gsl_integration_qawo_table,
    gsl_integration_qawo_table_free,
    "Table of Chebyshev moments according to parameters $ω$ and $L$ for [`qawo`] and [`qawf`]."
);

#[derive(Clone, PartialEq, PartialOrd, Debug, Copy)]
/// Choice of the type of oscillatory function.
pub enum QawOsc {
    Cosine,
    Sine,
}

#[doc(hidden)]
#[allow(clippy::from_over_into)]
impl Into<sys::gsl_integration_qawo_enum> for QawOsc {
    fn into(self) -> sys::gsl_integration_qawo_enum {
        match self {
            Self::Cosine => sys::gsl_integration_qawo_enum_GSL_INTEG_COSINE,
            Self::Sine => sys::gsl_integration_qawo_enum_GSL_INTEG_SINE,
        }
    }
}

#[doc(hidden)]
impl From<sys::gsl_integration_qawo_enum> for QawOsc {
    fn from(v: sys::gsl_integration_qawo_enum) -> QawOsc {
        match v {
            sys::gsl_integration_qawo_enum_GSL_INTEG_COSINE => Self::Cosine,
            sys::gsl_integration_qawo_enum_GSL_INTEG_SINE => Self::Sine,
            _ => panic!("Unknown QawOsc value"),
        }
    }
}

impl QawoTable {
    /// Default number of levels (see [`new`]).
    const N: usize = 50;

    /// This function allocates space for a [`QawoTable¯] and its
    /// associated workspace describing a sine or cosine weight
    /// function $w$ with the parameters $(ω, L)$,
    ///
    /// $w(x) = \sin(ω x)$,
    /// $w(x) = \cos(ω x)$.
    ///
    /// The parameter $L$ must be the length of the interval over
    /// which the function will be integrated $L = b - a$.  The choice
    /// of sine or cosine is made with the parameter `sine` which should
    /// be chosen from one of the values of [`QawOsc`].
    ///
    /// The [`QawoTable`] is a table of the trigonometric coefficients
    /// required in the integration process.  The parameter `n`
    /// determines the number of levels of coefficients that are
    /// computed.  Each level corresponds to one bisection of the
    /// interval $L$, so that `n` levels are sufficient for subintervals
    /// down to the length $L/2^n$.  The integration routine
    /// [`qawo`] returns the error [`Error::Table`] if the
    /// number of levels is insufficient for the requested accuracy.
    #[doc(alias = "gsl_integration_qawo_table_alloc")]
    pub fn new(omega: f64, l: f64, sine: QawOsc, n: usize) -> Option<QawoTable> {
        let tmp = unsafe { sys::gsl_integration_qawo_table_alloc(omega, l, sine.into(), n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// This function changes the parameters $ω, L$ and sine of the
    /// existing self workspace.
    #[doc(alias = "gsl_integration_qawo_table_set")]
    pub fn set(&mut self, omega: f64, l: f64, sine: QawOsc) -> Result<(), Error> {
        let ret = unsafe {
            sys::gsl_integration_qawo_table_set(self.unwrap_unique(), omega, l, sine.into())
        };
        Error::handle(ret, ())
    }

    /// This function allows the length parameter l of the self
    /// workspace to be changed.
    #[doc(alias = "gsl_integration_qawo_table_set_length")]
    pub fn set_length(&mut self, l: f64) -> Result<(), Error> {
        let ret = unsafe { sys::gsl_integration_qawo_table_set_length(self.unwrap_unique(), l) };
        Error::handle(ret, ())
    }
}

/// Adaptive integration for oscillatory functions (sin, cos).
///
/// This function uses an adaptive algorithm to compute the integral
/// of `f` over the interval $(a,b)$ with the weight function
/// $\sin(ωx)$ or $\cos(ωx)$ defined by the table `wf`,
///
/// $$I = ∫_a^b f(x) \sin(ω x) dx$$
/// $$I = ∫_a^b f(x) \cos(ω x) dx$$
///
/// The results are extrapolated using the epsilon-algorithm to
/// accelerate the convergence of the integral. The function returns
/// the final approximation from the extrapolation, result, and an
/// estimate of the absolute error, abserr. The subintervals and their
/// results are stored in the memory provided by workspace. The
/// maximum number of subintervals is given by limit, which may not
/// exceed the allocated size of the workspace.
///
/// Those subintervals with “large” widths d where $dω > 4$ are
/// computed using a 25-point Clenshaw-Curtis integration rule, which
/// handles the oscillatory behavior. Subintervals with a “small”
/// widths where $dω < 4$ are computed using a 15-point Gauss-Kronrod
/// integration.
#[doc(alias = "gsl_integration_qawo")]
pub fn qawo<F: Fn(f64) -> f64>(f: F, a: f64) -> QawoWeight<F> {
    QawoWeight { f, a }
}

/// Choice of the weight $\sin(ωx)$ or $\cos(ωx)$ for [`qawo`].
pub struct QawoWeight<F> {
    f: F,
    a: f64,
}

impl<F> QawoWeight<F> {
    pub fn table(self, wf: &mut QawoTable) -> Qawo<'static, '_, F> {
        Qawo::with_workspace(
            self.f,
            self.a,
            Some(BorrowedOrOwned::Borrowed(wf)),
            f64::NAN,
            f64::NAN,
            QawOsc::Cosine,
        ) // Last 3 do not matter
    }

    /// Set the weight to $\sin(ωx)$ on intervals of length `L`.
    pub fn sin(self, ω: f64, L: f64) -> Qawo<'static, 'static, F> {
        // Keep the info to create the table when `limit` is known.
        Qawo::with_workspace(self.f, self.a, None, ω, L, QawOsc::Sine)
    }

    /// Set the weight to $\cos(ωx)$ on intervals of length `L`.
    pub fn cos(self, ω: f64, L: f64) -> Qawo<'static, 'static, F> {
        Qawo::with_workspace(self.f, self.a, None, ω, L, QawOsc::Cosine)
    }
}

integ_builder_workspace! {
    /// Builder to compute integrals using [`qawo`].
    Qawo<'t,>,
    a f64, table Option<BorrowedOrOwned<'t, QawoTable>>,
    omega f64, l f64, sine QawOsc,
}

impl<F: Fn(f64) -> f64> Qawo<'_, '_, F> {
    /// Return $(∫_a^b f(x) w(x) dx, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        if self.table.is_none() {
            // If no table was provided, `'t` is `'static` which is
            // what we want for an owned table.
            let wf = QawoTable::new(self.omega, self.l, self.sine, QawoTable::N)
                .ok_or(Error::NoMemory)?;
            self.table = Some(BorrowedOrOwned::Owned(wf));
        }
        let w = get_workspace!(self);
        let mut function = wrap_callback!(self.f, F);
        let mut result = 0.;
        let mut abs_err = 0.;

        let ret = unsafe {
            sys::gsl_integration_qawo(
                &mut function,
                self.a,
                self.epsabs,
                self.epsrel,
                self.limit,
                w.unwrap_unique(),
                self.table.as_mut().unwrap().to_mut().unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

/// Adaptive integration for Fourier integrals.
///
/// This function attempts to compute a Fourier integral of the
/// function f over the semi-infinite interval $[a,+∞)$.
///
/// $$I = ∫_a^{+∞} f(x) \sin(ωx) dx$$
/// $$I = ∫_a^{+∞} f(x) \cos(ωx) dx$$
///
/// The parameter $ω$ and choice of $\sin$ or $\cos$ is taken from the
/// table set with the method [`QawfWeight::table`] or set with
/// [`QawfWeight::sin`] or [`QawfWeight::cos`] on the return value of
/// this function.  The integral is computed using the QAWO algorithm
/// over each of the subintervals,
///
/// $C_1 = [a, a + c]$
/// $C_2 = [a + c, a + 2 c]$
/// $... = ...$
/// $C_k = [a + (k-1) c, a + k c]$
///
/// where $c = (2 ⌊|ω|⌋ + 1) \pi/|ω|$.  The width $c$ is chosen to
/// cover an odd number of periods so that the contributions from the
/// intervals alternate in sign and are monotonically decreasing when
/// `f` is positive and monotonically decreasing.  The sum of this
/// sequence of contributions is accelerated using the
/// epsilon-algorithm.
///
/// This function works to an overall absolute tolerance of `abserr`
/// (the 2nd component returned by [`Qawf::val_err`]).  The following
/// strategy is used: on each interval $C_k$ the algorithm tries to
/// achieve the tolerance
///
/// $\TOL_k = u_k \abserr$
///
/// where $u_k = (1 - p)p^{k-1}$ and $p = 9/10$.  The sum of the
/// geometric series of contributions from each interval gives an
/// overall tolerance of abserr.
///
/// If the integration of a subinterval leads to difficulties then the
/// accuracy requirement for subsequent intervals is relaxed,
///
/// $$\TOL_k = u_k \max(\abserr, \max_{i<k}{E_i})$$
///
/// where $E_k$ is the estimated error on the interval $C_k$.
///
/// The subintervals and their results are stored in the memory
/// provided by workspace. The maximum number of subintervals is given
/// by limit, which may not exceed the allocated size of the
/// workspace. The integration over each subinterval uses the memory
/// provided by cycle_workspace as workspace for the QAWO algorithm.
#[doc(alias = "gsl_integration_qawf")]
pub fn qawf<F: Fn(f64) -> f64>(f: F, a: f64) -> QawfWeight<F> {
    QawfWeight { f, a }
}

/// Choice of the weight $\sin(ωx)$ or $\cos(ωx)$ for [`qawf`].
pub struct QawfWeight<F> {
    f: F,
    a: f64,
}

impl<F> QawfWeight<F> {
    /// Set the sine or cosine and the parameter $ω$ through the table
    /// `wf` (the value of $L$ does not matter, it is overridden to a
    /// value appropriate for Fourier integration).
    pub fn table(self, wf: &mut QawoTable) -> Qawf<'static, 'static, '_, F> {
        Qawf::with_workspace(
            self.f,
            self.a,
            Some(BorrowedOrOwned::Borrowed(wf)),
            None,
            QawOsc::Cosine,
            f64::NAN,
        )
        // Last 2 param do not matter
    }

    /// Set the weight to $\sin(ωx)$.
    pub fn sin(self, ω: f64) -> Qawf<'static, 'static, 'static, F> {
        Qawf::with_workspace(self.f, self.a, None, None, QawOsc::Sine, ω)
    }

    /// Set the weight to $\cos(ωx)$.
    pub fn cos(self, ω: f64) -> Qawf<'static, 'static, 'static, F> {
        Qawf::with_workspace(self.f, self.a, None, None, QawOsc::Cosine, ω)
    }
}

integ_builder_workspace! {no-epsrel;
    /// Builder to compute integrals using [`qawf`].
    Qawf<'cw, 't,>,
    a f64,
    table Option<BorrowedOrOwned<'t, QawoTable>>,
    cycle_workspace Option<BorrowedOrOwned<'cw, IntegrationWorkspace>>,
    sine QawOsc,
    omega f64,
}

impl<'w, 't, F: Fn(f64) -> f64> Qawf<'w, '_, 't, F> {
    pub fn cycle_workspace<'a>(self, w: &'a mut IntegrationWorkspace) -> Qawf<'w, 'a, 't, F> {
        Qawf {
            cycle_workspace: Some(BorrowedOrOwned::Borrowed(w)),
            ..self
        }
    }

    /// Return $(∫_a^b f(x) w(x) dx, \abserr)$.
    pub fn val_err(&mut self) -> Result<(f64, f64), Error> {
        let mut result = 0.;
        let mut abs_err = 0.;

        if self.table.is_none() {
            let wf =
                QawoTable::new(self.omega, 1., self.sine, QawoTable::N).ok_or(Error::NoMemory)?;
            self.table = Some(BorrowedOrOwned::Owned(wf));
        }
        if self.cycle_workspace.is_none() {
            let cw = IntegrationWorkspace::new(100).ok_or(Error::NoMemory)?;
            self.cycle_workspace = Some(BorrowedOrOwned::Owned(cw));
        }
        let w = get_workspace!(self);
        let wf = self.table.as_mut().unwrap().to_mut();
        let cw = self.cycle_workspace.as_mut().unwrap().to_mut();
        let mut function = wrap_callback!(self.f, F);
        let ret = unsafe {
            sys::gsl_integration_qawf(
                &mut function,
                self.a,
                self.epsabs,
                self.limit,
                w.unwrap_unique(),
                cw.unwrap_unique(),
                wf.unwrap_unique(),
                &mut result,
                &mut abs_err,
            )
        };
        Error::handle(ret, (result, abs_err))
    }
}

ffi_wrapper!(
    CquadWorkspace,
    *mut sys::gsl_integration_cquad_workspace,
    gsl_integration_cquad_workspace_free,
    "Workspace for the [`cquad`] integration algorithm."
);

impl CquadWorkspace {
    /// This function allocates a workspace sufficient to hold the
    /// data for `n` intervals.  The number `n` is not the maximum
    /// number of intervals that will be evaluated.  If the workspace
    /// is full, intervals with smaller error estimates will be
    /// discarded.  A minimum of 3 intervals is required and for most
    /// functions, a workspace of size 100 is sufficient.
    #[doc(alias = "gsl_integration_cquad_workspace_alloc")]
    pub fn new(n: usize) -> Option<CquadWorkspace> {
        let tmp = unsafe { sys::gsl_integration_cquad_workspace_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }
}

/// Doubly-adaptive integration.
///
/// CQUAD is a new doubly-adaptive general-purpose quadrature routine
/// which can handle most types of singularities, non-numerical
/// function values such as Inf or NaN, as well as some divergent
/// integrals. It generally requires more function evaluations than
/// the integration routines in QUADPACK, yet fails less often for
/// difficult integrands.
///
/// The CQUAD algorithm divides the integration region into
/// subintervals, and in each iteration, the subinterval with the
/// largest estimated error is processed.  The algorithm uses
/// Clenshaw-Curits quadrature rules of degree 4, 8, 16 and 32 over 5,
/// 9, 17 and 33 nodes respectively.  Each interval is initialized with
/// the lowest-degree rule.  When an interval is processed, the
/// next-higher degree rule is evaluated and an error estimate is
/// computed based on the $L^2$ -norm of the difference between the
/// underlying interpolating polynomials of both rules.  If the
/// highest-degree rule has already been used, or the interpolatory
/// polynomials differ significantly, the interval is bisected.
///
/// The subintervals and their results are stored in the memory
/// provided by workspace.
#[doc(alias = "gsl_integration_cquad")]
pub fn cquad<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> Cquad<'static, F> {
    Cquad::new(f, None, a, b)
}

integ_builder! {
    /// Builder to compute integrals using [`cquad`].
    Cquad<'w,>,
    workspace Option<BorrowedOrOwned<'w, CquadWorkspace>>,
    a f64, b f64,
}

impl<F: Fn(f64) -> f64> Cquad<'_, F> {
    /// Use the workspace `w` to compute the integral.
    pub fn workspace(self, w: &mut CquadWorkspace) -> Cquad<'_, F> {
        Cquad {
            workspace: Some(BorrowedOrOwned::Borrowed(w)),
            ..self
        }
    }

    /// Return $(∫_a^b f, \abserr, n)$.  The first component is an
    /// approximation of the integral, the second an estimate of the
    /// absolute error, and the third the number of function
    /// evaluations required.
    pub fn val_err_n(&mut self) -> Result<(f64, f64, usize), Error> {
        if self.workspace.is_none() {
            let w = CquadWorkspace::new(100).ok_or(Error::NoMemory)?;
            self.workspace = Some(BorrowedOrOwned::Owned(w));
        }
        let function = wrap_callback!(self.f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut n_evals = 0;

        let ret = unsafe {
            sys::gsl_integration_cquad(
                &function,
                self.a,
                self.b,
                self.epsabs,
                self.epsrel,
                self.workspace.as_mut().unwrap().to_mut().unwrap_unique(),
                &mut result,
                &mut abs_err,
                &mut n_evals,
            )
        };
        Error::handle(ret, (result, abs_err, n_evals))
    }
}

ffi_wrapper!(
    GLFixedTable,
    *mut sys::gsl_integration_glfixed_table,
    gsl_integration_glfixed_table_free,
    "Stores the Gauss-Legendre abscissae and weights"
);

impl GLFixedTable {
    /// This function determines the Gauss-Legendre abscissae and
    /// weights necessary for an `n`-point fixed order integration
    /// scheme.  If possible, high precision precomputed coefficients
    /// are used.  If precomputed weights are not available, lower
    /// precision coefficients are computed on the fly.
    #[doc(alias = "gsl_integration_glfixed_table_alloc")]
    pub fn new(n: usize) -> Option<GLFixedTable> {
        let tmp = unsafe { sys::gsl_integration_glfixed_table_alloc(n) };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// For `i` in {0, …, n - 1}, this function obtains the `i`-th
    /// Gauss-Legendre point `xi` and weight `wi` on the interval \[a,b\].
    /// The points and weights are ordered by increasing point value.
    /// A function f may be integrated on \[a,b\] by summing `wi` *
    /// f(`xi`) over `i`.
    ///
    /// Returns `(xi, wi)` if it succeeded.
    #[doc(alias = "gsl_integration_glfixed_point")]
    pub fn point(&self, a: f64, b: f64, i: usize) -> Result<(f64, f64), Error> {
        let mut xi = 0.;
        let mut wi = 0.;
        let ret = unsafe {
            sys::gsl_integration_glfixed_point(a, b, i, &mut xi, &mut wi, self.unwrap_shared())
        };
        Error::handle(ret, (xi, wi))
    }

    /// This function applies the Gauss-Legendre integration rule
    /// contained in table self and returns the result.
    #[doc(alias = "gsl_integration_glfixed")]
    pub fn glfixed<F: Fn(f64) -> f64>(&self, f: F, a: f64, b: f64) -> f64 {
        let function = wrap_callback!(f, F);
        unsafe { sys::gsl_integration_glfixed(&function, a, b, self.unwrap_shared()) }
    }
}

ffi_wrapper!(
    IntegrationFixedType,
    *const sys::gsl_integration_fixed_type,
    "Type of fixed point rule integration."
);

impl IntegrationFixedType {
    /// $w(x) = 1$ on the interval $(a,b)$.  Constraint $ b > a$.
    #[doc(alias = "gsl_integration_fixed_legendre")]
    pub fn legendre() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_legendre)
    }
    /// $w(x) = 1/ \sqrt{(b-x)(x-a)}$ on the interval $(a,b)$.
    /// Constraint $b > a$.
    #[doc(alias = "gsl_integration_fixed_chebyshev")]
    pub fn chebyshev() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_chebyshev)
    }
    /// $w(x) = \sqrt{(b - x) (x - a)}$ on the interval $(a,b)$.
    /// Constraint $b > a$.
    #[doc(alias = "gsl_integration_fixed_chebyshev2")]
    pub fn chebyshev2() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_chebyshev2)
    }
    /// $w(x) = ((b-x) (x-a))^α$ on the interval $(a,b)$.  Constraint
    /// $α > -1, b > a$.
    #[doc(alias = "gsl_integration_fixed_gegenbauer")]
    pub fn gegenbauer() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_gegenbauer)
    }
    /// $w(x) = (b-x)^α (x-a)^β$ on the interval $(a,b)$.  Constraint
    /// $α,β>-1, b > a$.
    #[doc(alias = "gsl_integration_fixed_jacobi")]
    pub fn jacobi() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_jacobi)
    }
    /// $w(x) = (x-a)^α \exp(-b(x-a))$ on the interval $(a,∞)$.
    /// Constraint $α > -1, b > 0$.
    #[doc(alias = "gsl_integration_fixed_laguerre")]
    pub fn laguerre() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_laguerre)
    }
    /// $w(x) = |x-a|^α \exp(-b(x-a)^2)$ on the interval $(-∞,∞)$.
    /// Constraint $α > -1, b > 0$.
    #[doc(alias = "gsl_integration_fixed_hermite")]
    pub fn hermite() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_hermite)
    }
    /// $w(x) = |x - (a + b)/2|^α$ on the interval $(a,b)$.
    /// Constraint $α > -1, b > a$.
    #[doc(alias = "gsl_integration_fixed_exponential")]
    pub fn exponential() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_exponential)
    }
    /// $w(x) = (x - a)^α (x + b)^β$ on the interval $(a,∞)$.
    /// Constraint $α > -1, α + β + 2n < 0, a + b > 0$.
    #[doc(alias = "gsl_integration_fixed_rational")]
    pub fn rational() -> IntegrationFixedType {
        ffi_wrap!(gsl_integration_fixed_rational)
    }
}

ffi_wrapper!(
    IntegrationFixedWorkspace,
    *mut sys::gsl_integration_fixed_workspace,
    gsl_integration_fixed_free,
    "Workspace for fixed point quadratures."
);

impl IntegrationFixedWorkspace {
    /// This function allocates a workspace for computing integrals
    /// with interpolating quadratures using n quadrature nodes.  The
    /// parameters `a`, `b`, `alpha`, and `beta` specify the
    /// integration interval and/or weighting function for the various
    /// quadrature types.  See the table at the beginning of
    /// [`integration`][crate::integration] or the methods in
    /// [`IntegrationFixedType`] for constraints on these parameters.
    /// The size of the workspace is O(4`n`).
    #[doc(alias = "gsl_integration_fixed_alloc")]
    pub fn new(
        type_: IntegrationFixedType,
        n: usize,
        a: f64,
        b: f64,
        alpha: f64,
        beta: f64,
    ) -> Option<IntegrationFixedWorkspace> {
        let tmp = unsafe {
            sys::gsl_integration_fixed_alloc(type_.unwrap_shared(), n, a, b, alpha, beta)
        };

        if tmp.is_null() {
            None
        } else {
            Some(Self::wrap(tmp))
        }
    }

    /// Return the number of quadrature nodes and weights.
    #[doc(alias = "gsl_integration_fixed_n")]
    pub fn n(&self) -> usize {
        unsafe { sys::gsl_integration_fixed_n(self.unwrap_shared()) }
    }

    /// Return a slice of size $n$ containing the quadrature nodes $x_i$.
    #[doc(alias = "gsl_integration_fixed_nodes")]
    pub fn nodes(&self) -> &[f64] {
        let tmp = unsafe { sys::gsl_integration_fixed_nodes(self.unwrap_shared()) };
        assert!(!tmp.is_null());
        unsafe { std::slice::from_raw_parts(tmp, self.n()) }
    }

    /// Return a slice of size $n$ containing the quadrature weights $w_i$.
    #[doc(alias = "gsl_integration_fixed_weights")]
    pub fn weights(&self) -> &[f64] {
        let tmp = unsafe { sys::gsl_integration_fixed_weights(self.unwrap_shared()) };
        assert!(!tmp.is_null());
        unsafe { std::slice::from_raw_parts(tmp, self.n()) }
    }

    /// Return the integral of the function `f` using previously
    /// computed fixed quadrature rules. The integral is approximated
    /// as
    ///
    /// $$\sum_{i=1}^n w_i f(x_i)$$
    ///
    /// where $w_i$ are the quadrature weights and $x_i$ are the
    /// quadrature nodes computed in `self`.
    #[doc(alias = "gsl_integration_fixed")]
    pub fn fixed<F: Fn(f64) -> f64>(&self, f: F) -> Result<f64, Error> {
        let mut result = 0.;
        let function = wrap_callback!(f, F);

        let ret =
            unsafe { sys::gsl_integration_fixed(&function, &mut result, self.unwrap_shared()) };
        Error::handle(ret, result)
    }
}

// TODO: integration on the unit sphere

/// Basic Gauss-Kronrod quadratures.
pub mod basic {
    /// 15-point Gauss-Kronrod quadrature.
    ///
    /// Gauss quadrature weights and kronrod quadrature abscissae and
    /// weights as evaluated with 80 decimal digit arithmetic by L. W.
    /// Fullerton, Bell Labs, Nov. 1981.
    ///
    /// Returns `(result, abs_err, resabs, resasc)`.
    #[doc(alias = "gsl_integration_qk15")]
    pub fn qk15<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
        let function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut resabs = 0.;
        let mut resasc = 0.;

        unsafe {
            sys::gsl_integration_qk15(
                &function,
                a,
                b,
                &mut result,
                &mut abs_err,
                &mut resabs,
                &mut resasc,
            )
        };
        (result, abs_err, resabs, resasc)
    }

    /// 21-point Gauss-Kronrod quadrature.
    ///
    /// Returns `(result, abs_err, resabs, resasc)`.
    #[doc(alias = "gsl_integration_qk21")]
    pub fn qk21<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
        let function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut resabs = 0.;
        let mut resasc = 0.;

        unsafe {
            sys::gsl_integration_qk21(
                &function,
                a,
                b,
                &mut result,
                &mut abs_err,
                &mut resabs,
                &mut resasc,
            )
        };
        (result, abs_err, resabs, resasc)
    }

    /// 31-point Gauss-Kronrod quadrature.
    ///
    /// Returns `(result, abs_err, resabs, resasc)`.
    #[doc(alias = "gsl_integration_qk31")]
    pub fn qk31<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
        let function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut resabs = 0.;
        let mut resasc = 0.;

        unsafe {
            sys::gsl_integration_qk31(
                &function,
                a,
                b,
                &mut result,
                &mut abs_err,
                &mut resabs,
                &mut resasc,
            )
        };
        (result, abs_err, resabs, resasc)
    }

    /// 41-point Gauss-Kronrod quadrature.
    ///
    /// Returns `(result, abs_err, resabs, resasc)`.
    #[doc(alias = "gsl_integration_qk41")]
    pub fn qk41<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
        let function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut resabs = 0.;
        let mut resasc = 0.;

        unsafe {
            sys::gsl_integration_qk41(
                &function,
                a,
                b,
                &mut result,
                &mut abs_err,
                &mut resabs,
                &mut resasc,
            )
        };
        (result, abs_err, resabs, resasc)
    }

    /// 51-point Gauss-Kronrod quadrature.
    ///
    /// Returns `(result, abs_err, resabs, resasc)`.
    #[doc(alias = "gsl_integration_qk51")]
    pub fn qk51<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
        let function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut resabs = 0.;
        let mut resasc = 0.;

        unsafe {
            sys::gsl_integration_qk51(
                &function,
                a,
                b,
                &mut result,
                &mut abs_err,
                &mut resabs,
                &mut resasc,
            )
        };
        (result, abs_err, resabs, resasc)
    }

    /// 61-point Gauss-Kronrod quadrature.
    ///
    /// Returns `(result, abs_err, resabs, resasc)`.
    #[doc(alias = "gsl_integration_qk61")]
    pub fn qk61<F: Fn(f64) -> f64>(f: F, a: f64, b: f64) -> (f64, f64, f64, f64) {
        let function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut resabs = 0.;
        let mut resasc = 0.;

        unsafe {
            sys::gsl_integration_qk61(
                &function,
                a,
                b,
                &mut result,
                &mut abs_err,
                &mut resabs,
                &mut resasc,
            )
        };
        (result, abs_err, resabs, resasc)
    }

    /// Returns `(result, abs_err, resabs, resasc)`.
    #[doc(alias = "gsl_integration_qk")]
    pub fn qk<F: Fn(f64) -> f64>(
        xgk: &[f64],
        wg: &[f64],
        wgk: &[f64],
        fv1: &mut [f64],
        fv2: &mut [f64],
        f: F,
        a: f64,
        b: f64,
    ) -> (f64, f64, f64, f64) {
        assert!(xgk.len() == wg.len());
        assert!(xgk.len() == wgk.len());
        assert!(xgk.len() == fv1.len());
        assert!(xgk.len() == fv2.len());

        let function = wrap_callback!(f, F);
        let mut result = 0.;
        let mut abs_err = 0.;
        let mut resabs = 0.;
        let mut resasc = 0.;

        unsafe {
            sys::gsl_integration_qk(
                xgk.len() as _,
                xgk.as_ptr(),
                wg.as_ptr(),
                wgk.as_ptr(),
                fv1.as_mut_ptr(),
                fv2.as_mut_ptr(),
                &function,
                a,
                b,
                &mut result,
                &mut abs_err,
                &mut resabs,
                &mut resasc,
            );
        }
        (result, abs_err, resabs, resasc)
    }
}
