//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
##Introduction

Each algorithm computes an approximation to a definite integral of the form,

I = \int_a^b f(x) w(x) dx
where w(x) is a weight function (for general integrands w(x)=1). The user provides absolute and relative error bounds (epsabs, epsrel) which
specify the following accuracy requirement,

|RESULT - I|  <= max(epsabs, epsrel |I|)

where RESULT is the numerical approximation obtained by the algorithm. The algorithms attempt to estimate the absolute error ABSERR = |RESULT
- I| in such a way that the following inequality holds,

|RESULT - I| <= ABSERR <= max(epsabs, epsrel |I|)

In short, the routines return the first approximation which has an absolute error smaller than epsabs or a relative error smaller than epsrel.

Note that this is an either-or constraint, not simultaneous. To compute to a specified absolute error, set epsrel to zero. To compute to a
specified relative error, set epsabs to zero. The routines will fail to converge if the error bounds are too stringent, but always return the
best approximation obtained up to that stage.

The algorithms in QUADPACK use a naming convention based on the following letters,

Q - quadrature routine

N - non-adaptive integrator
A - adaptive integrator

G - general integrand (user-defined)
W - weight function with integrand

S - singularities can be more readily integrated
P - points of special difficulty can be supplied
I - infinite range of integration
O - oscillatory weight function, cos or sin
F - Fourier integral
C - Cauchy principal value
The algorithms are built on pairs of quadrature rules, a higher order rule and a lower order rule. The higher order rule is used to compute the
best approximation to an integral over a small range. The difference between the results of the higher order rule and the lower order rule gives
an estimate of the error in the approximation.

 * [Integrands without weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-without-weight-functions.html#Integrands-without-weight-functions)
 * [Integrands with weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-with-weight-functions.html#Integrands-with-weight-functions)
 * [Integrands with singular weight functions](http://www.gnu.org/software/gsl/manual/html_node/Integrands-with-singular-weight-functions.html#Integrands-with-singular-weight-functions)

##QNG non-adaptive Gauss-Kronrod integration

The QNG algorithm is a non-adaptive procedure which uses fixed Gauss-Kronrod-Patterson abscissae to sample the integrand at a maximum of 87
points. It is provided for fast integration of smooth functions.

##QAG adaptive integration

The QAG algorithm is a simple adaptive integration procedure. The integration region is divided into subintervals, and on each iteration the
subinterval with the largest estimated error is bisected. This reduces the overall error rapidly, as the subintervals become concentrated
around local difficulties in the integrand. These subintervals are managed by a gsl_integration_workspace struct, which handles the memory
for the subinterval ranges, results and error estimates.

##QAGS adaptive integration with singularities

The presence of an integrable singularity in the integration region causes an adaptive routine to concentrate new subintervals around the
singularity. As the subintervals decrease in size the successive approximations to the integral converge in a limiting fashion. This
approach to the limit can be accelerated using an extrapolation procedure. The QAGS algorithm combines adaptive bisection with the Wynn
epsilon-algorithm to speed up the integration of many types of integrable singularities.

##References and Further Reading

The following book is the definitive reference for QUADPACK, and was written by the original authors. It provides descriptions of the
algorithms, program listings, test programs and examples. It also includes useful advice on numerical integration and many references
to the numerical integration literature used in developing QUADPACK.

R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K. Kahaner. QUADPACK A subroutine package for automatic integration Springer Verlag, 1983.
The CQUAD integration algorithm is described in the following paper:

P. Gonnet, “Increasing the Reliability of Adaptive Quadrature Using Explicit Interpolants”, ACM Transactions on Mathematical Software, Volume 37
(2010), Issue 3, Article 26.
!*/

use enums;
use ffi;
use std::ffi::CString;

fn rescale_error(err: f64, result_abs: f64, result_asc: f64) -> f64 {
    let mut t_err = unsafe { err.abs() };

    if result_asc != 0f64 && t_err != 0f64 {
        let scale = unsafe { (200f64 * t_err / result_asc).powf(1.5f64) };

        if scale < 1f64 {
            t_err = result_asc * scale;
        } else {
            t_err = result_asc;
        }
    }
    if result_abs > ::DBL_MIN / (50f64 * ::DBL_EPSILON) {
        let min_err = 50f64 * ::DBL_EPSILON * result_abs;

        if min_err > t_err {
            t_err = min_err;
        }
    }

    t_err
}

/// This function applies the Gauss-Kronrod 10-point, 21-point, 43-point and 87-point integration rules in succession until an estimate of the
/// integral of f over (a,b) is achieved within the desired absolute and relative error limits, eps_abs and eps_rel. The function returns the final
/// approximation, result, an estimate of the absolute error, abserr and the number of function evaluations used, neval. The Gauss-Kronrod rules
/// are designed in such a way that each rule uses all the results of its predecessors, in order to minimize the total number of function
/// evaluations.
pub fn qng<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    eps_abs: f64,
    eps_rel: f64,
    result: &mut f64,
    abs_err: &mut f64,
    n_eval: &mut usize,
) -> ::Value {
    let half_length = 0.5f64 * (b - a);
    let abs_half_length = unsafe { half_length.abs() };
    let center = 0.5f64 * (b + a);
    let f_center = f(center, arg);

    if eps_abs <= 0f64 && (eps_rel < 50f64 * ::DBL_EPSILON || eps_rel < 0.5e-28f64) {
        *result = 0f64;
        *abs_err = 0f64;
        *n_eval = 0usize;

        //to avoid the dead_code warning on gsl_error
        unsafe {
            let file = file!();
            let c_str = CString::new(
                "tolerance cannot be achieved with given eps_abs and eps_rel".as_bytes(),
            )
            .unwrap();
            let c_file = CString::new(file.as_bytes()).unwrap();

            unsafe {
                sys::gsl_error(
                    c_str.as_ptr(),
                    c_file.as_ptr(),
                    line!() as i32,
                    ::Value::BadTolerance.into(),
                );
            }
        }
        return ::Value::BadTolerance;
    }

    // w21b, weights of the 21-point formula for abscissae x2
    let w21b: [f64; 6] = [
        0.011694638867371874278064396062192f64,
        0.054755896574351996031381300244580f64,
        0.093125454583697605535065465083366f64,
        0.123491976262065851077958109831074f64,
        0.142775938577060080797094273138717f64,
        0.149445554002916905664936468389821f64,
    ];
    // w10, weights of the 10-point formula
    let w10: [f64; 5] = [
        0.066671344308688137593568809893332f64,
        0.149451349150580593145776339657697f64,
        0.219086362515982043995534934228163f64,
        0.269266719309996355091226921569469f64,
        0.295524224714752870173892994651338f64,
    ];
    // w21a, weights of the 21-point formula for abscissae x1 */
    let w21a: [f64; 5] = [
        0.032558162307964727478818972459390f64,
        0.075039674810919952767043140916190f64,
        0.109387158802297641899210590325805f64,
        0.134709217311473325928054001771707f64,
        0.147739104901338491374841515972068f64,
    ];
    // x1, abscissae common to the 10-, 21-, 43- and 87-point rule */
    let x1: [f64; 5] = [
        0.973906528517171720077964012084452f64,
        0.865063366688984510732096688423493f64,
        0.679409568299024406234327365114874f64,
        0.433395394129247190799265943165784f64,
        0.148874338981631210884826001129720f64,
    ];
    // x2, abscissae common to the 21-, 43- and 87-point rule
    let x2: [f64; 5] = [
        0.995657163025808080735527280689003f64,
        0.930157491355708226001207180059508f64,
        0.780817726586416897063717578345042f64,
        0.562757134668604683339000099272694f64,
        0.294392862701460198131126603103866f64,
    ];
    // x3, abscissae common to the 43- and 87-point rule */
    let x3: [f64; 11] = [
        0.999333360901932081394099323919911f64,
        0.987433402908088869795961478381209f64,
        0.954807934814266299257919200290473f64,
        0.900148695748328293625099494069092f64,
        0.825198314983114150847066732588520f64,
        0.732148388989304982612354848755461f64,
        0.622847970537725238641159120344323f64,
        0.499479574071056499952214885499755f64,
        0.364901661346580768043989548502644f64,
        0.222254919776601296498260928066212f64,
        0.074650617461383322043914435796506f64,
    ];
    // x4, abscissae of the 87-point rule */
    let x4: [f64; 22] = [
        0.999902977262729234490529830591582f64,
        0.997989895986678745427496322365960f64,
        0.992175497860687222808523352251425f64,
        0.981358163572712773571916941623894f64,
        0.965057623858384619128284110607926f64,
        0.943167613133670596816416634507426f64,
        0.915806414685507209591826430720050f64,
        0.883221657771316501372117548744163f64,
        0.845710748462415666605902011504855f64,
        0.803557658035230982788739474980964f64,
        0.757005730685495558328942793432020f64,
        0.706273209787321819824094274740840f64,
        0.651589466501177922534422205016736f64,
        0.593223374057961088875273770349144f64,
        0.531493605970831932285268948562671f64,
        0.466763623042022844871966781659270f64,
        0.399424847859218804732101665817923f64,
        0.329874877106188288265053371824597f64,
        0.258503559202161551802280975429025f64,
        0.185695396568346652015917141167606f64,
        0.111842213179907468172398359241362f64,
        0.037352123394619870814998165437704f64,
    ];
    // w43a, weights of the 43-point formula for abscissae x1, x3 */
    let w43a: [f64; 10] = [
        0.016296734289666564924281974617663f64,
        0.037522876120869501461613795898115f64,
        0.054694902058255442147212685465005f64,
        0.067355414609478086075553166302174f64,
        0.073870199632393953432140695251367f64,
        0.005768556059769796184184327908655f64,
        0.027371890593248842081276069289151f64,
        0.046560826910428830743339154433824f64,
        0.061744995201442564496240336030883f64,
        0.071387267268693397768559114425516f64,
    ];
    // w43b, weights of the 43-point formula for abscissae x3 */
    let w43b: [f64; 12] = [
        0.001844477640212414100389106552965f64,
        0.010798689585891651740465406741293f64,
        0.021895363867795428102523123075149f64,
        0.032597463975345689443882222526137f64,
        0.042163137935191811847627924327955f64,
        0.050741939600184577780189020092084f64,
        0.058379395542619248375475369330206f64,
        0.064746404951445885544689259517511f64,
        0.069566197912356484528633315038405f64,
        0.072824441471833208150939535192842f64,
        0.074507751014175118273571813842889f64,
        0.074722147517403005594425168280423f64,
    ];
    // w87a, weights of the 87-point formula for abscissae x1, x2, x3
    let w87a: [f64; 21] = [
        0.008148377384149172900002878448190f64,
        0.018761438201562822243935059003794f64,
        0.027347451050052286161582829741283f64,
        0.033677707311637930046581056957588f64,
        0.036935099820427907614589586742499f64,
        0.002884872430211530501334156248695f64,
        0.013685946022712701888950035273128f64,
        0.023280413502888311123409291030404f64,
        0.030872497611713358675466394126442f64,
        0.035693633639418770719351355457044f64,
        0.000915283345202241360843392549948f64,
        0.005399280219300471367738743391053f64,
        0.010947679601118931134327826856808f64,
        0.016298731696787335262665703223280f64,
        0.021081568889203835112433060188190f64,
        0.025370969769253827243467999831710f64,
        0.029189697756475752501446154084920f64,
        0.032373202467202789685788194889595f64,
        0.034783098950365142750781997949596f64,
        0.036412220731351787562801163687577f64,
        0.037253875503047708539592001191226f64,
    ];
    // w87b, weights of the 87-point formula for abscissae x4
    let w87b: [f64; 23] = [
        0.000274145563762072350016527092881f64,
        0.001807124155057942948341311753254f64,
        0.004096869282759164864458070683480f64,
        0.006758290051847378699816577897424f64,
        0.009549957672201646536053581325377f64,
        0.012329447652244853694626639963780f64,
        0.015010447346388952376697286041943f64,
        0.017548967986243191099665352925900f64,
        0.019938037786440888202278192730714f64,
        0.022194935961012286796332102959499f64,
        0.024339147126000805470360647041454f64,
        0.026374505414839207241503786552615f64,
        0.028286910788771200659968002987960f64,
        0.030052581128092695322521110347341f64,
        0.031646751371439929404586051078883f64,
        0.033050413419978503290785944862689f64,
        0.034255099704226061787082821046821f64,
        0.035262412660156681033782717998428f64,
        0.036076989622888701185500318003895f64,
        0.036698604498456094498018047441094f64,
        0.037120549269832576114119958413599f64,
        0.037334228751935040321235449094698f64,
        0.037361073762679023410321241766599f64,
    ];

    // Compute the integral using the 10- and 21-point formula.
    let mut res10 = 0f64;
    let mut res21 = w21b[5] * f_center;
    let mut resabs = unsafe { w21b[5] * f_center.abs() };
    let mut savfun: [f64; 21] = [0f64; 21];
    let mut fv1: [f64; 5] = [0f64; 5];
    let mut fv2: [f64; 5] = [0f64; 5];
    let mut fv3: [f64; 5] = [0f64; 5];
    let mut fv4: [f64; 5] = [0f64; 5];

    for k in 0usize..5usize {
        let abscissa = half_length * x1[k];
        let fval1 = f(center + abscissa, arg);
        let fval2 = f(center - abscissa, arg);
        let fval = fval1 + fval2;

        res10 += w10[k] * fval;
        res21 += w21a[k] * fval;
        resabs += unsafe { w21a[k] * (fval1.abs() + fval2.abs()) };
        savfun[k] = fval;
        fv1[k] = fval1;
        fv2[k] = fval2;
    }

    for k in 0usize..5usize {
        let abscissa = half_length * x2[k];
        let fval1 = f(center + abscissa, arg);
        let fval2 = f(center - abscissa, arg);
        let fval = fval1 + fval2;

        res21 += w21b[k] * fval;
        resabs += unsafe { w21b[k] * (fval1.abs() + fval2.abs()) };
        savfun[k + 5] = fval;
        fv3[k] = fval1;
        fv4[k] = fval2;
    }
    resabs *= abs_half_length;

    let mean = 0.5f64 * res21;
    let mut resasc = unsafe { w21b[5] * (f_center - mean).abs() };

    for k in 0usize..5usize {
        resasc += unsafe {
            w21a[k] * ((fv1[k] - mean).abs() + (fv2[k] - mean).abs())
                + w21b[k] * ((fv3[k] - mean).abs() + (fv4[k] - mean).abs())
        };
    }
    resasc *= abs_half_length;
    let mut result_kronrod = res21 * half_length;
    let mut err = rescale_error((res21 - res10) * half_length, resabs, resasc);

    // test for convergence.
    if err < eps_abs || err < eps_rel * unsafe { result_kronrod.abs() } {
        *result = result_kronrod;
        *abs_err = err;
        *n_eval = 21;
        return ::Value::Success;
    }

    // compute the integral using the 43-point formula.
    let mut res43 = w43b[11] * f_center;

    for k in 0usize..10usize {
        res43 += savfun[k] * w43a[k];
    }

    for k in 0usize..11usize {
        let abscissa = half_length * x3[k];
        let fval = f(center + abscissa, arg) + f(center - abscissa, arg);

        res43 += fval * w43b[k];
        savfun[k + 10] = fval;
    }

    // test for convergence
    result_kronrod = res43 * half_length;
    err = rescale_error((res43 - res21) * half_length, resabs, resasc);

    if err < eps_abs || err < eps_rel * unsafe { result_kronrod.abs() } {
        *result = result_kronrod;
        *abs_err = err;
        *n_eval = 43;
        return ::Value::Success;
    }

    // compute the integral using the 87-point formula.
    let mut res87 = w87b[22] * f_center;

    for k in 0..21 {
        res87 += savfun[k] * w87a[k];
    }
    for k in 0..22 {
        let abscissa = half_length * x4[k];

        res87 += w87b[k] * f(center + abscissa, arg) + f(center - abscissa, arg);
    }

    // test for convergence
    result_kronrod = res87 * half_length;
    err = rescale_error((res87 - res43) * half_length, resabs, resasc);

    if err < eps_abs || err < eps_rel * unsafe { result_kronrod.abs() } {
        *result = result_kronrod;
        *abs_err = err;
        *n_eval = 87;
        return ::Value::Success;
    }

    // failed to converge
    *result = result_kronrod;
    *abs_err = err;
    *n_eval = 87;
    rgsl_error!(
        "failed to reach tolerance with highest-order rule",
        ::Value::Tolerance
    );
    ::Value::Tolerance
}

/// Gauss quadrature weights and kronrod quadrature abscissae and weights as evaluated with 80 decimal digit arithmetic by L. W.
/// Fullerton, Bell Labs, Nov. 1981.
pub fn qk15<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    result: &mut f64,
    abserr: &mut f64,
    resabs: &mut f64,
    resasc: &mut f64,
) {
    // abscissae of the 15-point kronrod rule
    let xgk: [f64; 8] = [
        0.991455371120812639206854697526329f64,
        0.949107912342758524526189684047851f64,
        0.864864423359769072789712788640926f64,
        0.741531185599394439863864773280788f64,
        0.586087235467691130294144838258730f64,
        0.405845151377397166906606412076961f64,
        0.207784955007898467600689403773245f64,
        0.000000000000000000000000000000000f64,
    ];

    // xgk[1], xgk[3], ... abscissae of the 7-point gauss rule.
    // xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule

    // weights of the 7-point gauss rule
    let wg: [f64; 4] = [
        0.129484966168869693270611432679082f64,
        0.279705391489276667901467771423780f64,
        0.381830050505118944950369775488975f64,
        0.417959183673469387755102040816327f64,
    ];

    // weights of the 15-point kronrod rule
    let wgk: [f64; 8] = [
        0.022935322010529224963732008058970f64,
        0.063092092629978553290700663189204f64,
        0.104790010322250183839876322541518f64,
        0.140653259715525918745189590510238f64,
        0.169004726639267902826583426598550f64,
        0.190350578064785409913256402421014f64,
        0.204432940075298892414161999234649f64,
        0.209482141084727828012999174891714f64,
    ];

    let mut fv1: [f64; 8] = [0f64; 8];
    let mut fv2: [f64; 8] = [0f64; 8];

    qk(
        &xgk, &wg, &wgk, &mut fv1, &mut fv2, f, arg, a, b, result, abserr, resabs, resasc,
    );
}

pub fn qk21<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    result: &mut f64,
    abserr: &mut f64,
    resabs: &mut f64,
    resasc: &mut f64,
) {
    // abscissae of the 21-point kronrod rule
    let xgk: [f64; 11] = [
        0.995657163025808080735527280689003f64,
        0.973906528517171720077964012084452f64,
        0.930157491355708226001207180059508f64,
        0.865063366688984510732096688423493f64,
        0.780817726586416897063717578345042f64,
        0.679409568299024406234327365114874f64,
        0.562757134668604683339000099272694f64,
        0.433395394129247190799265943165784f64,
        0.294392862701460198131126603103866f64,
        0.148874338981631210884826001129720f64,
        0.000000000000000000000000000000000f64,
    ];

    // xgk[1], xgk[3], ... abscissae of the 10-point gauss rule.
    // xgk[0], xgk[2], ... abscissae to optimally extend the 10-point gauss rule

    // weights of the 10-point gauss rule
    let wg: [f64; 5] = [
        0.066671344308688137593568809893332f64,
        0.149451349150580593145776339657697f64,
        0.219086362515982043995534934228163f64,
        0.269266719309996355091226921569469f64,
        0.295524224714752870173892994651338f64,
    ];

    // weights of the 21-point kronrod rule
    let wgk: [f64; 11] = [
        0.011694638867371874278064396062192f64,
        0.032558162307964727478818972459390f64,
        0.054755896574351996031381300244580f64,
        0.075039674810919952767043140916190f64,
        0.093125454583697605535065465083366f64,
        0.109387158802297641899210590325805f64,
        0.123491976262065851077958109831074f64,
        0.134709217311473325928054001771707f64,
        0.142775938577060080797094273138717f64,
        0.147739104901338491374841515972068f64,
        0.149445554002916905664936468389821f64,
    ];

    let mut fv1: [f64; 11] = [0f64; 11];
    let mut fv2: [f64; 11] = [0f64; 11];

    qk(
        &xgk, &wg, &wgk, &mut fv1, &mut fv2, f, arg, a, b, result, abserr, resabs, resasc,
    );
}

pub fn qk31<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    result: &mut f64,
    abserr: &mut f64,
    resabs: &mut f64,
    resasc: &mut f64,
) {
    // abscissae of the 31-point kronrod rule
    let xgk: [f64; 16] = [
        0.998002298693397060285172840152271f64,
        0.987992518020485428489565718586613f64,
        0.967739075679139134257347978784337f64,
        0.937273392400705904307758947710209f64,
        0.897264532344081900882509656454496f64,
        0.848206583410427216200648320774217f64,
        0.790418501442465932967649294817947f64,
        0.724417731360170047416186054613938f64,
        0.650996741297416970533735895313275f64,
        0.570972172608538847537226737253911f64,
        0.485081863640239680693655740232351f64,
        0.394151347077563369897207370981045f64,
        0.299180007153168812166780024266389f64,
        0.201194093997434522300628303394596f64,
        0.101142066918717499027074231447392f64,
        0.000000000000000000000000000000000f64,
    ];

    // xgk[1], xgk[3], ... abscissae of the 15-point gauss rule.
    // xgk[0], xgk[2], ... abscissae to optimally extend the 15-point gauss rule

    // weights of the 15-point gauss rule
    let wg: [f64; 8] = [
        0.030753241996117268354628393577204f64,
        0.070366047488108124709267416450667f64,
        0.107159220467171935011869546685869f64,
        0.139570677926154314447804794511028f64,
        0.166269205816993933553200860481209f64,
        0.186161000015562211026800561866423f64,
        0.198431485327111576456118326443839f64,
        0.202578241925561272880620199967519f64,
    ];

    // weights of the 31-point kronrod rule
    let wgk: [f64; 16] = [
        0.005377479872923348987792051430128f64,
        0.015007947329316122538374763075807f64,
        0.025460847326715320186874001019653f64,
        0.035346360791375846222037948478360f64,
        0.044589751324764876608227299373280f64,
        0.053481524690928087265343147239430f64,
        0.062009567800670640285139230960803f64,
        0.069854121318728258709520077099147f64,
        0.076849680757720378894432777482659f64,
        0.083080502823133021038289247286104f64,
        0.088564443056211770647275443693774f64,
        0.093126598170825321225486872747346f64,
        0.096642726983623678505179907627589f64,
        0.099173598721791959332393173484603f64,
        0.100769845523875595044946662617570f64,
        0.101330007014791549017374792767493f64,
    ];

    let mut fv1: [f64; 16] = [0f64; 16];
    let mut fv2: [f64; 16] = [0f64; 16];

    qk(
        &xgk, &wg, &wgk, &mut fv1, &mut fv2, f, arg, a, b, result, abserr, resabs, resasc,
    );
}

pub fn qk41<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    result: &mut f64,
    abserr: &mut f64,
    resabs: &mut f64,
    resasc: &mut f64,
) {
    // abscissae of the 41-point kronrod rule
    let xgk: [f64; 21] = [
        0.998859031588277663838315576545863f64,
        0.993128599185094924786122388471320f64,
        0.981507877450250259193342994720217f64,
        0.963971927277913791267666131197277f64,
        0.940822633831754753519982722212443f64,
        0.912234428251325905867752441203298f64,
        0.878276811252281976077442995113078f64,
        0.839116971822218823394529061701521f64,
        0.795041428837551198350638833272788f64,
        0.746331906460150792614305070355642f64,
        0.693237656334751384805490711845932f64,
        0.636053680726515025452836696226286f64,
        0.575140446819710315342946036586425f64,
        0.510867001950827098004364050955251f64,
        0.443593175238725103199992213492640f64,
        0.373706088715419560672548177024927f64,
        0.301627868114913004320555356858592f64,
        0.227785851141645078080496195368575f64,
        0.152605465240922675505220241022678f64,
        0.076526521133497333754640409398838f64,
        0.000000000000000000000000000000000f64,
    ];

    // xgk[1], xgk[3], ... abscissae of the 20-point gauss rule.
    // xgk[0], xgk[2], ... abscissae to optimally extend the 20-point gauss rule

    // weights of the 20-point gauss rule
    let wg: [f64; 11] = [
        0.017614007139152118311861962351853f64,
        0.040601429800386941331039952274932f64,
        0.062672048334109063569506535187042f64,
        0.083276741576704748724758143222046f64,
        0.101930119817240435036750135480350f64,
        0.118194531961518417312377377711382f64,
        0.131688638449176626898494499748163f64,
        0.142096109318382051329298325067165f64,
        0.149172986472603746787828737001969f64,
        0.152753387130725850698084331955098f64,
        0.000000000000000000000000000000000f64,
    ];

    // weights of the 41-point kronrod rule
    let wgk: [f64; 21] = [
        0.003073583718520531501218293246031f64,
        0.008600269855642942198661787950102f64,
        0.014626169256971252983787960308868f64,
        0.020388373461266523598010231432755f64,
        0.025882133604951158834505067096153f64,
        0.031287306777032798958543119323801f64,
        0.036600169758200798030557240707211f64,
        0.041668873327973686263788305936895f64,
        0.046434821867497674720231880926108f64,
        0.050944573923728691932707670050345f64,
        0.055195105348285994744832372419777f64,
        0.059111400880639572374967220648594f64,
        0.062653237554781168025870122174255f64,
        0.065834597133618422111563556969398f64,
        0.068648672928521619345623411885368f64,
        0.071054423553444068305790361723210f64,
        0.073030690332786667495189417658913f64,
        0.074582875400499188986581418362488f64,
        0.075704497684556674659542775376617f64,
        0.076377867672080736705502835038061f64,
        0.076600711917999656445049901530102f64,
    ];

    let mut fv1: [f64; 21] = [0f64; 21];
    let mut fv2: [f64; 21] = [0f64; 21];

    qk(
        &xgk, &wg, &wgk, &mut fv1, &mut fv2, f, arg, a, b, result, abserr, resabs, resasc,
    );
}

pub fn qk51<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    result: &mut f64,
    abserr: &mut f64,
    resabs: &mut f64,
    resasc: &mut f64,
) {
    // abscissae of the 51-point kronrod rule
    let xgk: [f64; 26] = [
        0.999262104992609834193457486540341f64,
        0.995556969790498097908784946893902f64,
        0.988035794534077247637331014577406f64,
        0.976663921459517511498315386479594f64,
        0.961614986425842512418130033660167f64,
        0.942974571228974339414011169658471f64,
        0.920747115281701561746346084546331f64,
        0.894991997878275368851042006782805f64,
        0.865847065293275595448996969588340f64,
        0.833442628760834001421021108693570f64,
        0.797873797998500059410410904994307f64,
        0.759259263037357630577282865204361f64,
        0.717766406813084388186654079773298f64,
        0.673566368473468364485120633247622f64,
        0.626810099010317412788122681624518f64,
        0.577662930241222967723689841612654f64,
        0.526325284334719182599623778158010f64,
        0.473002731445714960522182115009192f64,
        0.417885382193037748851814394594572f64,
        0.361172305809387837735821730127641f64,
        0.303089538931107830167478909980339f64,
        0.243866883720988432045190362797452f64,
        0.183718939421048892015969888759528f64,
        0.122864692610710396387359818808037f64,
        0.061544483005685078886546392366797f64,
        0.000000000000000000000000000000000f64,
    ];

    // xgk[1], xgk[3], ... abscissae of the 25-point gauss rule.
    // xgk[0], xgk[2], ... abscissae to optimally extend the 25-point gauss rule

    // weights of the 25-point gauss rule
    let wg: [f64; 13] = [
        0.011393798501026287947902964113235f64,
        0.026354986615032137261901815295299f64,
        0.040939156701306312655623487711646f64,
        0.054904695975835191925936891540473f64,
        0.068038333812356917207187185656708f64,
        0.080140700335001018013234959669111f64,
        0.091028261982963649811497220702892f64,
        0.100535949067050644202206890392686f64,
        0.108519624474263653116093957050117f64,
        0.114858259145711648339325545869556f64,
        0.119455763535784772228178126512901f64,
        0.122242442990310041688959518945852f64,
        0.123176053726715451203902873079050f64,
    ];

    // weights of the 51-point kronrod rule
    let wgk: [f64; 26] = [
        0.001987383892330315926507851882843f64,
        0.005561932135356713758040236901066f64,
        0.009473973386174151607207710523655f64,
        0.013236229195571674813656405846976f64,
        0.016847817709128298231516667536336f64,
        0.020435371145882835456568292235939f64,
        0.024009945606953216220092489164881f64,
        0.027475317587851737802948455517811f64,
        0.030792300167387488891109020215229f64,
        0.034002130274329337836748795229551f64,
        0.037116271483415543560330625367620f64,
        0.040083825504032382074839284467076f64,
        0.042872845020170049476895792439495f64,
        0.045502913049921788909870584752660f64,
        0.047982537138836713906392255756915f64,
        0.050277679080715671963325259433440f64,
        0.052362885806407475864366712137873f64,
        0.054251129888545490144543370459876f64,
        0.055950811220412317308240686382747f64,
        0.057437116361567832853582693939506f64,
        0.058689680022394207961974175856788f64,
        0.059720340324174059979099291932562f64,
        0.060539455376045862945360267517565f64,
        0.061128509717053048305859030416293f64,
        0.061471189871425316661544131965264f64,
        0.061580818067832935078759824240066f64,
    ];

    let mut fv1: [f64; 26] = [0f64; 26];
    let mut fv2: [f64; 26] = [0f64; 26];

    qk(
        &xgk, &wg, &wgk, &mut fv1, &mut fv2, f, arg, a, b, result, abserr, resabs, resasc,
    );
}

pub fn qk61<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    result: &mut f64,
    abserr: &mut f64,
    resabs: &mut f64,
    resasc: &mut f64,
) {
    // abscissae of the 61-point kronrod rule
    let xgk: [f64; 31] = [
        0.999484410050490637571325895705811f64,
        0.996893484074649540271630050918695f64,
        0.991630996870404594858628366109486f64,
        0.983668123279747209970032581605663f64,
        0.973116322501126268374693868423707f64,
        0.960021864968307512216871025581798f64,
        0.944374444748559979415831324037439f64,
        0.926200047429274325879324277080474f64,
        0.905573307699907798546522558925958f64,
        0.882560535792052681543116462530226f64,
        0.857205233546061098958658510658944f64,
        0.829565762382768397442898119732502f64,
        0.799727835821839083013668942322683f64,
        0.767777432104826194917977340974503f64,
        0.733790062453226804726171131369528f64,
        0.697850494793315796932292388026640f64,
        0.660061064126626961370053668149271f64,
        0.620526182989242861140477556431189f64,
        0.579345235826361691756024932172540f64,
        0.536624148142019899264169793311073f64,
        0.492480467861778574993693061207709f64,
        0.447033769538089176780609900322854f64,
        0.400401254830394392535476211542661f64,
        0.352704725530878113471037207089374f64,
        0.304073202273625077372677107199257f64,
        0.254636926167889846439805129817805f64,
        0.204525116682309891438957671002025f64,
        0.153869913608583546963794672743256f64,
        0.102806937966737030147096751318001f64,
        0.051471842555317695833025213166723f64,
        0.000000000000000000000000000000000f64,
    ];

    // xgk[1], xgk[3], ... abscissae of the 30-point gauss rule.
    // xgk[0], xgk[2], ... abscissae to optimally extend the 30-point gauss rule

    // weights of the 30-point gauss rule
    let wg: [f64; 15] = [
        0.007968192496166605615465883474674f64,
        0.018466468311090959142302131912047f64,
        0.028784707883323369349719179611292f64,
        0.038799192569627049596801936446348f64,
        0.048402672830594052902938140422808f64,
        0.057493156217619066481721689402056f64,
        0.065974229882180495128128515115962f64,
        0.073755974737705206268243850022191f64,
        0.080755895229420215354694938460530f64,
        0.086899787201082979802387530715126f64,
        0.092122522237786128717632707087619f64,
        0.096368737174644259639468626351810f64,
        0.099593420586795267062780282103569f64,
        0.101762389748405504596428952168554f64,
        0.102852652893558840341285636705415f64,
    ];

    // weights of the 61-point kronrod rule
    let wgk: [f64; 31] = [
        0.001389013698677007624551591226760f64,
        0.003890461127099884051267201844516f64,
        0.006630703915931292173319826369750f64,
        0.009273279659517763428441146892024f64,
        0.011823015253496341742232898853251f64,
        0.014369729507045804812451432443580f64,
        0.016920889189053272627572289420322f64,
        0.019414141193942381173408951050128f64,
        0.021828035821609192297167485738339f64,
        0.024191162078080601365686370725232f64,
        0.026509954882333101610601709335075f64,
        0.028754048765041292843978785354334f64,
        0.030907257562387762472884252943092f64,
        0.032981447057483726031814191016854f64,
        0.034979338028060024137499670731468f64,
        0.036882364651821229223911065617136f64,
        0.038678945624727592950348651532281f64,
        0.040374538951535959111995279752468f64,
        0.041969810215164246147147541285970f64,
        0.043452539701356069316831728117073f64,
        0.044814800133162663192355551616723f64,
        0.046059238271006988116271735559374f64,
        0.047185546569299153945261478181099f64,
        0.048185861757087129140779492298305f64,
        0.049055434555029778887528165367238f64,
        0.049795683427074206357811569379942f64,
        0.050405921402782346840893085653585f64,
        0.050881795898749606492297473049805f64,
        0.051221547849258772170656282604944f64,
        0.051426128537459025933862879215781f64,
        0.051494729429451567558340433647099f64,
    ];

    let mut fv1: [f64; 31] = [0f64; 31];
    let mut fv2: [f64; 31] = [0f64; 31];

    qk(
        &xgk, &wg, &wgk, &mut fv1, &mut fv2, f, arg, a, b, result, abserr, resabs, resasc,
    );
}

pub fn qk<T>(
    xgk: &[f64],
    wg: &[f64],
    wgk: &[f64],
    fv1: &mut [f64],
    fv2: &mut [f64],
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    b: f64,
    result: &mut f64,
    abserr: &mut f64,
    resabs: &mut f64,
    resasc: &mut f64,
) {
    let n = fv1.len();

    let center = 0.5f64 * (a + b);
    let half_length = 0.5f64 * (b - a);
    let abs_half_length = unsafe { half_length.abs() };
    let f_center = f(center, arg);

    let mut result_gauss = 0f64;
    let mut result_kronrod = f_center * wgk[n - 1];

    let mut result_abs = unsafe { result_kronrod.abs() };

    if n % 2 == 0 {
        result_gauss = f_center * wg[n / 2 - 1];
    }

    for j in 0usize..((n - 1) / 2) {
        // in original fortran j=1,2,3 jtw=2,4,6
        let jtw = j * 2 + 1;
        let abscissa = half_length * xgk[jtw];
        let fval1 = f(center - abscissa, arg);
        let fval2 = f(center + abscissa, arg);
        let fsum = fval1 + fval2;

        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        result_gauss += wg[j] * fsum;
        result_kronrod += wgk[jtw] * fsum;
        result_abs += unsafe { wgk[jtw] * (fval1.abs() + fval2.abs()) };
    }

    for j in 0usize..(n / 2) {
        let jtwm1 = j * 2;
        let abscissa = half_length * xgk[jtwm1];
        let fval1 = f(center - abscissa, arg);
        let fval2 = f(center + abscissa, arg);

        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        result_kronrod += wgk[jtwm1] * (fval1 + fval2);
        result_abs += unsafe { wgk[jtwm1] * (fval1.abs() + fval2.abs()) };
    }

    let mean = result_kronrod * 0.5;

    let mut result_asc = unsafe { wgk[n - 1] * (f_center - mean).abs() };

    for j in 0usize..(n - 1) {
        result_asc += unsafe { wgk[j] * ((fv1[j] - mean).abs() + (fv2[j] - mean).abs()) };
    }

    // scale by the width of the integration region
    let err = (result_kronrod - result_gauss) * half_length;

    result_kronrod *= half_length;
    result_abs *= abs_half_length;
    result_asc *= abs_half_length;

    *result = result_kronrod;
    *resabs = result_abs;
    *resasc = result_asc;
    *abserr = rescale_error(err, result_abs, result_asc);
}

/// This function attempts to compute a Fourier integral of the function f over the semi-infinite interval [a,+\infty).
///
/// I = \int_a^{+\infty} dx f(x) sin(omega x)
/// I = \int_a^{+\infty} dx f(x) cos(omega x)
///
/// The parameter \omega and choice of \sin or \cos is taken from the table wf (the length L can take any value, since it is overridden by
/// this function to a value appropriate for the Fourier integration). The integral is computed using the QAWO algorithm over each of the
/// subintervals,
///
/// C_1 = [a, a + c]
/// C_2 = [a + c, a + 2 c]
/// ... = ...
/// C_k = [a + (k-1) c, a + k c]
///
/// where c = (2 floor(|\omega|) + 1) \pi/|\omega|. The width c is chosen to cover an odd number of periods so that the contributions from
/// the intervals alternate in sign and are monotonically decreasing when f is positive and monotonically decreasing. The sum of this sequence
/// of contributions is accelerated using the epsilon-algorithm.
///
/// This function works to an overall absolute tolerance of abserr. The following strategy is used: on each interval C_k the algorithm tries
/// to achieve the tolerance
///
/// TOL_k = u_k abserr
///
/// where u_k = (1 - p)p^{k-1} and p = 9/10. The sum of the geometric series of contributions from each interval gives an overall tolerance
/// of abserr.
///
/// If the integration of a subinterval leads to difficulties then the accuracy requirement for subsequent intervals is relaxed,
///
/// TOL_k = u_k max(abserr, max_{i<k}{E_i})
///
/// where E_k is the estimated error on the interval C_k.
///
/// The subintervals and their results are stored in the memory provided by workspace. The maximum number of subintervals is given by limit,
/// which may not exceed the allocated size of the workspace. The integration over each subinterval uses the memory provided by cycle_workspace
/// as workspace for the QAWO algorithm.
pub fn qawf<T>(
    f: ::function<T>,
    arg: &mut T,
    a: f64,
    epsabs: f64,
    limit: u64,
    workspace: &mut ::IntegrationWorkspace,
    cycle_workspace: &mut ::IntegrationWorkspace,
    wf: &mut ::IntegrationQawoTable,
    result: &mut f64,
    abserr: &mut f64,
) -> enums::Value {
    let mut total_error = 0f64;

    let mut ktmin: u64 = 0;
    let mut iteration: u64 = 0;

    unsafe {
        let mut table: ffi::extrapolation_table = ::std::mem::zeroed();

        let omega = (*ffi::FFI::unwrap_shared(wf)).omega;

        let p = 0.9f64;
        let mut factor = 1f64;
        let mut error_type = 0i32;

        /* Initialize results */
        workspace.initialise(a, a);

        *result = 0f64;
        *abserr = 0f64;

        if limit > (*ffi::FFI::unwrap_shared(workspace)).limit {
            rgsl_error!(
                "iteration limit exceeds available workspace",
                ::Value::Invalid
            );
        }

        /* Test on accuracy */
        if epsabs <= 0f64 {
            rgsl_error!(
                "absolute tolerance epsabs must be positive",
                ::Value::BadTolerance
            );
        }

        if omega == 0f64 {
            if (*ffi::FFI::unwrap_shared(wf)).sine == ::IntegrationQawo::Sine.into() {
                /* The function sin(w x) f(x) is always zero for w = 0 */
                *result = 0f64;
                *abserr = 0f64;

                return ::Value::Success;
            }
            /* The function cos(w x) f(x) is always f(x) for w = 0 */
            let limit = (*ffi::FFI::unwrap_shared(cycle_workspace)).limit;
            return cycle_workspace.qagiu(f, arg, a, epsabs, 0f64, limit as _, result, abserr);
        }

        let mut eps = if epsabs > ::DBL_MIN / (1f64 - p) {
            epsabs * (1f64 - p)
        } else {
            epsabs
        };

        let initial_eps = eps;

        let mut area = 0f64;
        let mut errsum = 0f64;

        let mut res_ext = 0f64;
        let mut err_ext = ::DBL_MAX;
        let mut correc = 0f64;

        let cycle = (2f64 * omega.abs().floor() + 1f64) * ::std::f64::consts::PI / omega.abs();

        wf.set_length(cycle);

        ::types::integration::initialise_table(&mut table);

        while iteration < limit {
            let mut area1 = 0f64;
            let mut error1 = 0f64;
            let mut reseps = 0f64;
            let mut erreps = 0f64;

            let a1 = a + iteration as f64 * cycle;
            let b1 = a1 + cycle;

            let epsabs1 = eps * factor;

            let status = wf.qawo(
                f,
                arg,
                a1,
                epsabs1,
                0f64,
                limit as _,
                cycle_workspace,
                &mut area1,
                &mut error1,
            );

            ::types::integration::append_interval(workspace, a1, b1, area1, error1);

            factor *= p;

            area = area + area1;
            errsum = errsum + error1;

            /* estimate the truncation error as 50 times the final term */
            let truncation_error = 50f64 * area1.abs();

            total_error = errsum + truncation_error;

            if total_error < epsabs && iteration > 4 {
                *result = area;
                *abserr = total_error;
                return ::types::integration::return_error(error_type);
            }

            if error1 > correc {
                correc = error1;
            }

            if status != ::Value::Success {
                eps = initial_eps.max(correc * (1f64 - p));
            }

            if status != ::Value::Success && total_error < 10f64 * correc && iteration > 3 {
                *result = area;
                *abserr = total_error;
                return ::types::integration::return_error(error_type);
            }

            ::types::integration::append_table(&mut table, area);

            if table.n < 2 {
                continue;
            }

            ::types::integration::intern_qelg(&mut table, &mut reseps, &mut erreps);

            ktmin += 1;

            if ktmin >= 15 && err_ext < 0.001f64 * total_error {
                error_type = 4;
            }

            if erreps < err_ext {
                ktmin = 0;
                err_ext = erreps;
                res_ext = reseps;

                if err_ext + 10f64 * correc <= epsabs {
                    break;
                }
                if err_ext <= epsabs && 10f64 * correc >= epsabs {
                    break;
                }
            }
            iteration += 1;
        }

        if iteration == limit {
            error_type = 1;
        }

        if err_ext == ::DBL_MAX {
            *result = area;
            *abserr = total_error;
            return ::types::integration::return_error(error_type);
        }

        err_ext = err_ext + 10f64 * correc;

        *result = res_ext;
        *abserr = err_ext;

        if error_type == 0 {
            return ::Value::Success;
        }

        if res_ext != 0f64 && area != 0f64 {
            if err_ext / res_ext.abs() > errsum / area.abs() {
                *result = area;
                *abserr = total_error;
                return ::types::integration::return_error(error_type);
            }
        } else if err_ext > errsum {
            *result = area;
            *abserr = total_error;
            return ::types::integration::return_error(error_type);
        } else if area == 0f64 {
            return ::types::integration::return_error(error_type);
        }

        /*if error_type == 4 {
            err_ext = err_ext + truncation_error;
        }*/

        ::types::integration::return_error(error_type)
    }
}
