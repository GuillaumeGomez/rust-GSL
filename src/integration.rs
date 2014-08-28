//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
The QNG algorithm is a non-adaptive procedure which uses fixed Gauss-Kronrod-Patterson abscissae to sample the integrand at a maximum of 87 
points. It is provided for fast integration of smooth functions.
!*/

use ffi;
use enums;
use std::intrinsics::{fabsf64, powf64};

fn rescale_error(err: f64, result_abs: f64, result_asc: f64) -> f64 {
    let mut t_err = unsafe { fabsf64(err) };

    if result_asc != 0f64 && t_err != 0f64 {
        let scale = unsafe { powf64((200f64 * t_err / result_asc), 1.5f64) };
        
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
pub fn qng<T>(f: ::function<T>, arg: &mut T, a: f64, b: f64, eps_abs: f64, eps_rel: f64, result: &mut f64, abs_err: &mut f64,
    n_eval: &mut u64) -> enums::Value {
    let half_length = 0.5f64 * (b - a);
    let abs_half_length = unsafe { fabsf64(half_length) };
    let center = 0.5f64 * (b + a);
    let f_center = f(center, arg);

    if eps_abs <= 0f64 && (eps_rel < 50f64 * ::DBL_EPSILON || eps_rel < 0.5e-28f64) {
        *result = 0f64;
        *abs_err = 0f64;
        *n_eval = 0u64;
        let file = file!();
        "tolerance cannot be achieved with given eps_abs and eps_rel".with_c_str(|c_str|{
            file.with_c_str(|c_file|{
                unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::BadTol as i32) }
            });
        });
        return enums::BadTol;
    }

    // w21b, weights of the 21-point formula for abscissae x2
    let w21b : [f64, ..6] = [
        0.011694638867371874278064396062192f64,
        0.054755896574351996031381300244580f64,
        0.093125454583697605535065465083366f64,
        0.123491976262065851077958109831074f64,
        0.142775938577060080797094273138717f64,
        0.149445554002916905664936468389821f64
    ];
    // w10, weights of the 10-point formula
    let w10 : [f64, ..5] = [
        0.066671344308688137593568809893332f64,
        0.149451349150580593145776339657697f64,
        0.219086362515982043995534934228163f64,
        0.269266719309996355091226921569469f64,
        0.295524224714752870173892994651338f64
    ];
    // w21a, weights of the 21-point formula for abscissae x1 */
    let w21a : [f64, ..5] = [
        0.032558162307964727478818972459390f64,
        0.075039674810919952767043140916190f64,
        0.109387158802297641899210590325805f64,
        0.134709217311473325928054001771707f64,
        0.147739104901338491374841515972068f64
    ];
    // x1, abscissae common to the 10-, 21-, 43- and 87-point rule */
    let x1 : [f64, ..5] = [
        0.973906528517171720077964012084452f64,
        0.865063366688984510732096688423493f64,
        0.679409568299024406234327365114874f64,
        0.433395394129247190799265943165784f64,
        0.148874338981631210884826001129720f64
    ];
    // x2, abscissae common to the 21-, 43- and 87-point rule
    let x2 : [f64, ..5] = [
        0.995657163025808080735527280689003f64,
        0.930157491355708226001207180059508f64,
        0.780817726586416897063717578345042f64,
        0.562757134668604683339000099272694f64,
        0.294392862701460198131126603103866f64
    ];
    // x3, abscissae common to the 43- and 87-point rule */
    let x3 : [f64, ..11] = [
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
        0.074650617461383322043914435796506f64
    ];
    // x4, abscissae of the 87-point rule */
    let x4 : [f64, ..22] = [
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
        0.037352123394619870814998165437704f64
    ];
    // w43a, weights of the 43-point formula for abscissae x1, x3 */
    let w43a : [f64, ..10] = [
        0.016296734289666564924281974617663f64,
        0.037522876120869501461613795898115f64,
        0.054694902058255442147212685465005f64,
        0.067355414609478086075553166302174f64,
        0.073870199632393953432140695251367f64,
        0.005768556059769796184184327908655f64,
        0.027371890593248842081276069289151f64,
        0.046560826910428830743339154433824f64,
        0.061744995201442564496240336030883f64,
        0.071387267268693397768559114425516f64
    ];
    // w43b, weights of the 43-point formula for abscissae x3 */
    let w43b : [f64, ..12] = [
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
        0.074722147517403005594425168280423f64
    ];
    // w87a, weights of the 87-point formula for abscissae x1, x2, x3
    let w87a : [f64, ..21] = [
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
        0.037253875503047708539592001191226f64
    ];
    // w87b, weights of the 87-point formula for abscissae x4
    let w87b : [f64, ..23] = [
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
        0.037361073762679023410321241766599f64
    ];

    // Compute the integral using the 10- and 21-point formula.
    let mut res10 = 0f64;
    let mut res21 = w21b[5] * f_center;
    let mut resabs = unsafe { w21b[5] * fabsf64(f_center) };
    let mut savfun : [f64, ..21] = [0f64, ..21];
    let mut fv1 : [f64, ..5] = [0f64, ..5];
    let mut fv2 : [f64, ..5] = [0f64, ..5];
    let mut fv3 : [f64, ..5] = [0f64, ..5];
    let mut fv4 : [f64, ..5] = [0f64, ..5];

    for k in range(0u, 5u) {
        let abscissa = half_length * x1[k];
        let fval1 = f(center + abscissa, arg);
        let fval2 = f(center - abscissa, arg);
        let fval = fval1 + fval2;

        res10 += w10[k] * fval;
        res21 += w21a[k] * fval;
        resabs += unsafe { w21a[k] * (fabsf64(fval1) + fabsf64(fval2)) };
        savfun[k] = fval;
        fv1[k] = fval1;
        fv2[k] = fval2;
    }

    for k in range(0u, 5u) {
        let abscissa = half_length * x2[k];
        let fval1 = f(center + abscissa, arg);
        let fval2 = f(center - abscissa, arg);
        let fval = fval1 + fval2;

        res21 += w21b[k] * fval;
        resabs += unsafe { w21b[k] * (fabsf64(fval1) + fabsf64(fval2)) };
        savfun[k + 5] = fval;
        fv3[k] = fval1;
        fv4[k] = fval2;
    }
    resabs *= abs_half_length;

    let mean = 0.5f64 * res21;
    let mut resasc = unsafe { w21b[5] * fabsf64(f_center - mean) };

    for k in range(0u, 5u) {
        resasc += unsafe { (w21a[k] * (fabsf64(fv1[k] - mean) + fabsf64(fv2[k] - mean)) + w21b[k] * (fabsf64(fv3[k] - mean) + fabsf64(fv4[k] - mean))) };
    }
    resasc *= abs_half_length;
    let mut result_kronrod = res21 * half_length;
    let mut err = rescale_error((res21 - res10) * half_length, resabs, resasc);

    // test for convergence.
    if err < eps_abs || err < eps_rel * unsafe { fabsf64(result_kronrod) } {
        *result = result_kronrod;
        *abs_err = err;
        *n_eval = 21;
        return enums::Success;
    }

    // compute the integral using the 43-point formula.
    let mut res43 = w43b[11] * f_center;

    for k in range(0u, 10u) {
        res43 += savfun[k] * w43a[k];
    }

    for k in range(0u, 11u) {
        let abscissa = half_length * x3[k];
        let fval = f(center + abscissa, arg) + f(center - abscissa, arg);
      
        res43 += fval * w43b[k];
        savfun[k + 10] = fval;
    }

    // test for convergence
    result_kronrod = res43 * half_length;
    err = rescale_error((res43 - res21) * half_length, resabs, resasc);

    if err < eps_abs || err < eps_rel * unsafe { fabsf64(result_kronrod) } {
        *result = result_kronrod;
        *abs_err = err;
        *n_eval = 43;
        return enums::Success;
    }

    // compute the integral using the 87-point formula.
    let mut res87 = w87b[22] * f_center;

    for k in range(0, 21) {
        res87 += savfun[k] * w87a[k];
    }
    for k in range(0, 22) {
        let abscissa = half_length * x4[k];

        res87 += w87b[k] * f(center + abscissa, arg) + f(center - abscissa, arg);
    }

    // test for convergence
    result_kronrod = res87 * half_length ;
    err = rescale_error((res87 - res43) * half_length, resabs, resasc);
  
    if err < eps_abs || err < eps_rel * unsafe { fabsf64(result_kronrod) } {
        *result = result_kronrod;
        *abs_err = err;
        *n_eval = 87;
        return enums::Success;
    }

    // failed to converge
    *result = result_kronrod;
    *abs_err = err;
    *n_eval = 87;
    let file = file!();
    "failed to reach tolerance with highest-order rule".with_c_str(|c_str|{
        file.with_c_str(|c_file|{
            unsafe { ffi::gsl_error(c_str, c_file, line!() as i32, enums::Tol as i32) }
        });
    });
    enums::Tol
}