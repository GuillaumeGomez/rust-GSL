//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

pub fn subinterval_too_small(a1: f64, a2: f64, b2: f64) -> bool {
    let e = crate::DBL_EPSILON;
    let u = crate::DBL_MIN;

    let tmp = (1f64 + 100f64 * e) * (a2.abs() + 1000f64 * u);

    a1.abs() <= tmp && b2.abs() <= tmp
}
