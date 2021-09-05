//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# Polynomials

This chapter describes functions for evaluating and solving polynomials. There are routines for finding real and complex roots of quadratic
and cubic equations using analytic methods. An iterative polynomial solver is also available for finding the roots of general polynomials
with real coefficients (of any order).

## References and Further Reading

The balanced-QR method and its error analysis are described in the following papers,

R.S. Martin, G. Peters and J.H. Wilkinson, “The QR Algorithm for Real Hessenberg Matrices”, Numerische Mathematik, 14 (1970), 219–231.
B.N. Parlett and C. Reinsch, “Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors”, Numerische Mathematik, 13 (1969), 293–304.
A. Edelman and H. Murakami, “Polynomial roots from companion matrix eigenvalues”, Mathematics of Computation, Vol. 64, No. 210 (1995), 763–776.
The formulas for divided differences are given in the following texts,

Abramowitz and Stegun, Handbook of Mathematical Functions, Sections 25.1.4 and 25.2.26.
R. L. Burden and J. D. Faires, Numerical Analysis, 9th edition, ISBN 0-538-73351-9, 2011.
!*/

/// The functions described here evaluate the polynomial
/// `P(x) = c[0] + c[1] x + c[2] x^2 + \dots + c[len-1] x^{len-1}` using Horner’s method for
/// stability.
pub mod evaluation {
    use crate::Value;
    use std::mem::transmute;
    use types::ComplexF64;

    /// This function evaluates a polynomial with real coefficients for the real variable x.
    #[doc(alias = "gsl_poly_eval")]
    pub fn poly_eval(c: &[f64], x: f64) -> f64 {
        unsafe { sys::gsl_poly_eval(c.as_ptr(), c.len() as i32, x) }
    }

    /// This function evaluates a polynomial with real coefficients for the complex variable z.
    #[doc(alias = "gsl_poly_complex_eval")]
    pub fn poly_complex_eval(c: &[f64], z: &ComplexF64) -> ComplexF64 {
        unsafe {
            transmute(sys::gsl_poly_complex_eval(
                c.as_ptr(),
                c.len() as i32,
                transmute(*z),
            ))
        }
    }

    /// This function evaluates a polynomial with complex coefficients for the complex variable z.
    #[doc(alias = "gsl_complex_poly_complex_eval")]
    pub fn complex_poly_complex_eval(c: &[ComplexF64], z: &ComplexF64) -> ComplexF64 {
        let mut tmp = Vec::new();

        for it in c.iter() {
            unsafe { tmp.push(transmute(*it)) };
        }
        unsafe {
            transmute(sys::gsl_complex_poly_complex_eval(
                tmp.as_ptr(),
                tmp.len() as i32,
                transmute(*z),
            ))
        }
    }

    /// This function evaluates a polynomial and its derivatives storing the results in the array res of size lenres. The output array contains
    /// the values of d^k P/d x^k for the specified value of x starting with k = 0.
    #[doc(alias = "gsl_poly_eval_derivs")]
    pub fn poly_eval_derivs(c: &[f64], x: f64, res: &mut [f64]) -> Value {
        Value::from(unsafe {
            sys::gsl_poly_eval_derivs(
                c.as_ptr(),
                c.len() as _,
                x,
                res.as_mut_ptr(),
                res.len() as _,
            )
        })
    }
}

/// The functions described here manipulate polynomials stored in Newton’s divided-difference representation. The use of divided-differences
/// is described in Abramowitz & Stegun sections 25.1.4 and 25.2.26, and Burden and Faires, chapter 3, and discussed briefly below.
///
/// Given a function f(x), an nth degree interpolating polynomial P_{n}(x) can be constructed which agrees with f at n+1 distinct points x_0,
/// x_1,...,x_{n}. This polynomial can be written in a form known as Newton’s divided-difference representation:
///
/// P_n(x) = f(x_0) + \sum_(k=1)^n [x_0,x_1,...,x_k] (x-x_0)(x-x_1)...(x-x_(k-1))
///
/// where the divided differences [x_0,x_1,...,x_k] are defined in section 25.1.4 of Abramowitz and Stegun. Additionally, it is possible to
/// construct an interpolating polynomial of degree 2n+1 which also matches the first derivatives of f at the points x_0,x_1,...,x_n. This is
/// called the Hermite interpolating polynomial and is defined as
///
/// H_(2n+1)(x) = f(z_0) + \sum_(k=1)^(2n+1) [z_0,z_1,...,z_k] (x-z_0)(x-z_1)...(x-z_(k-1))
///
/// where the elements of z = \{x_0,x_0,x_1,x_1,...,x_n,x_n\} are defined by z_{2k} = z_{2k+1} = x_k. The divided-differences [z_0,z_1,...,z_k]
/// are discussed in Burden and Faires, section 3.4.
pub mod divided_difference_representation {
    use crate::Value;

    /// This function computes a divided-difference representation of the interpolating polynomial for the points (x, y) stored in the arrays
    /// xa and ya of length size. On output the divided-differences of (xa,ya) are stored in the array dd, also of length size. Using the
    /// notation above, `dd[k] = [x_0,x_1,...,x_k]`.
    #[doc(alias = "gsl_poly_dd_init")]
    pub fn poly_dd_init(dd: &mut [f64], xa: &[f64], ya: &[f64]) -> Value {
        Value::from(unsafe {
            sys::gsl_poly_dd_init(dd.as_mut_ptr(), xa.as_ptr(), ya.as_ptr(), dd.len() as _)
        })
    }

    /// This function evaluates the polynomial stored in divided-difference form in the arrays dd and xa of length size at the point x.
    #[doc(alias = "gsl_poly_dd_eval")]
    pub fn poly_dd_eval(dd: &[f64], xa: &[f64], x: f64) -> f64 {
        unsafe { sys::gsl_poly_dd_eval(dd.as_ptr(), xa.as_ptr(), dd.len() as _, x) }
    }

    /// This function converts the divided-difference representation of a polynomial to a Taylor expansion. The divided-difference representation
    /// is supplied in the arrays dd and xa of length size. On output the Taylor coefficients of the polynomial expanded about the point xp are
    /// stored in the array c also of length size. A workspace of length size must be provided in the array w.
    #[doc(alias = "gsl_poly_dd_taylor")]
    pub fn poly_dd_taylor(c: &mut [f64], xp: f64, dd: &[f64], xa: &[f64], w: &mut [f64]) -> Value {
        Value::from(unsafe {
            sys::gsl_poly_dd_taylor(
                c.as_mut_ptr(),
                xp,
                dd.as_ptr(),
                xa.as_ptr(),
                dd.len() as _,
                w.as_mut_ptr(),
            )
        })
    }

    /// This function computes a divided-difference representation of the interpolating Hermite polynomial for the points (x, y) stored in the
    /// arrays xa and ya of length size. Hermite interpolation constructs polynomials which also match first derivatives dy/dx which are provided
    /// in the array dya also of length size. The first derivatives can be incorported into the usual divided-difference algorithm by forming a
    /// new dataset z = \{x_0,x_0,x_1,x_1,...\}, which is stored in the array za of length 2*size on output. On output the divided-differences
    /// of the Hermite representation are stored in the array dd, also of length 2*size. Using the notation above, `dd[k] = [z_0,z_1,...,z_k]`.
    /// The resulting Hermite polynomial can be evaluated by calling gsl_poly_dd_eval and using za for the input argument xa.
    #[doc(alias = "gsl_poly_dd_hermite_init")]
    pub fn poly_dd_hermite_init(
        dd: &mut [f64],
        za: &mut [f64],
        xa: &[f64],
        ya: &[f64],
        dya: &[f64],
    ) -> Value {
        Value::from(unsafe {
            sys::gsl_poly_dd_hermite_init(
                dd.as_mut_ptr(),
                za.as_mut_ptr(),
                xa.as_ptr(),
                ya.as_ptr(),
                dya.as_ptr(),
                dd.len() as _,
            )
        })
    }
}

pub mod quadratic_equations {
    use std::mem::transmute;
    use types::ComplexF64;

    /// This function finds the real roots of the quadratic equation,
    ///
    /// a x^2 + b x + c = 0
    ///
    /// The number of real roots (either zero, one or two) is returned, and their locations are
    /// stored in x0 and x1. If no real roots are found then x0 and x1 are not modified. If one real
    /// root is found (i.e. if a=0) then it is stored in x0. When two real roots are found they
    /// are stored in x0 and x1 in ascending order. The case of coincident roots is not considered
    /// special. For example (x-1)^2=0 will have two roots, which happen to have exactly equal
    /// values.
    ///
    /// The number of roots found depends on the sign of the discriminant b^2 - 4 a c. This will be
    /// subject to rounding and cancellation errors when computed in double precision, and will also
    /// be subject to errors if the coefficients of the polynomial are inexact. These errors
    /// may cause a discrete change in the number of roots. However, for polynomials with small
    /// integer coefficients the discriminant can always be computed exactly.
    ///
    /// Returns `(Value, x0, x1)`.
    #[doc(alias = "gsl_poly_solve_quadratic")]
    pub fn poly_solve_quadratic(a: f64, b: f64, c: f64) -> (::Value, f64, f64) {
        let mut x0 = 0.;
        let mut x1 = 0.;
        let ret = unsafe { sys::gsl_poly_solve_quadratic(a, b, c, &mut x0, &mut x1) };
        (::Value::from(ret), x0, x1)
    }

    /// This function finds the complex roots of the quadratic equation,
    ///
    /// a z^2 + b z + c = 0
    ///
    /// The number of complex roots is returned (either one or two) and the locations of the roots are stored in z0 and z1. The roots are returned
    /// in ascending order, sorted first by their real components and then by their imaginary components. If only one real root is found (i.e. if
    /// a=0) then it is stored in z0.
    #[doc(alias = "gsl_poly_complex_solve_quadratic")]
    pub fn poly_complex_solve_quadratic(
        a: f64,
        b: f64,
        c: f64,
        z0: &mut ComplexF64,
        z1: &mut ComplexF64,
    ) -> ::Value {
        ::Value::from(unsafe {
            sys::gsl_poly_complex_solve_quadratic(a, b, c, transmute(z0), transmute(z1))
        })
    }
}

pub mod cubic_equations {
    use std::mem::transmute;
    use types::ComplexF64;

    /// This function finds the real roots of the cubic equation,
    ///
    /// x^3 + a x^2 + b x + c = 0
    ///
    /// with a leading coefficient of unity. The number of real roots (either one or three) is
    /// returned, and their locations are stored in x0, x1 and x2. If one real root is found then
    /// only x0 is modified. When three real roots are found they are stored in x0, x1 and x2 in
    /// ascending order. The case of coincident roots is not considered special. For example, the
    /// equation (x-1)^3=0 will have three roots with exactly equal values. As in the quadratic
    /// case, finite precision may cause equal or closely-spaced real roots to move off the
    /// real axis into the complex plane, leading to a discrete change in the number of real roots.
    ///
    /// Returns `(Value, x0, x1, x2)`.
    #[doc(alias = "gsl_poly_solve_cubic")]
    pub fn poly_solve_cubic(a: f64, b: f64, c: f64) -> (::Value, f64, f64, f64) {
        let mut x0 = 0.;
        let mut x1 = 0.;
        let mut x2 = 0.;
        let ret = unsafe { sys::gsl_poly_solve_cubic(a, b, c, &mut x0, &mut x1, &mut x2) };
        (::Value::from(ret), x0, x1, x2)
    }

    /// This function finds the complex roots of the cubic equation,
    ///
    /// z^3 + a z^2 + b z + c = 0
    ///
    /// The number of complex roots is returned (always three) and the locations of the roots are
    /// stored in z0, z1 and z2. The roots are returned in ascending order, sorted first by their
    /// real components and then by their imaginary components.
    #[doc(alias = "gsl_poly_complex_solve_cubic")]
    pub fn poly_complex_solve_cubic(
        a: f64,
        b: f64,
        c: f64,
        z0: &mut ComplexF64,
        z1: &mut ComplexF64,
        z2: &mut ComplexF64,
    ) -> ::Value {
        ::Value::from(unsafe {
            sys::gsl_poly_complex_solve_cubic(a, b, c, transmute(z0), transmute(z1), transmute(z2))
        })
    }
}
