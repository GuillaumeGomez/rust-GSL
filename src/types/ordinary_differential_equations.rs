//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*! Numerical ODE solvers.

#Ordinary Differential Equations

This chapter describes functions for solving ordinary differential equation (ODE) initial value problems. The library provides a variety of
low-level methods, such as Runge-Kutta and Bulirsch-Stoer routines, and higher-level components for adaptive step-size control. The
components can be combined by the user to achieve the desired solution, with full access to any intermediate steps. A driver object can be
used as a high level wrapper for easy use of low level functions.

##Defining the ODE System

The routines solve the general n-dimensional first-order system,

dy_i(t)/dt = f_i(t, y_1(t), ..., y_n(t))
for i = 1, \dots, n. The stepping functions rely on the vector of derivatives f_i and the Jacobian matrix, J_{ij} = df_i(t,y(t)) / dy_j. A
system of equations is defined using the gsl_odeiv2_system datatype.

##Stepping Functions

The lowest level components are the stepping functions which advance a solution from time t to t+h for a fixed step-size h and estimate
the resulting local error.

##Adaptive Step-size Control

The control function examines the proposed change to the solution produced by a stepping function and attempts to determine the optimal
step-size for a user-specified level of error.

##Evolution

The evolution function combines the results of a stepping function and control function to reliably advance the solution forward one
step using an acceptable step-size.

##Driver

The driver object is a high level wrapper that combines the evolution, control and stepper objects for easy use.

##References and Further Reading

Ascher, U.M., Petzold, L.R., Computer Methods for Ordinary Differential and Differential-Algebraic Equations, SIAM, Philadelphia, 1998.
Hairer, E., Norsett, S. P., Wanner, G., Solving Ordinary Differential Equations I: Nonstiff Problems, Springer, Berlin, 1993.
Hairer, E., Wanner, G., Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems, Springer, Berlin, 1996.
Many of the basic Runge-Kutta formulas can be found in the Handbook of Mathematical Functions,

Abramowitz & Stegun (eds.), Handbook of Mathematical Functions, Section 25.5.
The implicit Bulirsch-Stoer algorithm bsimp is described in the following paper,

G. Bader and P. Deuflhard, “A Semi-Implicit Mid-Point Rule for Stiff Systems of Ordinary Differential Equations.”, Numer. Math. 41,
373–398, 1983.
The Adams and BDF multistep methods msadams and msbdf are based on the following articles,

G. D. Byrne and A. C. Hindmarsh, “A Polyalgorithm for the Numerical Solution of Ordinary Differential Equations.”, ACM Trans. Math.
Software, 1, 71–96, 1975.
P. N. Brown, G. D. Byrne and A. C. Hindmarsh, “VODE: A Variable-coefficient ODE Solver.”, SIAM J. Sci. Stat. Comput. 10, 1038–1051, 1989.
A. C. Hindmarsh, P. N. Brown, K. E. Grant, S. L. Lee, R. Serban, D. E. Shumaker and C. S. Woodward, “SUNDIALS: Suite of Nonlinear and
Differential/Algebraic Equation Solvers.”, ACM Trans. Math. Software 31, 363–396, 2005.
!*/

use ffi;
use enums::{self, GSLResult};
use libc::c_void;

/// Description of a system of ODEs.
///
/// `y' = f(t,y) = dydt(t, y)`
///
/// The system is specified by giving the right-hand-side of the equation and possibly a jacobian function.
///
/// Some methods require the jacobian function, which calculates the matrix dfdy and the vector dfdt. The matrix dfdy conforms
/// to the GSL standard, being a continuous range of floating point values, in row-order.
pub struct ODEiv2System<'a> {
    function: &'a mut FnMut(f64, &[f64], &mut [f64]) -> GSLResult<()>,
    jacobian: Option<&'a mut FnMut(f64, &[f64], &mut [f64], &mut [f64]) -> GSLResult<()>>,
    dimension: usize,
}

impl<'a> ODEiv2System<'a> {
    /// Returns a new ODEiv2System with a given dimension and right-hand side.
    pub fn new(dimension: usize,
               function: &'a mut FnMut(f64, &[f64], &mut [f64])
                        -> GSLResult<()>) -> ODEiv2System<'a> {
        ODEiv2System {
            function: function,
            jacobian: None,
            dimension: dimension,
        }
    }

    /// Returns a new ODEiv2System with a jacobian function provided.
    pub fn with_jacobian(dimension: usize,
               function: &'a mut FnMut(f64, &[f64], &mut [f64])
                        -> GSLResult<()>,
               jacobian: &'a mut FnMut(f64, &[f64], &mut [f64], &mut [f64])
                        -> GSLResult<()>) -> ODEiv2System<'a> {
        ODEiv2System {
            function: function,
            jacobian: Some(jacobian),
            dimension: dimension,
        }
    }

    /// Return `ffi::gsl_odeiv2_system` structure.
    fn to_raw(&mut self) -> ffi::gsl_odeiv2_system {
        ffi::gsl_odeiv2_system {
            function: function_handler,
            jacobian: if self.jacobian.is_some() { Some(jacobian_handler) } else { None },
            dimension: self.dimension,
            params: self as *mut _ as *mut c_void,
        }
    }
}

/// Default handler for calling the function closure.
extern fn function_handler(t: f64, t_y: *const f64, t_f: *mut f64, params: *mut c_void)
        -> enums::Value {
    let mut sys = unsafe { &mut *(params as *mut ODEiv2System) };
    let n = sys.dimension as usize;
    let t_y = unsafe { ::std::slice::from_raw_parts(t_y, n) };
    let t_f = unsafe { ::std::slice::from_raw_parts_mut(t_f, n) };

    match (sys.function)(t, t_y, t_f) {
        Ok(()) => enums::Value::Success,
        Err(e) => e,
    }
}

/// Default handler for calling the jacobian closure.
extern fn jacobian_handler(t: f64, t_y: *const f64, t_dfdy: *mut f64, t_dfdt: *mut f64,
                           params: *mut c_void) -> enums::Value {
    let mut sys = unsafe { &mut *(params as *mut ODEiv2System) };
    let n = sys.dimension as usize;
    let t_y = unsafe { ::std::slice::from_raw_parts(t_y, n) };
    let t_dfdy = unsafe { ::std::slice::from_raw_parts_mut(t_dfdy, n * n) };
    let t_dfdt = unsafe { ::std::slice::from_raw_parts_mut(t_dfdt, n) };

    match sys.jacobian {
        Some(ref mut j) =>
            match j(t, t_y, t_dfdy, t_dfdt) {
                Ok(()) => enums::Value::Success,
                Err(e) => e,
            },
        None => enums::Value::BadFunction,
    }
}

pub struct ODEiv2Step {
    s: *mut ffi::gsl_odeiv2_step
}

impl ODEiv2Step {
    /// This function returns a pointer to a newly allocated instance of a stepping function of type T for a system of dim dimensions.
    /// Please note that if you use a stepper method that requires access to a driver object, it is advisable to use a driver allocation
    /// method, which automatically allocates a stepper, too.
    pub fn new(t: &ODEiv2StepType, dim: usize) -> Option<ODEiv2Step> {
        let tmp = unsafe { ffi::gsl_odeiv2_step_alloc(ffi::FFI::unwrap(t), dim) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Step {
                s: tmp,
            })
        }
    }

    /// This function resets the stepping function s. It should be used whenever the next use of s will not be a continuation of a previous
    /// step.
    pub fn reset(&mut self) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_step_reset(self.s) })
    }

    /// This function returns a pointer to the name of the stepping function. For example,
    ///
    /// ```Rust
    /// println!("step method is '{}'", s.name().unwrap());
    /// ```
    /// would print something like step method is 'rkf45'.
    pub fn name(&self) -> Option<String> {
        let tmp = unsafe { ffi::gsl_odeiv2_step_name(self.s) };

        if tmp.is_null() {
            None
        } else {
            unsafe { Some(String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()) }
        }
    }

    /// This function returns the order of the stepping function on the previous step. The order can vary if the stepping function itself is
    /// adaptive.
    pub fn order(&self) -> u32 {
        unsafe { ffi::gsl_odeiv2_step_order(self.s) }
    }

    /// This function sets a pointer of the driver object d for stepper s, to allow the stepper to access control (and evolve) object through
    /// the driver object. This is a requirement for some steppers, to get the desired error level for internal iteration of stepper.
    /// Allocation of a driver object calls this function automatically.
    pub fn set_driver(&mut self, d: &ODEiv2Driver) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_step_set_driver(self.s, d.d) })
    }

    /// This function applies the stepping function s to the system of equations defined by sys, using the step-size h to advance the system
    /// from time t and state y to time t+h. The new state of the system is stored in y on output, with an estimate of the absolute error
    /// in each component stored in yerr. If the argument dydt_in is not null it should point an array containing the derivatives for the
    /// system at time t on input. This is optional as the derivatives will be computed internally if they are not provided, but allows
    /// the reuse of existing derivative information. On output the new derivatives of the system at time t+h will be stored in dydt_out
    /// if it is not null.
    ///
    /// The stepping function returns enums::value::Failure if it is unable to compute the requested step. Also, if the user-supplied functions defined
    /// in the system sys return a status other than ::Value::Success the step will be aborted. In that case, the elements of y will be restored
    /// to their pre-step values and the error code from the user-supplied function will be returned. Failure may be due to a singularity in
    /// the system or too large step-size h. In that case the step should be attempted again with a smaller step-size, e.g. h/2.
    ///
    /// If the driver object is not appropriately set via gsl_odeiv2_step_set_driver for those steppers that need it, the stepping function
    /// returns ::Fault. If the user-supplied functions defined in the system sys returns enums::value::BadFunc, the function returns
    /// immediately with the same return code. In this case the user must call gsl_odeiv2_step_reset before calling this function again.
    pub fn apply(&mut self, t: f64, h: f64, y: &mut [f64], yerr: &mut [f64], dydt_in: &[f64], dydt_out: &mut [f64],
        sys: &mut ODEiv2System) -> GSLResult<()> {
        let sys_raw = sys.to_raw();
        let r = unsafe { ffi::gsl_odeiv2_step_apply(self.s, t, h, y.as_mut_ptr(), yerr.as_mut_ptr(), dydt_in.as_ptr(), dydt_out.as_mut_ptr(),
            &sys_raw as *const ffi::gsl_odeiv2_system) };
        GSLResult::from(r)
    }
}

impl Drop for ODEiv2Step {
    fn drop(&mut self) {
        unsafe { ffi::gsl_odeiv2_step_free(self.s) };
        self.s = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_odeiv2_step> for ODEiv2Step {
    fn wrap(s: *mut ffi::gsl_odeiv2_step) -> ODEiv2Step {
        ODEiv2Step {
            s: s
        }
    }

    fn unwrap(s: &ODEiv2Step) -> *mut ffi::gsl_odeiv2_step {
        s.s
    }
}

#[derive(Clone, Copy)]
pub struct ODEiv2StepType {
    t: *const ffi::gsl_odeiv2_step_type
}

impl ODEiv2StepType {
    /// Explicit embedded Runge-Kutta (2, 3) method.
    pub fn rk2() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rk2
            }
        }
    }

    /// Explicit 4th order (classical) Runge-Kutta. Error estimation is carried out by the step doubling method. For more efficient
    /// estimate of the error, use the embedded methods described below.
    pub fn rk4() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rk4
            }
        }
    }

    /// Explicit embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good general-purpose integrator.
    pub fn rk45() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rkf45
            }
        }
    }

    /// Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.
    pub fn rkck() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rkck
            }
        }
    }

    /// Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.
    pub fn rk8pd() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rk8pd
            }
        }
    }

    /// Implicit Gaussian first order Runge-Kutta. Also known as implicit Euler or backward Euler method. Error estimation is carried out by
    /// the step doubling method. This algorithm requires the Jacobian and access to the driver object via gsl_odeiv2_step_set_driver.
    pub fn rk1imp() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rk1imp
            }
        }
    }

    /// Implicit Gaussian second order Runge-Kutta. Also known as implicit mid-point rule. Error estimation is carried out by the step doubling
    /// method. This stepper requires the Jacobian and access to the driver object via gsl_odeiv2_step_set_driver.
    pub fn rk2imp() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rk2imp
            }
        }
    }

    /// Implicit Gaussian 4th order Runge-Kutta. Error estimation is carried out by the step doubling method. This algorithm requires the
    /// Jacobian and access to the driver object via gsl_odeiv2_step_set_driver.
    pub fn rk4imp() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_rk4imp
            }
        }
    }

    /// Implicit Bulirsch-Stoer method of Bader and Deuflhard. The method is generally suitable for stiff problems. This stepper requires
    /// the Jacobian.
    pub fn bsimp() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_bsimp
            }
        }
    }

    /// A variable-coefficient linear multistep Adams method in Nordsieck form. This stepper uses explicit Adams-Bashforth (predictor) and
    /// implicit Adams-Moulton (corrector) methods in P(EC)^m functional iteration mode. Method order varies dynamically between 1 and 12.
    /// This stepper requires the access to the driver object via gsl_odeiv2_step_set_driver.
    pub fn msadams() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_msadams
            }
        }
    }

    /// A variable-coefficient linear multistep backward differentiation formula (BDF) method in Nordsieck form. This stepper uses the explicit
    /// BDF formula as predictor and implicit BDF formula as corrector. A modified Newton iteration method is used to solve the system of
    /// non-linear equations. Method order varies dynamically between 1 and 5. The method is generally suitable for stiff problems. This
    /// stepper requires the Jacobian and the access to the driver object via gsl_odeiv2_step_set_driver.
    pub fn msbdf() -> ODEiv2StepType {
        unsafe {
            ODEiv2StepType {
                t: ffi::gsl_odeiv2_step_msbdf
            }
        }
    }
}

impl ffi::FFI<ffi::gsl_odeiv2_step_type> for ODEiv2StepType {
    fn wrap(t: *mut ffi::gsl_odeiv2_step_type) -> ODEiv2StepType {
        ODEiv2StepType {
            t: t
        }
    }

    fn unwrap(t: &ODEiv2StepType) -> *mut ffi::gsl_odeiv2_step_type {
        t.t as *mut ffi::gsl_odeiv2_step_type
    }
}

pub struct ODEiv2Control {
    c: *mut ffi::gsl_odeiv2_control
}

impl ODEiv2Control {
    /// The standard control object is a four parameter heuristic based on absolute and relative errors eps_abs and eps_rel, and scaling
    /// factors a_y and a_dydt for the system state y(t) and derivatives y'(t) respectively.
    ///
    /// The step-size adjustment procedure for this method begins by computing the desired error level D_i for each component,
    ///
    /// D_i = eps_abs + eps_rel * (a_y |y_i| + a_dydt h |y\prime_i|)
    /// and comparing it with the observed error E_i = |yerr_i|. If the observed error E exceeds the desired error level D by more than
    /// 10% for any component then the method reduces the step-size by an appropriate factor,
    ///
    /// h_new = h_old * S * (E/D)^(-1/q)
    /// where q is the consistency order of the method (e.g. q=4 for 4(5) embedded RK), and S is a safety factor of 0.9. The ratio E/D is
    /// taken to be the maximum of the ratios E_i/D_i.
    ///
    /// If the observed error E is less than 50% of the desired error level D for the maximum ratio E_i/D_i then the algorithm takes the
    /// opportunity to increase the step-size to bring the error in line with the desired level,
    ///
    /// h_new = h_old * S * (E/D)^(-1/(q+1))
    ///
    /// This encompasses all the standard error scaling methods. To avoid uncontrolled changes in the stepsize, the overall scaling factor
    /// is limited to the range 1/5 to 5.
    pub fn standard_new(eps_abs: f64, eps_rel: f64, a_y: f64, a_dydt: f64) -> Option<ODEiv2Control> {
        let tmp = unsafe { ffi::gsl_odeiv2_control_standard_new(eps_abs, eps_rel, a_y, a_dydt) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Control {
                c: tmp
            })
        }
    }

    /// This function creates a new control object which will keep the local error on each step within an absolute error of eps_abs and relative
    /// error of eps_rel with respect to the solution y_i(t). This is equivalent to the standard control object with a_y=1 and a_dydt=0.
    pub fn y_new(eps_abs: f64, eps_rel: f64) -> Option<ODEiv2Control> {
        let tmp = unsafe { ffi::gsl_odeiv2_control_y_new(eps_abs, eps_rel) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Control {
                c: tmp
            })
        }
    }

    /// This function creates a new control object which will keep the local error on each step within an absolute error of eps_abs and relative
    /// error of eps_rel with respect to the derivatives of the solution y'_i(t). This is equivalent to the standard control object with
    /// a_y=0 and a_dydt=1.
    pub fn yp_new(eps_abs: f64, eps_rel: f64) -> Option<ODEiv2Control> {
        let tmp = unsafe { ffi::gsl_odeiv2_control_yp_new(eps_abs, eps_rel) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Control {
                c: tmp
            })
        }
    }

    /// This function creates a new control object which uses the same algorithm as gsl_odeiv2_control_standard_new but with an absolute error
    /// which is scaled for each component by the array scale_abs. The formula for D_i for this control object is,
    ///
    /// D_i = eps_abs * s_i + eps_rel * (a_y |y_i| + a_dydt h |y\prime_i|)
    ///
    /// where s_i is the i-th component of the array scale_abs. The same error control heuristic is used by the Matlab ODE suite.
    pub fn scaled_new(eps_abs: f64, eps_rel: f64, a_y: f64, a_dydt: f64, scale_abs: &[f64]) -> Option<ODEiv2Control> {
        let tmp = unsafe { ffi::gsl_odeiv2_control_scaled_new(eps_abs, eps_rel, a_y, a_dydt, scale_abs.as_ptr(), scale_abs.len() as usize) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Control {
                c: tmp
            })
        }
    }

    /// This function returns a pointer to a newly allocated instance of a control function of type T. This function is only needed for
    /// defining new types of control functions. For most purposes the standard control functions described above should be sufficient.
    pub fn alloc(t: &ODEiv2ControlType) -> Option<ODEiv2Control> {
        let tmp = unsafe { ffi::gsl_odeiv2_control_alloc(ffi::FFI::unwrap(t)) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Control {
                c: tmp
            })
        }
    }

    /// This function initializes the control function c with the parameters eps_abs (absolute error), eps_rel (relative error), a_y
    /// (scaling factor for y) and a_dydt (scaling factor for derivatives).
    pub fn init(&mut self, eps_abs: f64, eps_rel: f64, a_y: f64, a_dydt: f64) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_control_init(self.c, eps_abs, eps_rel, a_y, a_dydt) })
    }

    /// This function adjusts the step-size h using the control function c, and the current values of y, yerr and dydt. The stepping function
    /// step is also needed to determine the order of the method. If the error in the y-values yerr is found to be too large then the step-size
    /// h is reduced and the function returns ODEiv::Dec. If the error is sufficiently small then h may be increased and
    /// ODEiv::Inc is returned. The function returns ODEiv::Nil if the step-size is unchanged. The goal of the function is to estimate
    /// the largest step-size which satisfies the user-specified accuracy requirements for the current point.
    pub fn hadjust(&mut self, s: &ODEiv2Step, y: &[f64], yerr: &[f64], dydt: &[f64], h: &mut f64) -> ::ODEiv {
        unsafe { ffi::gsl_odeiv2_control_hadjust(self.c, s.s, y.as_ptr(), yerr.as_ptr(), dydt.as_ptr(), h) }
    }

    /// This function returns a pointer to the name of the control function. For example,
    ///
    /// ```Rust
    /// println!("control method is '{}'", c.name());
    /// ```
    /// would print something like control method is 'standard'
    pub fn name(&self) -> Option<String> {
        let tmp = unsafe { ffi::gsl_odeiv2_control_name(self.c) };

        if tmp.is_null() {
            None
        } else {
            unsafe { Some(String::from_utf8_lossy(::std::ffi::CStr::from_ptr(tmp).to_bytes()).to_string()) }
        }
    }

    /// This function calculates the desired error level of the ind-th component to errlev. It requires the value (y) and value of the derivative
    /// (dydt) of the component, and the current step size h.
    pub fn errlevel(&mut self, y: f64, dydt: f64, h: f64, ind: usize, errlev: &mut f64) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_control_errlevel(self.c, y, dydt, h, ind, errlev) })
    }

    /// This function sets a pointer of the driver object d for control object c.
    pub fn set_driver(&mut self, d: &ODEiv2Driver) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_control_set_driver(self.c, d.d) })
    }
}

impl Drop for ODEiv2Control {
    fn drop(&mut self) {
        unsafe { ffi::gsl_odeiv2_control_free(self.c) };
        self.c = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_odeiv2_control> for ODEiv2Control {
    fn wrap(c: *mut ffi::gsl_odeiv2_control) -> ODEiv2Control {
        ODEiv2Control {
            c: c
        }
    }

    fn unwrap(c: &ODEiv2Control) -> *mut ffi::gsl_odeiv2_control {
        c.c
    }
}

#[derive(Clone, Copy)]
pub struct ODEiv2ControlType {
    t: *const ffi::gsl_odeiv2_control_type
}

impl ODEiv2ControlType {
    pub fn scaled() -> ODEiv2ControlType {
        unsafe {
            ODEiv2ControlType {
                t: ffi::gsl_odeiv2_control_scaled
            }
        }
    }

    pub fn standard() -> ODEiv2ControlType {
        unsafe {
            ODEiv2ControlType {
                t: ffi::gsl_odeiv2_control_standard
            }
        }
    }
}

impl ffi::FFI<ffi::gsl_odeiv2_control_type> for ODEiv2ControlType {
    fn wrap(t: *mut ffi::gsl_odeiv2_control_type) -> ODEiv2ControlType {
        ODEiv2ControlType {
            t: t
        }
    }

    fn unwrap(t: &ODEiv2ControlType) -> *mut ffi::gsl_odeiv2_control_type {
        t.t as *mut ffi::gsl_odeiv2_control_type
    }
}

pub struct ODEiv2Evolve {
    e: *mut ffi::gsl_odeiv2_evolve
}

impl ODEiv2Evolve {
    /// This function returns a pointer to a newly allocated instance of an evolution function for a system of dim dimensions.
    pub fn new(dim: usize) -> Option<ODEiv2Evolve> {
        let tmp = unsafe { ffi::gsl_odeiv2_evolve_alloc(dim) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Evolve {
                e: tmp
            })
        }
    }

    /// This function advances the system (e, sys) from time t and position y using the stepping function step. The new time and position
    /// are stored in t and y on output.
    ///
    /// The initial step-size is taken as h. The control function con is applied to check whether the local error estimated by the stepping
    /// function step using step-size h exceeds the required error tolerance. If the error is too high, the step is retried by calling step
    /// with a decreased step-size. This process is continued until an acceptable step-size is found. An estimate of the local error for
    /// the step can be obtained from the components of the array e->yerr[].
    ///
    /// If the user-supplied functions defined in the system sys returns enums::value::BadFunc, the function returns immediately with the same
    /// return code. In this case the user must call gsl_odeiv2_step_reset and gsl_odeiv2_evolve_reset before calling this function again.
    ///
    /// Otherwise, if the user-supplied functions defined in the system sys or the stepping function step return a status other than
    /// ::Value::Success, the step is retried with a decreased step-size. If the step-size decreases below machine precision, a status of
    /// ::Failuer is returned if the user functions returned ::Value::Success. Otherwise the value returned by user function is returned.
    /// If no acceptable step can be made, t and y will be restored to their pre-step values and h contains the final attempted step-size.
    ///
    /// If the step is successful the function returns a suggested step-size for the next step in h. The maximum time t1 is guaranteed not
    /// to be exceeded by the time-step. On the final time-step the value of t will be set to t1 exactly.
    pub fn apply(&mut self, c: &ODEiv2Control, s: &ODEiv2Step, sys: &mut ODEiv2System, t: &mut f64, t1: f64, h: &mut f64,
        y: &mut [f64]) -> GSLResult<()> {
        let sys_raw = sys.to_raw();
        let psys = &sys_raw as *const _;
        GSLResult::from(unsafe { ffi::gsl_odeiv2_evolve_apply(self.e, c.c, s.s, psys, t, t1, h, y.as_mut_ptr()) })
    }

    /// This function advances the ODE-system (e, sys, con) from time t and position y using the stepping function step by a specified step
    /// size h. If the local error estimated by the stepping function exceeds the desired error level, the step is not taken and the function
    /// returns enums::value::Failure. Otherwise the value returned by user function is returned.
    pub fn apply_fixed_step(&mut self, c: &ODEiv2Control, s: &ODEiv2Step, sys: &mut ODEiv2System, t: &mut f64, h: f64,
        y: &mut [f64]) -> GSLResult<()> {
        let sys_raw = sys.to_raw();
        let psys = &sys_raw as *const _;
        GSLResult::from(unsafe { ffi::gsl_odeiv2_evolve_apply_fixed_step(self.e, c.c, s.s, psys, t, h, y.as_mut_ptr()) })
    }

    /// This function resets the evolution function e. It should be used whenever the next use of e will not be a continuation of a previous
    /// step.
    pub fn reset(&mut self) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_evolve_reset(self.e) })
    }

    /// This function sets a pointer of the driver object d for evolve object e.
    ///
    /// If a system has discontinuous changes in the derivatives at known points, it is advisable to evolve the system between each discontinuity
    /// in sequence. For example, if a step-change in an external driving force occurs at times t_a, t_b and t_c then evolution should be carried
    /// out over the ranges (t_0,t_a), (t_a,t_b), (t_b,t_c), and (t_c,t_1) separately and not directly over the range (t_0,t_1).
    pub fn set_driver(&mut self, d: &ODEiv2Driver) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_evolve_set_driver(self.e, d.d) })
    }
}

impl Drop for ODEiv2Evolve {
    fn drop(&mut self) {
        unsafe { ffi::gsl_odeiv2_evolve_free(self.e) };
        self.e = ::std::ptr::null_mut();
    }
}

impl ffi::FFI<ffi::gsl_odeiv2_evolve> for ODEiv2Evolve {
    fn wrap(e: *mut ffi::gsl_odeiv2_evolve) -> ODEiv2Evolve {
        ODEiv2Evolve {
            e: e
        }
    }

    fn unwrap(e: &ODEiv2Evolve) -> *mut ffi::gsl_odeiv2_evolve {
        e.e
    }
}

pub struct ODEiv2Driver<'a> {
    d: *mut ffi::gsl_odeiv2_driver,
    /// `ffi::gsl_odeiv2_system` provided when constructing `d`.
    #[allow(dead_code)]
    raw_system: Box<ffi::gsl_odeiv2_system>,
    /// `PhantomData` to bind lifetime of this struct to the lifetime
    /// of the provided ODEiv2System.
    phantom: ::std::marker::PhantomData<&'a ODEiv2System<'a>>,
}

impl<'a> ODEiv2Driver<'a> {
    /// These functions return a pointer to a newly allocated instance of a driver object. The functions automatically allocate and initialise
    /// the evolve, control and stepper objects for ODE system sys using stepper type T. The initial step size is given in hstart. The rest
    /// of the arguments follow the syntax and semantics of the control functions with same name (gsl_odeiv2_control_*_new).
    pub fn alloc_y_new(sys: &'a mut ODEiv2System, t: &ODEiv2StepType, hstart: f64, epsabs: f64, epsrel: f64) -> Option<ODEiv2Driver<'a>> {
        let sys_raw = Box::new(sys.to_raw());
        let psys = &*sys_raw as *const _;
        let tmp = unsafe { ffi::gsl_odeiv2_driver_alloc_y_new(psys, t.t, hstart, epsabs, epsrel) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Driver {
                d: tmp,
                raw_system: sys_raw,
                phantom: ::std::marker::PhantomData,
            })
        }
    }

    /// These functions return a pointer to a newly allocated instance of a driver object. The functions automatically allocate and initialise
    /// the evolve, control and stepper objects for ODE system sys using stepper type T. The initial step size is given in hstart. The rest
    /// of the arguments follow the syntax and semantics of the control functions with same name (gsl_odeiv2_control_*_new).
    pub fn alloc_yp_new(sys: &'a mut ODEiv2System, t: &ODEiv2StepType, hstart: f64, epsabs: f64, epsrel: f64) -> Option<ODEiv2Driver<'a>> {
        let sys_raw = Box::new(sys.to_raw());
        let psys = &*sys_raw as *const _;
        let tmp = unsafe { ffi::gsl_odeiv2_driver_alloc_yp_new(psys, t.t, hstart, epsabs, epsrel) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Driver {
                d: tmp,
                raw_system: sys_raw,
                phantom: ::std::marker::PhantomData,
            })
        }
    }

    /// These functions return a pointer to a newly allocated instance of a driver object. The functions automatically allocate and initialise
    /// the evolve, control and stepper objects for ODE system sys using stepper type T. The initial step size is given in hstart. The rest
    /// of the arguments follow the syntax and semantics of the control functions with same name (gsl_odeiv2_control_*_new).
    pub fn alloc_standard_new(sys: &'a mut ODEiv2System, t: &ODEiv2StepType, hstart: f64, epsabs: f64, epsrel: f64, a_y: f64,
        a_dydt: f64) -> Option<ODEiv2Driver<'a>> {
        let sys_raw = Box::new(sys.to_raw());
        let psys = &*sys_raw as *const _;
        let tmp = unsafe { ffi::gsl_odeiv2_driver_alloc_standard_new(psys, t.t, hstart, epsabs, epsrel, a_y, a_dydt) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Driver {
                d: tmp,
                raw_system: sys_raw,
                phantom: ::std::marker::PhantomData,
            })
        }
    }

    /// These functions return a pointer to a newly allocated instance of a driver object. The functions automatically allocate and initialise
    /// the evolve, control and stepper objects for ODE system sys using stepper type T. The initial step size is given in hstart. The rest
    /// of the arguments follow the syntax and semantics of the control functions with same name (gsl_odeiv2_control_*_new).
    pub fn alloc_scaled_new(sys: &'a mut ODEiv2System, t: &ODEiv2StepType, hstart: f64, epsabs: f64, epsrel: f64, a_y: f64, a_dydt: f64,
        scale_abs: &[f64]) -> Option<ODEiv2Driver<'a>> {
        let sys_raw = Box::new(sys.to_raw());
        let psys = &*sys_raw as *const _;
        let tmp = unsafe { ffi::gsl_odeiv2_driver_alloc_scaled_new(psys, t.t, hstart, epsabs, epsrel, a_y, a_dydt,
            scale_abs.as_ptr()) };

        if tmp.is_null() {
            None
        } else {
            Some(ODEiv2Driver {
                d: tmp,
                raw_system: sys_raw,
                phantom: ::std::marker::PhantomData,
            })
        }
    }

    /// The function sets a minimum for allowed step size hmin for driver self. Default value is 0.
    pub fn set_hmin(&mut self, hmin: f64) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_driver_set_hmin(self.d, hmin) })
    }

    /// The function sets a maximum for allowed step size hmax for driver self. Default value is ::DBL_MAX.
    pub fn set_hmax(&mut self, hmax: f64) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_driver_set_hmax(self.d, hmax) })
    }

    /// The function sets a maximum for allowed number of steps nmax for driver self. Default value of 0 sets no limit for steps.
    pub fn set_nmax(&mut self, nmax: usize) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_driver_set_nmax(self.d, nmax) })
    }

    /// This function evolves the driver system d from t to t1. Initially vector y should contain the values of dependent variables at
    /// point t. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply is returned, and t and
    /// y contain the values from last successful step.
    ///
    /// If maximum number of steps is reached, a value of enums::Value::MaxIteration is returned. If the step size drops below minimum value, the
    /// function returns with ::NoProg. If the user-supplied functions defined in the system sys returns enums::value::BadFunc, the function
    /// returns immediately with the same return code. In this case the user must call gsl_odeiv2_driver_reset before calling this
    /// function again.
    pub fn apply(&mut self, t: &mut f64, t1: f64, y: &mut [f64]) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_driver_apply(self.d, t, t1, y.as_mut_ptr()) })
    }

    /// This function evolves the driver system d from t with n steps of size h. If the function is unable to complete the calculation, an
    /// error code from gsl_odeiv2_evolve_apply_fixed_step is returned, and t and y contain the values from last successful step.
    pub fn apply_fixed_step(&mut self, t: &mut f64, h: f64, n: usize, y: &mut [f64]) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_driver_apply_fixed_step(self.d, t, h, n, y.as_mut_ptr()) })
    }

    /// This function resets the evolution and stepper objects.
    pub fn reset(&mut self) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_driver_reset(self.d) })
    }

    /// The routine resets the evolution and stepper objects and sets new initial step size to hstart. This function can be used e.g. to
    /// change the direction of integration.
    pub fn reset_hstart(&mut self, hstart: f64) -> GSLResult<()> {
        GSLResult::from(unsafe { ffi::gsl_odeiv2_driver_reset_hstart(self.d, hstart) })
    }
}

impl<'a> Drop for ODEiv2Driver<'a> {
    fn drop(&mut self) {
        unsafe { ffi::gsl_odeiv2_driver_free(self.d) };
        self.d = ::std::ptr::null_mut();
    }
}

// We cannot wrap a driver object since we need a boxed gsl_odeiv2_system.
// impl<'a> ffi::FFI<ffi::gsl_odeiv2_driver> for ODEiv2Driver<'a> {
//     fn wrap(d: *mut ffi::gsl_odeiv2_driver) -> ODEiv2Driver<'a> {
//         ODEiv2Driver {
//             d: d,
//             phantom: ::std::marker::PhantomData,
//         }
//     }
//
//     fn unwrap(d: &ODEiv2Driver) -> *mut ffi::gsl_odeiv2_driver {
//         d.d
//     }
// }
