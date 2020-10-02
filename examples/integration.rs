//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

extern crate rgsl;

struct FParams {
    // Amplitude
    a: f64,
    // Phase
    phi: f64,
}

fn f(x: f64, p: &mut FParams) -> f64 {
    (p.a * x + p.phi).sin()
}

#[allow(unused_variables)]
fn cqf1(x: f64, p: &mut f64) -> f64 {
    x.exp()
}

fn qags_fn(x: f64, alpha: &mut f64) -> f64 {
    (*alpha * x).ln() / x.sqrt()
}

/* f458(x) = 1/(1 + log(x)^2)^2 */
/* integ(log(x) f458(x),x,0,1) = (Ci(1) sin(1) + (pi/2 - Si(1)) cos(1))/pi
= -0.1892752 */
#[allow(unused_variables)]
fn f458<T>(x: f64, params: &mut T) -> f64 {
    if x == 0f64 {
        0f64
    } else {
        let u = x.ln();
        let v = 1f64 + u * u;

        1f64 / (v * v)
    }
}

fn main() {
    let mut params = FParams { a: 1f64, phi: 0f64 };
    let mut result = 0f64;
    let mut error = 0f64;
    let mut n_eval = 0;

    let xlow = 0f64;
    let xhigh = 10f64;
    let eps_abs = 1e-4f64;
    let eps_rel = 1e-4f64;

    println!("=== integration::qng ===");
    match rgsl::integration::qng(
        f,
        xlow,
        xhigh,
        eps_abs,
        eps_rel,
        &mut result,
        &mut error,
        &mut n_eval,
    ) {
        rgsl::Value::Success => {
            println!(
                "Result {} +/- {} from {} evaluations",
                result, error, n_eval
            );
        }
        e => {
            println!("There was a problem with integration: {:?}", e);
        }
    };

    println!("\n=== IntegrationWorkspace.qag ===");
    let mut iw = rgsl::IntegrationWorkspace::new(5).unwrap();

    match iw.qag(
        f,
        &mut params,
        xlow,
        xhigh,
        eps_abs,
        eps_rel,
        1,
        rgsl::GaussKonrodRule::Gauss15,
        &mut result,
        &mut error,
    ) {
        rgsl::Value::Success => {
            println!("Result {} +/- {}", result, error);
        }
        e => {
            println!("There was a problem with integration: {:?}", e);
        }
    };

    println!("\n=== IntegrationWorkspace.qagi ===");
    let limit = iw.limit();
    match iw.qagi(
        f,
        &mut params,
        1.0e-7f64,
        0f64,
        limit,
        &mut result,
        &mut error,
    ) {
        rgsl::Value::Success => {
            println!("Result {} +/- {}", result, error);
        }
        e => {
            println!("There was a problem with integration: {:?}", e);
        }
    }

    {
        println!("\n=== IntegrationQawsTable.qaws ===");
        let mut t = rgsl::IntegrationQawsTable::new(0f64, 0f64, 1, 0).unwrap();
        let mut w = rgsl::IntegrationWorkspace::new(1000).unwrap();

        match t.qaws(
            f458,
            &mut 1f64,
            0f64,
            1f64,
            0f64,
            1.0e-7f64,
            w.limit(),
            &mut w,
            &mut result,
            &mut error,
        ) {
            rgsl::Value::Success => {
                println!("Result {} +/- {}", result, error);
            }
            e => {
                println!("There was a problem with integration: {:?}", e);
            }
        }
    }

    println!("\n=== CquadWorkspace.cquad ===");
    let mut t = rgsl::CquadWorkspace::new(200).unwrap();

    match t.cquad(
        cqf1,
        &mut 1f64,
        0f64,
        1f64,
        0f64,
        1.0e-12f64,
        &mut result,
        &mut error,
        &mut n_eval,
    ) {
        rgsl::Value::Success => {
            println!("Result {} +/- {} -> {}", result, error, n_eval);
        }
        e => {
            println!("There was a problem with integration: {:?}", e);
        }
    }

    {
        println!("\n=== IntegrationWorkspace.qags ===");
        let mut w = rgsl::IntegrationWorkspace::new(1000).unwrap();

        let expected = -4f64;
        let mut alpha = 1f64;

        w.qags(
            qags_fn,
            &mut alpha,
            0f64,
            1f64,
            0f64,
            1e-7f64,
            1000,
            &mut result,
            &mut error,
        );

        println!("result          = {:.18}", result);
        println!("exact result    = {:.18}", expected);
        println!("estimated error = {:.18}", error);
        println!("actual error    = {:.18}", result - expected);
        println!("intervals =  {}", w.size());
    }
}
