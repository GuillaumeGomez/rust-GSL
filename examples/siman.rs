#![allow(non_snake_case)]

extern crate rgsl;

// TODO: add TSP example

/*
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_siman.h>

/* set up parameters for this simulated annealing run */

/* how many points do we try before stepping */
#define N_TRIES 200
/* how many iterations for each T? */
#define ITERS_FIXED_T 1000
/* max step size in random walk */
#define STEP_SIZE 1.0
/* Boltzmann constant */
#define K 1.0
/* initial temperature */
#define T_INITIAL 0.008
/* damping factor for temperature */
#define MU_T 1.003
#define T_MIN 2.0e-6

gsl_siman_params_t params
= {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
K, T_INITIAL, MU_T, T_MIN};

/* now some functions to test in one dimension */
double E1(void *xp)
{
    double x = * ((double *) xp);
    return exp(-pow((x-1.0),2.0))*sin(8*x);
}

double M1(void *xp, void *yp)
{
    double x = *((double *) xp);
    double y = *((double *) yp);
    return fabs(x - y);
}

void S1(const gsl_rng * r, void *xp, double step_size)
{
    double old_x = *((double *) xp);
    double new_x;
    double u = gsl_rng_uniform(r);
    new_x = u * 2 * step_size - step_size + old_x;
    memcpy(xp, &new_x, sizeof(new_x));
}

void P1(void *xp)
{
    printf ("%12g", *((double *) xp));
}

int
main(int argc, char *argv[])
{
    const gsl_rng_type * T;
    gsl_rng * r;
    double x_initial = 15.5;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_siman_solve(r, &x_initial, E1, S1, M1, P1,
                    NULL, NULL, NULL,
                    sizeof(double), params);
    gsl_rng_free (r);
    return 0;
}
*/

fn E1(x: &f64) -> f64 {
    (-(x - 1.0).powi(2)).exp() * (x * 8.0).sin()
}

fn M1(x: &f64, y: &f64) -> f64 {
    (x - y).abs()
}

fn S1(r: &mut rgsl::Rng, x: &mut f64, step_size: f64) {
    let new_x = r.uniform() * 2.0 * step_size - step_size + *x;
    *x = new_x;
}

fn P1(x: &f64) {
    print!("{: >15.10}", x);
}

fn main() {
    rgsl::RngType::env_setup();
    let t: rgsl::RngType = rgsl::rng::default();
    let mut r = rgsl::Rng::new(&t).unwrap();

    let params = rgsl::SimAnnealingParams::new(200, 1000, 1.0, 1.0, 0.008, 1.003, 2.0e-6);
    // single run
    let siman = rgsl::SimAnnealing::new(15.5, E1, S1, M1, Some(P1), params);
    siman.solve(&mut r);
    println!(); // ignore this line, just there to split the log of both method calls
    siman.solve_many(&mut r);
}
