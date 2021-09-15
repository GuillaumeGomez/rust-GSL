//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

/*!
# 25 Simulated Annealing

Stochastic search techniques are used when the structure of a space is not well understood or
is not smooth, so that techniques like Newton’s method (which requires calculating Jacobian
derivative matrices) cannot be used. In particular, these techniques are frequently used to
solve combinatorial optimization problems, such as the traveling salesman problem.

The goal is to find a point in the space at which a real valued energy function (or cost
function) is minimized. Simulated annealing is a minimization technique which has given
good results in avoiding local minima; it is based on the idea of taking a random walk
through the space at successively lower temperatures, where the probability of taking a
step is given by a Boltzmann distribution.

The functions described in this chapter are declared in the header file gsl_siman.h.

## Simulated Annealing algorithm

The simulated annealing algorithm takes random walks through the problem space, looking
for points with low energies; in these random walks, the probability of taking a step is
determined by the Boltzmann distribution,

> p = e −(E i+1 −E i )/(kT )
> if E i+1 > E i , and p = 1 when E i+1 ≤ E i

In other words, a step will occur if the new energy is lower. If the new energy is higher,
the transition can still occur, and its likelihood is proportional to the temperature T and
inversely proportional to the energy difference E i+1 − E i .

The temperature T is initially set to a high value, and a random walk is carried out
at that temperature. Then the temperature is lowered very slightly according to a cooling
schedule, for example: T → T /μ T where μ T is slightly greater than 1.

The slight probability of taking a step that gives higher energy is what allows simulated
annealing to frequently get out of local minima.
!*/

const GSL_LOG_DBL_MIN: f64 = -7.0839641853226408e+02;

pub struct SimAnnealing<T: Clone> {
    x0_p: T,
    params: SimAnnealingParams,
    Efunc_t: gsl_siman_Efunc_t<T>,
    step_t: gsl_siman_step_t<T>,
    metric_t: gsl_siman_metric_t<T>,
    print_t: Option<gsl_siman_print_t<T>>,
}

type gsl_siman_Efunc_t<T> = fn(&T) -> f64;
type gsl_siman_step_t<T> = fn(&mut ::Rng, &mut T, f64);
type gsl_siman_metric_t<T> = fn(&T, &T) -> f64;
type gsl_siman_print_t<T> = fn(&T);

fn boltzmann(E: f64, new_E: f64, T: f64, params: &SimAnnealingParams) -> f64 {
    let x = -(new_E - E) / (params.k * T);
    // avoid underflow errors for large uphill steps
    if x < GSL_LOG_DBL_MIN {
        0.0
    } else {
        x.exp()
    }
}

impl<T> SimAnnealing<T>
where
    T: Clone,
{
    pub fn new(
        x0_p: T,
        ef: gsl_siman_Efunc_t<T>,
        take_step: gsl_siman_step_t<T>,
        distance: gsl_siman_metric_t<T>,
        print_pos: Option<gsl_siman_print_t<T>>,
        params: SimAnnealingParams,
    ) -> SimAnnealing<T> {
        SimAnnealing {
            x0_p,
            params,
            Efunc_t: ef,
            step_t: take_step,
            metric_t: distance,
            print_t: print_pos,
        }
    }

    /*
    /* implementation of a basic simulated annealing algorithm */

    void
    gsl_siman_solve (const gsl_rng * r, void *x0_p, gsl_siman_Efunc_t Ef,
                    gsl_siman_step_t take_step,
                    gsl_siman_metric_t distance,
                    gsl_siman_print_t print_position,
                    gsl_siman_copy_t copyfunc,
                    gsl_siman_copy_construct_t copy_constructor,
                    gsl_siman_destroy_t destructor,
                    size_t element_size,
                    gsl_siman_params_t params)
    {
        void *x, *new_x, *best_x;
        double E, new_E, best_E;
        int i;
        double T, T_factor;
        int n_evals = 1, n_iter = 0, n_accepts, n_rejects, n_eless;

        /* this function requires that either the dynamic functions (copy,
            copy_constructor and destrcutor) are passed, or that an element
            size is given */
        assert((copyfunc != NULL && copy_constructor != NULL && destructor != NULL)
                || (element_size != 0));

        distance = 0 ; /* This parameter is not currently used */
        E = Ef(x0_p);

        if (copyfunc) {
            x = copy_constructor(x0_p);
            new_x = copy_constructor(x0_p);
            best_x = copy_constructor(x0_p);
        } else {
            x = (void *) malloc (element_size);
            memcpy (x, x0_p, element_size);
            new_x = (void *) malloc (element_size);
            best_x =  (void *) malloc (element_size);
            memcpy (best_x, x0_p, element_size);
        }

        best_E = E;

        T = params.t_initial;
        T_factor = 1.0 / params.mu_t;

        if (print_position) {
            printf ("#-iter  #-evals   temperature     position   energy\n");
        }

        while (1) {

            n_accepts = 0;
            n_rejects = 0;
            n_eless = 0;

            for (i = 0; i < params.iters_fixed_T; ++i) {

                copy_state(x, new_x, element_size, copyfunc);

                take_step (r, new_x, params.step_size);
                new_E = Ef (new_x);

                if(new_E <= best_E){
                    if (copyfunc) {
                        copyfunc(new_x,best_x);
                    } else {
                        memcpy (best_x, new_x, element_size);
                    }
                    best_E=new_E;
                }

                ++n_evals;                /* keep track of Ef() evaluations */
                /* now take the crucial step: see if the new point is accepted
                    or not, as determined by the boltzmann probability */
                if (new_E < E) {
                    if (new_E < best_E) {
                        copy_state(new_x, best_x, element_size, copyfunc);
                        best_E = new_E;
                    }

                    /* yay! take a step */
                    copy_state(new_x, x, element_size, copyfunc);
                    E = new_E;
                    ++n_eless;

                } else if (gsl_rng_uniform(r) < boltzmann(E, new_E, T, &params)) {
                    /* yay! take a step */
                    copy_state(new_x, x, element_size, copyfunc);
                        E = new_E;
                        ++n_accepts;

                } else {
                    ++n_rejects;
                }
            }

            if (print_position) {
                /* see if we need to print stuff as we go */
                /*       printf("%5d %12g %5d %3d %3d %3d", n_iter, T, n_evals, */
                /*           100*n_eless/n_steps, 100*n_accepts/n_steps, */
                /*           100*n_rejects/n_steps); */
                printf ("%5d   %7d  %12g", n_iter, n_evals, T);
                print_position (x);
                printf ("  %12g  %12g\n", E, best_E);
            }

            /* apply the cooling schedule to the temperature */
            /* FIXME: I should also introduce a cooling schedule for the iters */
            T *= T_factor;
            ++n_iter;
            if (T < params.t_min) {
                break;
            }
        }

        /* at the end, copy the result onto the initial point, so we pass it
            back to the caller */
        copy_state(best_x, x0_p, element_size, copyfunc);

        if (copyfunc) {
            destructor(x);
            destructor(new_x);
            destructor(best_x);
        } else {
            free (x);
            free (new_x);
            free (best_x);
        }
    }
    */
    /// This function performs a simulated annealing search through a given space. The space
    /// is specified by providing the functions Ef and distance. The simulated annealing steps
    /// are generated using the random number generator `rng` and the function `take_step`.
    ///
    /// The starting configuration of the system should be given by x0\_p.
    ///
    /// The params structure (described below) controls the run by providing the temperature
    /// schedule and other tunable parameters to the algorithm.
    ///
    /// On exit the best result achieved during the search is returned. If the annealing
    /// process has been successful this should be a good approximation to the optimal point
    /// in the space.
    ///
    /// If the argument `print_pos` is not None, a debugging log will be printed to
    /// stdout with the following columns: ```#-iter #-evals temperature position energy best_energy```
    /// and the output of the function print position itself.
    pub fn solve(&self, rng: &mut ::Rng) -> T {
        let mut x = self.x0_p.clone();
        let mut new_x = self.x0_p.clone();
        let mut best_x = self.x0_p.clone();

        let mut n_evals = 0_usize;

        let mut E = (self.Efunc_t)(&self.x0_p);
        let mut best_E = E;

        let mut Temp = self.params.t_initial;
        let Temp_factor = 1.0 / self.params.mu_t;

        if self.print_t.is_some() {
            println!(
                "{i:^6} | {e:^7} | {t:^12} | {p:^15} | {E:^13}",
                i = "#-iter",
                e = "#-evals",
                t = "temperature",
                p = "position",
                E = "energy"
            );
        }

        let mut iter = 0;
        loop {
            //let mut n_accepts = 0;
            //let mut n_rejects = 0;
            //let mut n_eless = 0;

            for _ in 0..self.params.iters_fixed_T {
                x = new_x.clone();

                (self.step_t)(rng, &mut new_x, self.params.step_size);
                let new_E = (self.Efunc_t)(&new_x);
                n_evals += 1; // keep track of Ef() evaluations

                if new_E <= best_E {
                    best_x = new_x.clone();
                    best_E = new_E;
                }

                // now take the crucial step: see if the new point is accepted
                // or not, as determined by the boltzmann probability
                if new_E < E {
                    if new_E < best_E {
                        best_x = new_x.clone();
                        best_E = new_E;
                    }
                    // yay! take a step
                    x = new_x.clone();
                    E = new_E;
                //n_eless += 1;
                } else if rng.uniform() < boltzmann(E, new_E, Temp, &self.params) {
                    // yay! take a step
                    x = new_x.clone();
                    E = new_E;
                    //n_accepts += 1;
                }
            }

            if let Some(ref printf) = self.print_t {
                print!("{:>06} | {:>07} | {:>12.10} | ", iter, n_evals, Temp);
                printf(&x);
                println!(" | {:+>13.12}", E);
            }

            Temp *= Temp_factor;
            iter += 1;
            if Temp < self.params.t_min {
                break;
            }
        }

        best_x
    }

    /*
    /* implementation of a simulated annealing algorithm with many tries */

    void
    gsl_siman_solve_many (const gsl_rng * r, void *x0_p, gsl_siman_Efunc_t Ef,
                        gsl_siman_step_t take_step,
                        gsl_siman_metric_t distance,
                        gsl_siman_print_t print_position,
                        size_t element_size,
                        gsl_siman_params_t params)
    {
        /* the new set of trial points, and their energies and probabilities */
        void *x, *new_x;
        double *energies, *probs, *sum_probs;
        double Ex;                    /* energy of the chosen point */
        double T, T_factor;           /* the temperature and a step multiplier */
        int i;
        double u;                     /* throw the die to choose a new "x" */
        int n_iter;

        if (print_position) {
            printf ("#-iter    temperature       position");
            printf ("         delta_pos        energy\n");
        }

        x = (void *) malloc (params.n_tries * element_size);
        new_x = (void *) malloc (params.n_tries * element_size);
        energies = (double *) malloc (params.n_tries * sizeof (double));
        probs = (double *) malloc (params.n_tries * sizeof (double));
        sum_probs = (double *) malloc (params.n_tries * sizeof (double));

        T = params.t_initial;
        T_factor = 1.0 / params.mu_t;

        memcpy (x, x0_p, element_size);

        n_iter = 0;
        while (1)
            {
            Ex = Ef (x);
            for (i = 0; i < params.n_tries - 1; ++i)
                {                       /* only go to N_TRIES-2 */
                /* center the new_x[] around x, then pass it to take_step() */
                sum_probs[i] = 0;
                memcpy ((char *)new_x + i * element_size, x, element_size);
                take_step (r, (char *)new_x + i * element_size, params.step_size);
                energies[i] = Ef ((char *)new_x + i * element_size);
                probs[i] = boltzmann(Ex, energies[i], T, &params);
                }
            /* now add in the old value of "x", so it is a contendor */
            memcpy ((char *)new_x + (params.n_tries - 1) * element_size, x, element_size);
            energies[params.n_tries - 1] = Ex;
            probs[params.n_tries - 1] = boltzmann(Ex, energies[i], T, &params);

            /* now throw biased die to see which new_x[i] we choose */
            sum_probs[0] = probs[0];
            for (i = 1; i < params.n_tries; ++i)
                {
                sum_probs[i] = sum_probs[i - 1] + probs[i];
                }
            u = gsl_rng_uniform (r) * sum_probs[params.n_tries - 1];
            for (i = 0; i < params.n_tries; ++i)
                {
                if (u < sum_probs[i])
                    {
                    memcpy (x, (char *) new_x + i * element_size, element_size);
                    break;
                    }
                }
            if (print_position)
                {
                printf ("%5d\t%12g\t", n_iter, T);
                print_position (x);
                printf ("\t%12g\t%12g\n", distance (x, x0_p), Ex);
                }
            T *= T_factor;
            ++n_iter;
            if (T < params.t_min)
            {
            break;
                }
            }

        /* now return the value via x0_p */
        memcpy (x0_p, x, element_size);

        /*  printf("the result is: %g (E=%g)\n", x, Ex); */

        free (x);
        free (new_x);
        free (energies);
        free (probs);
        free (sum_probs);
    }
    */
    /// Like the function solve, but performs multiple runs and returns the best result.
    pub fn solve_many(&self, rng: &mut ::Rng) -> T {
        let mut x = self.x0_p.clone();
        let mut new_x = Vec::with_capacity(self.params.n_tries);

        let mut energies = Vec::with_capacity(self.params.n_tries);
        let mut probs = Vec::with_capacity(self.params.n_tries);
        let mut sum_probs = Vec::with_capacity(self.params.n_tries);

        let mut Temp = self.params.t_initial;
        let Temp_factor = 1.0 / self.params.mu_t;

        if self.print_t.is_some() {
            println!(
                "{i:^6} | {t:^12} | {p:^15} | {d:^15} | {E:^13}",
                i = "#-iter",
                t = "temperature",
                p = "position",
                d = "delta_pos",
                E = "energy"
            );
        }

        let mut iter = 0;
        loop {
            let Ex = (self.Efunc_t)(&x);
            for i in 0..self.params.n_tries - 1 {
                // only go to N_TRIES-2
                sum_probs.push(0.0);
                new_x.push(x.clone());

                (self.step_t)(rng, &mut new_x[i], self.params.step_size);
                energies.push((self.Efunc_t)(&new_x[i]));
                probs.push(boltzmann(Ex, energies[i], Temp, &self.params));
            }

            // now add in the old value of "x", so it is a contendor
            new_x.push(x.clone());
            energies.push(Ex);
            probs.push(boltzmann(
                Ex,
                energies[self.params.n_tries - 2],
                Temp,
                &self.params,
            ));
            sum_probs.push(0.0);

            // now throw biased die to see which new_x[i] we choose
            sum_probs[0] = probs[0];
            for i in 1..self.params.n_tries {
                sum_probs[i] = sum_probs[i - 1] + probs[i];
            }
            let u = rng.uniform() * *sum_probs.last().unwrap();
            for i in 0..self.params.n_tries {
                if u < sum_probs[i] {
                    x = new_x[i].clone();
                    break;
                }
            }

            if let Some(ref printf) = self.print_t {
                print!("{:>06} | {:>12.10} | ", iter, Temp);
                printf(&x);
                println!(
                    " | {: >15.12} | {: >13.12}",
                    (self.metric_t)(&x, &self.x0_p),
                    Ex
                );
            }

            Temp *= Temp_factor;
            iter += 1;
            if Temp < self.params.t_min {
                break;
            }
        }
        x
    }
}

pub struct SimAnnealingParams {
    n_tries: usize,
    iters_fixed_T: usize,
    step_size: f64,
    k: f64,
    t_initial: f64,
    mu_t: f64,
    t_min: f64,
}

impl SimAnnealingParams {
    /// These are the parameters that control a run of the simulated annealing algorithm.
    /// This structure contains all the information needed to control the search,
    /// beyond the energy function, the step function and the initial guess.
    ///
    /// - n_tries: The number of points to try for each step.
    /// - iters: The number of iterations at each temperature.
    /// - step_size: The maximum step size in the random walk.
    /// - k, t_initial, mu_t, t_min: The parameters of the Boltzmann distribution and
    ///   cooling schedule.
    pub fn new(
        n_tries: usize,
        iters: usize,
        step_size: f64,
        k: f64,
        t_initial: f64,
        mu_t: f64,
        t_min: f64,
    ) -> SimAnnealingParams {
        SimAnnealingParams {
            n_tries,
            iters_fixed_T: iters,
            step_size,
            k,
            t_initial,
            mu_t,
            t_min,
        }
    }
}
