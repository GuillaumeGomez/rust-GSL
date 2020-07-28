//! Random distributions

use libc::{c_double, c_int, c_uint, c_void, size_t};

use super::{gsl_ran_discrete_t, gsl_rng};

extern "C" {
    // Random Number Distributions
    // The Gaussian Distribution
    pub fn gsl_ran_gaussian(r: *mut gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_pdf(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_ziggurat(r: *mut gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_ratio_method(r: *mut gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_ugaussian(r: *mut gsl_rng) -> c_double;
    pub fn gsl_ran_ugaussian_pdf(x: c_double) -> c_double;
    pub fn gsl_ran_ugaussian_ratio_method(r: *mut gsl_rng) -> c_double;
    pub fn gsl_cdf_gaussian_P(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_gaussian_Q(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_gaussian_Pinv(P: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_gaussian_Qinv(Q: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_P(x: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_Q(x: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_Pinv(P: c_double) -> c_double;
    pub fn gsl_cdf_ugaussian_Qinv(Q: c_double) -> c_double;
    // The Gaussian Tail Distribution
    pub fn gsl_ran_gaussian_tail(r: *mut gsl_rng, a: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_gaussian_tail_pdf(x: c_double, a: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_ugaussian_tail(r: *mut gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_ugaussian_tail_pdf(x: c_double, a: c_double) -> c_double;
    // The Bivariate Gaussian Distribution
    pub fn gsl_ran_bivariate_gaussian(
        r: *mut gsl_rng,
        sigma_x: c_double,
        sigma_y: c_double,
        rho: c_double,
        x: *mut c_double,
        y: *mut c_double,
    );
    pub fn gsl_ran_bivariate_gaussian_pdf(
        x: c_double,
        y: c_double,
        sigma_x: c_double,
        sigma_y: c_double,
        rho: c_double,
    ) -> c_double;
    // The Exponential Distribution
    pub fn gsl_ran_exponential(r: *mut gsl_rng, mu: c_double) -> c_double;
    pub fn gsl_ran_exponential_pdf(x: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_P(x: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_Q(x: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_Pinv(P: c_double, mu: c_double) -> c_double;
    pub fn gsl_cdf_exponential_Qinv(Q: c_double, mu: c_double) -> c_double;
    // The Laplace Distribution
    pub fn gsl_ran_laplace(r: *mut gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_laplace_pdf(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_P(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_Q(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_Pinv(P: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_laplace_Qinv(Q: c_double, a: c_double) -> c_double;
    // The Exponential Power Distribution
    pub fn gsl_ran_exppow(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_exppow_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_exppow_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_exppow_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    // The Cauchy Distribution
    pub fn gsl_ran_cauchy(r: *mut gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_cauchy_pdf(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_P(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_Q(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_Pinv(P: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_cauchy_Qinv(Q: c_double, a: c_double) -> c_double;
    // The Rayleigh Distribution
    pub fn gsl_ran_rayleigh(r: *mut gsl_rng, sigma: c_double) -> c_double;
    pub fn gsl_ran_rayleigh_pdf(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_P(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_Q(x: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_Pinv(P: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_rayleigh_Qinv(Q: c_double, sigma: c_double) -> c_double;
    // The Rayleigh Tail Distribution
    pub fn gsl_ran_rayleigh_tail(r: *mut gsl_rng, a: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_rayleigh_tail_pdf(x: c_double, a: c_double, sigma: c_double) -> c_double;
    // The Landau Distribution
    pub fn gsl_ran_landau(r: *mut gsl_rng) -> c_double;
    pub fn gsl_ran_landau_pdf(x: c_double) -> c_double;
    // The Levy alpha-Stable Distributions
    pub fn gsl_ran_levy(r: *mut gsl_rng, c: c_double, alpha: c_double) -> c_double;
    // The Levy skew alpha-Stable Distribution
    pub fn gsl_ran_levy_skew(
        r: *mut gsl_rng,
        c: c_double,
        alpha: c_double,
        beta: c_double,
    ) -> c_double;
    // The Gamma Distribution
    pub fn gsl_ran_gamma(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gamma_knuth(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gamma_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gamma_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Flat (Uniform) Distribution
    pub fn gsl_ran_flat(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_flat_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_flat_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Lognormal Distribution
    pub fn gsl_ran_lognormal(r: *mut gsl_rng, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_ran_lognormal_pdf(x: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_P(x: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_Q(x: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_Pinv(P: c_double, zeta: c_double, sigma: c_double) -> c_double;
    pub fn gsl_cdf_lognormal_Qinv(Q: c_double, zeta: c_double, sigma: c_double) -> c_double;
    // The Chi-squared Distribution
    pub fn gsl_ran_chisq(r: *mut gsl_rng, nu: c_double) -> c_double;
    pub fn gsl_ran_chisq_pdf(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_P(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_Q(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_Pinv(P: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_chisq_Qinv(Q: c_double, nu: c_double) -> c_double;
    // The F-distribution
    pub fn gsl_ran_fdist(r: *mut gsl_rng, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_ran_fdist_pdf(x: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_P(x: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_Q(x: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_Pinv(P: c_double, nu1: c_double, nu2: c_double) -> c_double;
    pub fn gsl_cdf_fdist_Qinv(Q: c_double, nu1: c_double, nu2: c_double) -> c_double;
    // The t-distribution
    pub fn gsl_ran_tdist(r: *mut gsl_rng, nu: c_double) -> c_double;
    pub fn gsl_ran_tdist_pdf(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_P(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_Q(x: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_Pinv(P: c_double, nu: c_double) -> c_double;
    pub fn gsl_cdf_tdist_Qinv(Q: c_double, nu: c_double) -> c_double;
    // The Beta Distribution
    pub fn gsl_ran_beta(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_beta_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_beta_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Logistic Distribution
    pub fn gsl_ran_logistic(r: *mut gsl_rng, a: c_double) -> c_double;
    pub fn gsl_ran_logistic_pdf(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_P(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_Q(x: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_Pinv(P: c_double, a: c_double) -> c_double;
    pub fn gsl_cdf_logistic_Qinv(Q: c_double, a: c_double) -> c_double;
    // The Pareto Distribution
    pub fn gsl_ran_pareto(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_pareto_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_pareto_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // Spherical Vector Distributions
    pub fn gsl_ran_dir_2d(r: *mut gsl_rng, x: *mut c_double, y: *mut c_double);
    pub fn gsl_ran_dir_2d_trig_method(r: *mut gsl_rng, x: *mut c_double, y: *mut c_double);
    pub fn gsl_ran_dir_3d(r: *mut gsl_rng, x: *mut c_double, y: *mut c_double, z: *mut c_double);
    pub fn gsl_ran_dir_nd(r: *mut gsl_rng, n: size_t, x: *mut c_double);
    // The Weibull Distribution
    pub fn gsl_ran_weibull(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_weibull_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_weibull_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Type-1 Gumbel Distribution
    pub fn gsl_ran_gumbel1(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gumbel1_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel1_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Type-2 Gumbel Distribution
    pub fn gsl_ran_gumbel2(r: *mut gsl_rng, a: c_double, b: c_double) -> c_double;
    pub fn gsl_ran_gumbel2_pdf(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_P(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_Q(x: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_Pinv(P: c_double, a: c_double, b: c_double) -> c_double;
    pub fn gsl_cdf_gumbel2_Qinv(Q: c_double, a: c_double, b: c_double) -> c_double;
    // The Dirichlet Distribution
    pub fn gsl_ran_dirichlet(
        r: *mut gsl_rng,
        K: size_t,
        alpha: *const c_double,
        theta: *mut c_double,
    );
    pub fn gsl_ran_dirichlet_pdf(
        K: size_t,
        alpha: *const c_double,
        theta: *const c_double,
    ) -> c_double;
    pub fn gsl_ran_dirichlet_lnpdf(
        K: size_t,
        alpha: *const c_double,
        theta: *const c_double,
    ) -> c_double;
    // General Discrete Distributions
    pub fn gsl_ran_discrete_preproc(K: size_t, P: *const c_double) -> *mut gsl_ran_discrete_t;
    pub fn gsl_ran_discrete(r: *mut gsl_rng, g: *const gsl_ran_discrete_t) -> size_t;
    pub fn gsl_ran_discrete_pdf(k: size_t, g: *const gsl_ran_discrete_t) -> c_double;
    pub fn gsl_ran_discrete_free(g: *mut gsl_ran_discrete_t);
    // The Poisson Distribution
    pub fn gsl_ran_poisson(r: *mut gsl_rng, mu: c_double) -> c_uint;
    pub fn gsl_ran_poisson_pdf(k: c_uint, mu: c_double) -> c_double;
    pub fn gsl_cdf_poisson_P(k: c_uint, mu: c_double) -> c_double;
    pub fn gsl_cdf_poisson_Q(k: c_uint, mu: c_double) -> c_double;
    // The Bernoulli Distribution
    pub fn gsl_ran_bernoulli(r: *mut gsl_rng, p: c_double) -> c_uint;
    pub fn gsl_ran_bernoulli_pdf(k: c_uint, p: c_double) -> c_double;
    // The Binomial Distribution
    pub fn gsl_ran_binomial(r: *mut gsl_rng, p: c_double, n: c_uint) -> c_uint;
    pub fn gsl_ran_binomial_pdf(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_binomial_P(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_binomial_Q(k: c_uint, p: c_double, n: c_uint) -> c_double;
    // The Multinomial Distribution
    pub fn gsl_ran_multinomial(
        r: *mut gsl_rng,
        K: size_t,
        N: c_uint,
        p: *const c_double,
        n: *mut c_uint,
    );
    pub fn gsl_ran_multinomial_pdf(K: size_t, p: *const c_double, n: *const c_uint) -> c_double;
    pub fn gsl_ran_multinomial_lnpdf(K: size_t, p: *const c_double, n: *const c_uint) -> c_double;
    // The Negative Binomial Distribution
    pub fn gsl_ran_negative_binomial(r: *mut gsl_rng, p: c_double, n: c_double) -> c_uint;
    pub fn gsl_ran_negative_binomial_pdf(k: c_uint, p: c_double, n: c_double) -> c_double;
    pub fn gsl_cdf_negative_binomial_P(k: c_uint, p: c_double, n: c_double) -> c_double;
    pub fn gsl_cdf_negative_binomial_Q(k: c_uint, p: c_double, n: c_double) -> c_double;
    // The Pascal Distribution
    pub fn gsl_ran_pascal(r: *mut gsl_rng, p: c_double, n: c_uint) -> c_uint;
    pub fn gsl_ran_pascal_pdf(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_pascal_P(k: c_uint, p: c_double, n: c_uint) -> c_double;
    pub fn gsl_cdf_pascal_Q(k: c_uint, p: c_double, n: c_uint) -> c_double;
    // The Geometric Distribution
    pub fn gsl_ran_geometric(r: *mut gsl_rng, p: c_double) -> c_uint;
    pub fn gsl_ran_geometric_pdf(k: c_uint, p: c_double) -> c_double;
    pub fn gsl_cdf_geometric_P(k: c_uint, p: c_double) -> c_double;
    pub fn gsl_cdf_geometric_Q(k: c_uint, p: c_double) -> c_double;
    // The Hypergeometric Distribution
    pub fn gsl_ran_hypergeometric(r: *mut gsl_rng, n1: c_uint, n2: c_uint, t: c_uint) -> c_uint;
    pub fn gsl_ran_hypergeometric_pdf(k: c_uint, n1: c_uint, n2: c_uint, t: c_uint) -> c_double;
    pub fn gsl_cdf_hypergeometric_P(k: c_uint, n1: c_uint, n2: c_uint, t: c_uint) -> c_double;
    pub fn gsl_cdf_hypergeometric_Q(k: c_uint, n1: c_uint, n2: c_uint, t: c_uint) -> c_double;
    // The Logarithmic Distribution
    pub fn gsl_ran_logarithmic(r: *mut gsl_rng, p: c_double) -> c_uint;
    pub fn gsl_ran_logarithmic_pdf(k: c_uint, p: c_double) -> c_double;
    // Shuffling and Sampling
    pub fn gsl_ran_shuffle(r: *mut gsl_rng, base: *mut c_void, n: size_t, size: size_t);
    pub fn gsl_ran_choose(
        r: *mut gsl_rng,
        dest: *mut c_void,
        k: size_t,
        src: *mut c_void,
        n: size_t,
        size: size_t,
    ) -> c_int;
    pub fn gsl_ran_sample(
        r: *mut gsl_rng,
        dest: *mut c_void,
        k: size_t,
        src: *mut c_void,
        n: size_t,
        size: size_t,
    ) -> c_int;
}
