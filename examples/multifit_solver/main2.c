The following example program fits a weighted exponential model with background to experimental data, Y = A \exp(-\lambda t) + b. The first part of the program sets up the functions expb_f and expb_df to calculate the model and its Jacobian. The appropriate fitting function is given by,

f_i = ((A \exp(-\lambda t_i) + b) - y_i)/\sigma_i
where we have chosen t_i = i. The Jacobian matrix J is the derivative of these functions with respect to the three parameters (A, \lambda, b). It is given by,

J_{ij} = d f_i / d x_j
where x_0 = A, x_1 = \lambda and x_2 = b.

/* expfit.c -- model functions for exponential + background */

struct data {
  size_t n;
  double * y;
  double * sigma;
};

int
expb_f (const gsl_vector * x, void *data, 
        gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *) data)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);
  double b = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * i) + b */
      double t = i;
      double Yi = A * exp (-lambda * t) + b;
      gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
    }

  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *data, 
         gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *) data)->sigma;

  double A = gsl_vector_get (x, 0);
  double lambda = gsl_vector_get (x, 1);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      double t = i;
      double s = sigma[i];
      double e = exp(-lambda * t);
      gsl_matrix_set (J, i, 0, e/s); 
      gsl_matrix_set (J, i, 1, -t * A * e/s);
      gsl_matrix_set (J, i, 2, 1/s);
    }
  return GSL_SUCCESS;
}

int
expb_fdf (const gsl_vector * x, void *data,
          gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, data, f);
  expb_df (x, data, J);

  return GSL_SUCCESS;
}
The main part of the program sets up a Levenberg-Marquardt solver and some simulated random data. The data uses the known parameters (1.0,5.0,0.1) combined with Gaussian noise (standard deviation = 0.1) over a range of 40 timesteps. The initial guess for the parameters is chosen as (0.0, 1.0, 0.0).

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "expfit.c"

#define N 40

void print_state (size_t iter, gsl_multifit_fdfsolver * s);

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const size_t n = N;
  const size_t p = 3;

  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double y[N], sigma[N];
  struct data d = { n, y, sigma};
  gsl_multifit_function_fdf f;
  double x_init[3] = { 1.0, 0.0, 0.0 };
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  /* This is the data to be fitted */

  for (i = 0; i < n; i++)
    {
      double t = i;
      y[i] = 1.0 + 5 * exp (-0.1 * t) 
                 + gsl_ran_gaussian (r, 0.1);
      sigma[i] = 0.1;
      printf ("data: %u %g %g\n", i, y[i], sigma[i]);
    };

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      printf ("status = %s\n", gsl_strerror (status));

      print_state (iter, s);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x,
                                        1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  { 
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

    printf ("A      = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
    printf ("lambda = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
    printf ("b      = %.5f +/- %.5f\n", FIT(2), c*ERR(2));
  }

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_rng_free (r);
  return 0;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2), 
          gsl_blas_dnrm2 (s->f));
}
The iteration terminates when the change in x is smaller than 0.0001, as both an absolute and relative change. Here are the results of running the program:

iter: 0 x=1.00000000 0.00000000 0.00000000 |f(x)|=117.349
status=success
iter: 1 x=1.64659312 0.01814772 0.64659312 |f(x)|=76.4578
status=success
iter: 2 x=2.85876037 0.08092095 1.44796363 |f(x)|=37.6838
status=success
iter: 3 x=4.94899512 0.11942928 1.09457665 |f(x)|=9.58079
status=success
iter: 4 x=5.02175572 0.10287787 1.03388354 |f(x)|=5.63049
status=success
iter: 5 x=5.04520433 0.10405523 1.01941607 |f(x)|=5.44398
status=success
iter: 6 x=5.04535782 0.10404906 1.01924871 |f(x)|=5.44397
chisq/dof = 0.800996
A      = 5.04536 +/- 0.06028
lambda = 0.10405 +/- 0.00316
b      = 1.01925 +/- 0.03782
status = success
The approximate values of the parameters are found correctly, and the chi-squared value indicates a good fit (the chi-squared per degree of freedom is approximately 1). In this case the errors on the parameters can be estimated from the square roots of the diagonal elements of the covariance matrix.

If the chi-squared value shows a poor fit (i.e. chi^2/dof >> 1) then the error estimates obtained from the covariance matrix will be too small. In the example program the error estimates are multiplied by \sqrt{\chi^2/dof} in this case, a common way of increasing the errors for a poor fit. Note that a poor fit will result from the use an inappropriate model, and the scaled error estimates may then be outside the range of validity for Gaussian errors.