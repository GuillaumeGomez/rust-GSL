#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>

typedef double (*opt_function)(void*, double, void*);

extern double rust_callback_f(opt_function, const gsl_vector*, size_t, double, const gsl_vector*, size_t);
extern double rust_callback_dfs(opt_function*, const gsl_vector*, size_t, size_t, double, const gsl_vector*, size_t);

struct data {
  opt_function func_f;
  opt_function* func_dfs;
  size_t params_len;
  double* ts;
  double* ys;
  size_t vars_len;
  double* args;
  size_t args_len;
};

int
call_f (const gsl_vector * x, void *data,
        gsl_vector * f)
{
  size_t params_len = ((struct data *)data)->params_len;
  double *t = ((struct data *)data)->ts;
  double *y = ((struct data *)data)->ys;
  size_t n = ((struct data *)data)->vars_len;

  opt_function func_f = ((struct data *)data)->func_f;

  double *args = ((struct data *)data)->args;
  size_t args_len = ((struct data *)data)->args_len;
  gsl_vector_view args_vec = gsl_vector_view_array(args, args_len);

  for (size_t i = 0; i < n; i++)
    {
      /* Model Yi = A * exp(-lambda * t_i) + b */
      double Yi = rust_callback_f(func_f, x, params_len, t[i], &args_vec.vector, args_len); 
      gsl_vector_set (f, i, Yi - y[i]);
    }

  return GSL_SUCCESS;
}

int
call_dfs (const gsl_vector * x, void *data,
        gsl_matrix * jacobian)
{
  size_t params_len = ((struct data *)data)->params_len;
  double *t = ((struct data *)data)->ts;
  size_t n = ((struct data *)data)->vars_len;

  opt_function* func_dfs = ((struct data *)data)->func_dfs;

  double *args = ((struct data *)data)->args;
  size_t args_len = ((struct data *)data)->args_len;
  gsl_vector_view args_vec = gsl_vector_view_array(args, args_len);

  for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < params_len; j++) {
          gsl_matrix_set (jacobian, i, j, rust_callback_dfs(func_dfs, x, params_len, j, t[i], &args_vec.vector, args_len));
      }
  }

  return GSL_SUCCESS;
}


void
solver_progress(const size_t iter, const size_t vector_len,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *x = w->x;
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  printf("iter %2zu: ", iter);
  for (size_t i = 0; i < vector_len; i++) {
      printf("par[%lu] = %.4f, ", i, gsl_vector_get(x, i));
  }
  printf("cond(J) = %8.4f\n", 1.0 / rcond);
}

void run_gsl_multifit_nlinear_df(
    opt_function func_f,
    opt_function* func_dfs,
    double* params,
    double* parerr,
    size_t params_len,
    double* ts,
    double* ys,
    size_t vars_len,
    double* args,
    size_t args_len,
    size_t max_iters,
    int* status
) {

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    gsl_vector *f;
    gsl_matrix *jacobian;
    gsl_matrix *covar = gsl_matrix_alloc(params_len, params_len);
  
    struct data d = {
       func_f,
       func_dfs,
       params_len,
       ts,
       ys,
       vars_len,
       args,
       args_len
    };

    gsl_vector_view x = gsl_vector_view_array(params, params_len);
    double dof = vars_len - params_len;
    double chisq, chisq0, chisq_dof;

    void* call_dfs_ptr = NULL;
    if (func_dfs != NULL) {
        call_dfs_ptr = &call_dfs;
    }

    const double xtol = 1e-08;
    const double gtol = 1e-08;
    const double ftol = 1e-08;

    fdf.f = call_f;
    fdf.df = call_dfs_ptr;
    fdf.fvv = NULL;
    fdf.n = vars_len;
    fdf.p = params_len;
    fdf.params = &d;

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc(T, &fdf_params, vars_len, params_len);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_init (&x.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system within a maximum of max_iters iterations */
    for (size_t i = 0; i < max_iters; i++) {

        double rcond;

        gsl_multifit_nlinear_iterate(w);
        gsl_multifit_nlinear_test(xtol, gtol, ftol, status, w);

        if (*status != 0) {

            if (*status == 1) {
                printf("Reason for stopping: Small Step Size\n");
            } else if (*status == 2) {
                printf("Reason for stopping: Small Gradient\n");
            }

            break;
        }

        gsl_multifit_nlinear_rcond(&rcond, w);
        if (isnan(rcond) && i != 0) {

            printf("Reason for stopping: Invalid Status\n");
            *status = -1;

            break;
        }

        //solver_progress(i, params_len, w);
    }

    /* compute covariance of best fit parameters */
    jacobian = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (jacobian, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

    for (size_t i = 0; i < params_len; i++) {
        chisq_dof = GSL_MAX_DBL(1, sqrt(chisq / dof));
        params[i] = gsl_vector_get(w->x, i);
        parerr[i] = chisq_dof * sqrt(gsl_matrix_get(covar, i, i));
    }

    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
}

void run_gsl_multifit_nlinear(
    opt_function func_f,
    double* params,
    double* parerr,
    size_t params_len,
    double* ts,
    double* ys,
    size_t vars_len,
    double* args,
    size_t args_len,
    size_t max_iters,
    int* status
) {
    run_gsl_multifit_nlinear_df(
        func_f,
        NULL,
        params,
        parerr,
        params_len,
        ts,
        ys,
        vars_len,
        args,
        args_len,
        max_iters,
        status
    );
}

