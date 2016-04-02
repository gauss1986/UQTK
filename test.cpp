#include <stdio.h>
#include <cmath>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "KL_fn.h"
#define PI 3.14159265
#define D 1.e-7

double
even (double x, void *params)
{//compute data associated with eq: c*tan(a*x)+x=0

  struct KL_params *p 
    = (struct KL_params *) params;

  double a = p->a;
  double c = p->c;

  double temp = tan(a*x);

  return (c*temp+x);  
}

double
odd (double x, void *params)
{//compute data associated with eq: c*tan(a*x)+x=0

  struct KL_params *p 
    = (struct KL_params *) params;

  double a = p->a;
  double c = p->c;

  double temp = tan(a*x);

  return (c-x*temp);  
}

double
even_deriv (double x, void *params)
{
  struct KL_params *p 
    = (struct KL_params *) params;

  double a = p->a;
  double c = p->c;

  return (a*c/pow(cos(a*x),2)+1);
}

double
odd_deriv (double x, void *params)
{
  struct KL_params *p 
    = (struct KL_params *) params;

  double a = p->a;

  double temp = a*x;

  return (-temp/pow(cos(temp),2)-tan(temp));
}

void
even_fdf (double x, void *params, 
               double *y, double *dy)
{
  struct KL_params *p 
    = (struct KL_params *) params;

  double a = p->a;
  double c = p->c;

  double temp = a*x;

  *y = c*tan(temp)+x;
  *dy = (a*c/pow(cos(a*x),2)+1);
}

void
odd_fdf (double x, void *params, 
               double *y, double *dy)
{
  struct KL_params *p 
    = (struct KL_params *) params;

  double a = p->a;
  double c = p->c;

  double temp = a*x;

  *y = -x*tan(temp)+c;
  *dy = -temp/pow(cos(temp),2)-tan(temp);
}

int
main (void){
  double a=5;
  double clen=0.01;
  double c=1/clen;
  double N = 10.0;

  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0, r_expected = 0.0;
  struct KL_params params = {a,c};

  for (int i=0;i<=ceil(N/2);i++){

  gsl_function F;
  F.function = &odd;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  printf ("using %s method\n", 
          gsl_root_fsolver_name (s));

  double x_lo = PI/2/a+D, x_hi = 3*PI/2/a-D;
  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
          "iter", "lower", "upper", "root", 
          "err", "err(est)");

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
                                       0, 0.001);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
              iter, x_lo, x_hi,
              r, r - r_expected, 
              x_hi - x_lo);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (s);
}
  /*Newton's method*/

  const gsl_root_fdfsolver_type *solver_type;
  gsl_root_fdfsolver *solver;
  gsl_function_fdf FDF;
  double x0;
  
  /* Starting from a guess of 5.0 */
  double x = 10;
  
  FDF.f = &odd;
  FDF.df = &odd_deriv;
  FDF.fdf = &odd_fdf;
  FDF.params = &params;

  /* Allocate a newton solver and set it to use FDF */
    solver_type = gsl_root_fdfsolver_newton;
    solver = gsl_root_fdfsolver_alloc(solver_type);
    gsl_root_fdfsolver_set(solver, &FDF, x);

    printf("using %s method\n", gsl_root_fdfsolver_name(solver));

    printf ("%-5s %10s %10s %10s\n", "iter", "root", "err", "err(est)");

    status = GSL_CONTINUE;
    for (int i = 1; i <= max_iter && status == GSL_CONTINUE; ++i) {
        /* Iterate one step of the solver */
        status = gsl_root_fdfsolver_iterate(solver);
        if (status != GSL_SUCCESS)
            break;

        /* Get the new approximate solution */
        x0 = x;
        x = gsl_root_fdfsolver_root(solver);

        /* Check to see if the solution is within 0.001 */
        status = gsl_root_test_delta(x, x0, 0, 0.001);
        if (status == GSL_SUCCESS)
            printf("Converged:\n");

        printf ("%5d %10.7f %+10.7f %10.7f\n", i, x, x - r_expected, x - x0);
    }

    /* Free the solver */
    gsl_root_fdfsolver_free(solver);

    if (status == GSL_CONTINUE) {
        printf("error: too many iterations");
    } else if (status != GSL_SUCCESS) {
        printf("error: %s\n", gsl_strerror(status));
    }

  return status;

}
