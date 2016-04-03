#include <stdio.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "KL_fn.h"
#define PI 3.14159265
#define D 1.e-7
#define ERR 1.e-7

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

int
main (void){
  double a=5;
  double clen=0.1;
  double c=1/clen;
  int N = 10;

  int status;
  int max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r_expected = 0.0;
  struct KL_params params = {a,c};

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  printf ("using %s method\n", gsl_root_fsolver_name (s));

  double intv = PI/2/a;
  for (int i=0;i<N;i++){
    std::cout << "No." << i << std::endl;
    double x_lo=0, x_hi=0;

    gsl_function F;
    if (i%2 == 0){
        x_lo = std::max((i-1)*intv,0.0);
        x_hi = (i+1)*intv;
        F.function = &odd;
    }
    else {
        x_lo = std::max(i*intv,0.0);
        x_hi = (i+2)*intv;
        F.function = &even;
    }
    F.params = &params;

    gsl_root_fsolver_set (s, &F, x_lo+D, x_hi-D);

    printf ("%5s [%9s, %9s] %9s %10s %9s\n","iter", "lower", "upper", "root", "err", "err(est)");

    int iter = 0;
    double r = 0;
    do{
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,0, ERR);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
              iter, x_lo, x_hi,
              r, r - r_expected, 
              x_hi - x_lo);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

  }

  gsl_root_fsolver_free (s);
  return status;

}
