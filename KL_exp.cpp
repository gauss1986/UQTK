#include <stdio.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include "arraytools.h"
#include "Array2D.h"
#include "Array1D.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "KL_exp.h"
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

Array1D<double> KL_exp (double sigma, double a, double clen, int N, Array2D<double>& KLmodes, Array1D<double>& xgrid){
  Array1D<double> eigs(N,0.e0);
  double c=1/clen;

  int status;
  int max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
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

    int iter = 0;
    double r = 0;
    do{
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,0, ERR);

      if (status == GSL_SUCCESS){
        printf ("iter=%5d, result=%.7f, range=%.7f\n", iter, r, x_hi-x_lo); 
      }  
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    // compute the eigenvalues using the w just found
    eigs(i)=KL_lambda(sigma,c,r);

    for (unsigned int it=0;it<KLmodes.XSize();it++){
        if (i%2 == 0){
           KLmodes(it,i)=cos(r*(xgrid(it)-a))/sqrt(a+sin(2*r*a)/2/r); 
        }
        else{
           KLmodes(it,i)=sin(r*(xgrid(it)-a))/sqrt(a-sin(2*r*a)/2/r); 
        }
    }
  }

  gsl_root_fsolver_free (s);
  
  return eigs;
}

double KL_lambda(double sigma, double c, double w){
    double lambda = 2*sigma*sigma*c/(w*w+c*c);
    return lambda;
}
