#ifndef __KL_w__
#define __KL_w__

struct KL_params {
    double a, c;
};

/* Function of even terms*/
double even (double x, void *params);

/* Function of odd terms*/
double odd (double x, void *params);

/* Function to find w values of KL eigensystems*/
Array1D<double> KL_exp (double sigma, double a, double clen, int N, Array2D<double>& KLmodes, Array1D<double>& xgrid);

/* Compute the eighenvalues*/
double KL_lambda(double sigma, double c, double w);
#endif
