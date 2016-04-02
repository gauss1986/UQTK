#ifndef __KL_H__
#define __KL_H__

struct KL_params {
    double a, c;
};

/* Function of even terms*/
double even (double x, void *params);

/* Function of odd terms*/
double odd (double x, void *params);

/* Returns the derivative of the even fn. at x */
double even_deriv (double x, void *params);

/* Returns the derivative of the odd fn. at x */
double odd_deriv (double x, void *params);

/* Evaluates both the quadratic and its derivative at x and stores the
 *  * value in *y and the derivative in *dy. */
void even_fdf (double x, void *params, double *y, double *dy);

/* Evaluates both the quadratic and its derivative at x and stores the
 *  * value in *y and the derivative in *dy. */
void odd_fdf (double x, void *params, double *y, double *dy);

#endif
