#ifndef __KL_H__
#define __KL_H__

struct KL_params {
    double a, c;
};

/* Function of even terms*/
double even (double x, void *params);

/* Function of odd terms*/
double odd (double x, void *params);

#endif
