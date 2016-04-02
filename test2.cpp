#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "demo_fn.h"

#define MAX_ITERATIONS 100

int main()
{
    int status;
    int i;
    const gsl_root_fsolver_type *solver_type;
    gsl_root_fsolver *solver;
    gsl_function F;
    double r;

    /* We want to solve x^2 - 5 */
    struct quadratic_params params = { 1.0, 0.0, -5.0 };
    double r_expected = sqrt(5.0);
    /* On the interval [0, 5] */
    double x_lo = 0.0, x_hi = 5.0;

    /* Set up the function */
    F.function = &quadratic;
    F.params = &params;

    /* Allocate a bisection solver and set it to use F */
    solver_type = gsl_root_fsolver_bisection;
    solver = gsl_root_fsolver_alloc(solver_type);
    gsl_root_fsolver_set(solver, &F, x_lo, x_hi);

    printf("using %s method\n", gsl_root_fsolver_name(solver));

    printf("%5s [%9s, %9s] %9s %10s %9s\n",
           "iter", "lower", "upper", "root", "err", "err(est)");

    status = GSL_CONTINUE;
    for (i = 1; i <= MAX_ITERATIONS && status == GSL_CONTINUE; ++i) {
        /* iterate one step of the solver */
        status = gsl_root_fsolver_iterate(solver);
        if (status != GSL_SUCCESS)
            break;

        /* get the solver's current best solution and bounds */
        r = gsl_root_fsolver_root(solver);
        x_lo = gsl_root_fsolver_x_lower(solver);
        x_hi = gsl_root_fsolver_x_upper(solver);

        /* Check to see if the solution is within 0.001 */
        status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
        if (status == GSL_SUCCESS)
            printf("Converged:\n");

        printf("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
               i, x_lo, x_hi, r, r - r_expected, x_hi - x_lo);
    }

    /* Free the solver */
    gsl_root_fsolver_free(solver);

    if (status == GSL_CONTINUE) {
        printf("error: too many iterations");
    } else if (status != GSL_SUCCESS) {
        printf("error: %s\n", gsl_strerror(status));
    }

    return status;
}
