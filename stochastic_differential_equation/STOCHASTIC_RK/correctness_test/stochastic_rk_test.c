#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "stochastic_rk.h"

int main(void);
void test01(void);
double fi(double x);
double gi(double x);

/******************************************************************************/

int main(void)
{
    test01();

    return 0;
}
/******************************************************************************/

void test01(void)
{
    double h;
    int i;
    int n;
    double q;
    int seed;
    double t;
    double t0 = 0.0;
    double tn = 1.0;
    double *x;

    n = 1000;
    x = (double *)malloc((n + 1) * sizeof(double));
    h = (tn - t0) / (double)(n);
    q = 1.0;
    seed = 123456789;

    i = 0;
    t = t0;
    x[i] = 0.0;

    fprintf(stdout, "  %8d  %14f  %14f\n", i, t, x[i]);

    for (i = 1; i <= n; i++)
    {
        t = ((double)(n - i) * t0 + (double)(i)*tn) / (double)(n);

        x[i] = rk2_ti_step(x[i - 1], t, h, q, fi, gi, &seed);

        fprintf(stdout, "  %8d  %14f  %14f\n", i, t, x[i]);
    }

    free(x);
    return;
}
/******************************************************************************/

double fi(double x)

/******************************************************************************/
/*
  Purpose:
    FI is a time invariant deterministic right hand side.

  Parameters:
    Input, double X, the argument.
    Output, double FI, the value.
*/
{
    double value;

    value = 3.0;

    return value;
}
/******************************************************************************/

double gi(double x)

/******************************************************************************/
/*
  Purpose:
    GI is a time invariant stochastic right hand side.

  Parameters:
    Input, double X, the argument.
    Output, double GI, the value.
*/
{
    double value;

    value = 1.4;

    return value;
}
