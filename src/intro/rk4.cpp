#include <vector>
#include <valarray>
#include <iostream>
using namespace std;

typedef valarray<double> dim1;

dim1 ode(const dim1& x)
{
    dim1 dxdt(1);
    dxdt[0] = 1.0 + x[0] * x[0];
    return dxdt;
}


dim1 rk4(dim1 &x,
         const double h,
         dim1 (*f)(const dim1 &))
{
    dim1 k1 = h * f(x);
    dim1 k2 = h * f(x + 0.5 * k1);
    dim1 k3 = h * f(x + 0.5 * k2);
    dim1 k4 = h * f(x + k3);

    x = x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

int main()
{
    const int N = 1;
    const double dt = 0.1;
    const double tSimulation = 1.4;
    const size_t steps = int(tSimulation / dt);

    dim1 x{0.0};
    for (size_t i = 0; i < steps; ++i)
    {
        rk4(x, dt, ode);

        printf("%15.6f", ((i+1) * dt));
        for (int j = 0; j < N; ++j)
            printf("%15.9f", x[j]);
        printf("\n");
    }


    return 0;
}