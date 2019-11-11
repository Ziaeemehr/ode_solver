#include <iostream>
#include <assert.h>
#include <valarray>
#include <vector>

typedef std::valarray<double> dim1;
//-------------------------------------------------------//
dim1 harmonicOscillator(const dim1 &x)
{
    dim1 dxdt(2);

    dxdt[0] = x[1];
    dxdt[1] = -x[0] - 0.05 * x[1];

    return dxdt;
}
//-------------------------------------------------------//
void eulerIntegrator(dim1 &y, const double dt)
{
    size_t n = y.size();
    dim1 dxdt(n);

    dxdt = harmonicOscillator(y);
    y =  y + dxdt * dt;
}
//-------------------------------------------------------//
void rungeKutta4Integrator(dim1 &y, const double dt)
{
    int n = y.size();
    
    dim1 k1 = dt * harmonicOscillator(y); 
    dim1 k2 = dt * harmonicOscillator(y + 0.5*k1); 
    dim1 k3 = dt * harmonicOscillator(y + 0.5*k2); 
    dim1 k4 = dt * harmonicOscillator(y + k3); 
  
    // Update next value of y 
    y = y + (k1 + 2.0 * k2 + 2.0 * k3 + k4)/6.0;
}
//-------------------------------------------------------//
int main(int argc, char **argv)
{
    const int N = 2;
    const double dt = 0.01;
    const double tSimulation = 100.0;
    const size_t steps = int(tSimulation / dt);

    dim1 x{0.0, 1.0};
    for (size_t i = 0; i < steps; ++i)
    {
        // eulerIntegrator(x, dt);
        rungeKutta4Integrator(x, dt);

        printf("%15.6f", ((i+1) * dt));
        for (int j = 0; j < N; ++j)
            printf("%15.9f", x[j]);
        printf("\n");
    }
    
    return 0;
}
//-------------------------------------------------------//