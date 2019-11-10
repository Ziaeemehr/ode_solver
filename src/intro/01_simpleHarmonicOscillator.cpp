#include <iostream>
#include <assert.h>
#include <random>
#include <vector>

typedef std::vector<double> dim1;
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
    for (int i = 0; i < n; ++i)
        y[i] += dxdt[i] * dt;
}

void rungeKutta4Integrator(dim1 &y, const double dt)
{
    int n = y.size();
    dim1 k1(n);
    dim1 k2(n);
    dim1 k3(n);
    dim1 k4(n);
    dim1 f(n);

    k1 = harmonicOscillator(y);
    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k1[i];

    k2 = harmonicOscillator(f);
    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k2[i];
    
    k3 = harmonicOscillator(f);
    for (int i = 0; i < n; i++)
        f[i] = y[i] + dt * k3[i];
    
    k4 = harmonicOscillator(f);
    
    for (int i = 0; i < n; i++)
        y[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dt / 6.0;
}
//-------------------------------------------------------//

int main(int argc, char **argv)
{
    const int N = 2;
    const double dt = 0.01;
    const double tSimulation = 100.0;
    const size_t steps = int(tSimulation / dt);

    FILE *outFile = fopen("resutl.txt", "w");

    dim1 x{0.0, 1.0};

    for (size_t i = 0; i < steps; ++i)
    {
        // eulerIntegrator(x, dt);
        rungeKutta4Integrator(x, dt);

        fprintf(outFile, "%15.6f", (i * dt));
        for (int j = 0; j < N; ++j)
            fprintf(outFile, "%15.9f", x[j]);
        fprintf(outFile, "\n");
    }

    fclose(outFile);
}
