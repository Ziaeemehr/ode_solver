#include <iostream>
#include <assert.h>
#include <random>
#include <vector>

class harmonicOscillator;
typedef std::vector<double> dim1;
typedef dim1 (harmonicOscillator::*ode_sys)(const dim1 &);
typedef void (harmonicOscillator::*Integrator)(dim1 &, ode_sys);

//-------------------------------------------------------//
class harmonicOscillator
{

private:
    const int N;
    const double dt;
    const double gamma;
    dim1 state;
    double tSimulation;

public:
    // constructor
    harmonicOscillator(int N,
                       const double dt,
                       const double gamma,
                       const double tSimulation,
                       const dim1 &x0) : N(N),
                                         dt(dt),
                                         gamma(gamma)
    {
        this->tSimulation = tSimulation;
        state = x0;
    }
    //---------------------------------------------------//
    void integrate(Integrator integrator, ode_sys ode_system)
    {

        FILE *outFile = fopen("resutlClass.txt", "w");
        const size_t numberOfSteps = int(tSimulation / dt);

        for (int step = 0; step < numberOfSteps; ++step)
        {
            (this->*integrator)(state, ode_system);

            // print state to file
            fprintf(outFile, "%15.6f", (step * dt));
            for (int j = 0; j < N; ++j)
                fprintf(outFile, "%15.9f", state[j]);
            fprintf(outFile, "\n");
        }

        fclose(outFile);
    }
    //---------------------------------------------------//
    dim1 dampedOscillator(const dim1 &x)
    {
        dim1 dxdt(2);

        dxdt[0] = x[1];
        dxdt[1] = -x[0] - gamma * x[1];

        return dxdt;
    }
    //---------------------------------------------------//
    void eulerIntegrator(dim1 &y, ode_sys dydt)
    {
        dim1 f(N);
        f = (this->*dydt)(y);
        for (int i = 0; i < y.size(); i++)
            y[i] += f[i] * dt;
    }
    //---------------------------------------------------//
    void rungeKutta4Integrator(dim1 &y, ode_sys dydt)
    {
        int n = y.size();
        dim1 k1(n), k2(n), k3(n), k4(n);
        dim1 f(n);

        k1 = (this->*dydt)(y);
        for (int i = 0; i < n; i++)
            f[i] = y[i] + 0.5 * dt * k1[i];

        k2 = (this->*dydt)(f);

        for (int i = 0; i < n; i++)
            f[i] = y[i] + 0.5 * dt * k2[i];
        k3 = (this->*dydt)(f);

        for (int i = 0; i < n; i++)
            f[i] = y[i] + dt * k3[i];
        k4 = (this->*dydt)(f);

        for (int i = 0; i < n; i++)
            y[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dt / 6.0;
    }
};

//-------------------------------------------------------//

int main(int argc, char **argv)
{
    const int N = 2;
    const double dt = 0.01;
    const double tSimulation = 100.0;

    dim1 x0{0.0, 1.0};

    harmonicOscillator ho(N, dt, 0.05, tSimulation, x0);
    ho.integrate(&harmonicOscillator::rungeKutta4Integrator,
                 &harmonicOscillator::dampedOscillator);

    return 0;
}