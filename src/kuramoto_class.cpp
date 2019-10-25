#include <iostream>
#include <algorithm>
#include <assert.h>
#include <numeric>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include "omp.h"

typedef std::vector<double> dim1;
typedef std::vector<std::vector<double>> dim2;
//-------------------------------------------------------------------

class ODE
{
private:
    int N;
    double dt;
    dim1 omega;
    int n_steps;
    double t_sim;
    double t_cut;
    double coupling;
    dim1 initial_phases;

public:
    virtual ~ODE() {}
    //---------------------------------------------------------------
    void set_params(int N,
                    double dt,
                    double t_sim,
                    double t_cut,
                    double coupling,
                    dim1 initial_phase,
                    dim1 omega)
    {
        this->N = N;
        this->dt = dt;
        this->t_sim = t_sim;
        this->t_cut = t_cut;
        this->coupling = coupling;
        this->initial_phases = initial_phase;
        this->omega = omega;
        n_steps = (int)(t_sim / dt);
    }
    //---------------------------------------------------------------
    dim1 integrate()
    {
        dim1 r(n_steps);

        dim1 y = initial_phases;
        for (int step = 0; step < n_steps; ++step)
        {
            euler_integrator(y);
            r[step] = order_parameter(y);
        }
        return r;
    }
    //---------------------------------------------------------------
    dim1 kuramoto_model(const dim1 &x)
    {
        dim1 dydt(N);

        for (int i = 0; i < N; i++)
        {
            double sumx = 0;
            for (int j = 0; j < N; j++)
                if (i != j)
                    sumx += sin(x[j] - x[i]);

            dydt[i] = omega[i] + coupling * sumx;
        }

        return dydt;
    }
    //---------------------------------------------------------------
    void euler_integrator(dim1 &y)
    {
        dim1 f(N);
        f = kuramoto_model(y);
        for (int i = 0; i < y.size(); i++)
            y[i] += f[i] * dt;
    }
    //---------------------------------------------------------------
    double order_parameter(const dim1 &x)
    {
        int n = x.size();
        double real_R = 0.;
        double imag_R = 0.;

        for (int i = 0; i < n; i++)
        {
            real_R += cos(x[i]);
            imag_R += sin(x[i]);
        }
        real_R /= (double)n;
        imag_R /= (double)n;
        double r = sqrt(real_R * real_R + imag_R * imag_R);

        return r;
    }
};

int main()
{
    size_t N = 10;
    dim1 initial_phases = {0.1, 0.3, -0.5, 1.4, -1.5, 0.1, 0.6, 0.0, 0.2, -0.5};
    dim1 omega(N, 1.0);
    ODE ode;
    ode.set_params(N, 0.01, 100.0, 0.0, 0.01, initial_phases, omega);
    dim1 r = ode.integrate();

    for (auto i : r)
        std::cout << i << std::endl;

    return 0;
}