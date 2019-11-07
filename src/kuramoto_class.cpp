#include <iostream>
#include <algorithm>
#include <assert.h>
#include <numeric>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include "omp.h"

class ODE;

typedef std::vector<double> dim1;
typedef std::vector<std::vector<double>> dim2;
typedef dim1 (ODE::*ode_sys)(const dim1 &);
typedef void (ODE::*Integrator)(dim1&, ode_sys);
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
    // dim1 integrate(Integrator integrator)
    dim1 integrate(Integrator integrator, ode_sys ode_system)
    {
        dim1 r(n_steps);

        dim1 y = initial_phases;
        for (int step = 0; step < n_steps; ++step)
        {
            (this->*integrator)(y, ode_system);
            // (this->*integrator)(y, &ODE::kuramoto_model);
            // eulerIntegrator(y, &ODE::kuramoto_model);
            r[step] = order_parameter(y);
        }
        return r;
    }
    //---------------------------------------------------------------
    dim1 kuramoto_model(const dim1 &x)
    {
        dim1 dy(N);

        for (int i = 0; i < N; i++)
        {
            double sumx = 0;
            for (int j = 0; j < N; j++)
                if (i != j)
                    sumx += sin(x[j] - x[i]);

            dy[i] = omega[i] + coupling * sumx;
        }

        return dy;
    }
    //---------------------------------------------------------------
    void eulerIntegrator(dim1 &y, ode_sys dydt)
    {
        dim1 f(N);
        f = (this->*dydt)(y);
        for (int i = 0; i < y.size(); i++)
            y[i] += f[i] * dt;
    }
    //---------------------------------------------------------------
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
    dim1 initial_phases = {
        0.1, 0.3, -0.5, 1.4, -1.5,
        0.1, 0.6, 0.0, 0.2, -0.5};
    dim1 omega(N, 1.0);
    ODE ode;
    ode.set_params(
        N,     //number of nodes
        0.01,  //dt
        100.0, //t_simulation
        0.0,   //t_transition
        0.01,  //coupling
        initial_phases,
        omega);
    // dim1 r = ode.integrate(&ODE::eulerIntegrator);
    dim1 r = ode.integrate(&ODE::eulerIntegrator, &ODE::kuramoto_model);

    for (auto i : r)
        std::cout << i << std::endl;

    return 0;
}