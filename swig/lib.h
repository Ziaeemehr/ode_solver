#include <iostream>
#include <cstdio>
#include <vector>
#include <omp.h>
#include <cmath>

class ODE;

typedef std::vector<double> dim1;
typedef std::vector<std::vector<double>> dim2;
typedef dim1 (ODE::*ode_sys)(const dim1 &);

class my_class 
{
    private:    
        int N;
    public:
        my_class(){ }
        std::vector<double> half(const std::vector<double>& ); 
};

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
    void eulerIntegrator(dim1 &y, ode_sys dydt)
    {
        dim1 f(N);
        f = (this->*dydt)(y);
        for (int i = 0; i < y.size(); i++)
            y[i] += f[i] * dt;
    }
};