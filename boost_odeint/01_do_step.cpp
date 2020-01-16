
#include <vector>
#include <iostream>
#include <boost/numeric/odeint.hpp>

/*-----------------------------------------------------------------*/
/* rhs_function, The type of container used to hold the state vector */
typedef std::vector<double> dim1;
typedef std::vector<dim1> dim2;

/* The rhs of x' = f(x) */
void harmonic_oscillator(const dim1 &x, dim1 &dxdt, const double /* t */)
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - 0.05 * x[1];
}
/*-----------------------------------------------------------------*/

// to record the state of system in each time step
struct push_back_state_and_time
{
    dim2 &m_states;
    dim1 &m_times;

    push_back_state_and_time(dim2 &states, dim1 &times)
        : m_states(states), m_times(times) {}

    void operator()(const dim1 &x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

/*-----------------------------------------------------------------*/
int main(int /* argc */, char ** /* argv */)
{
    using namespace std;
    using namespace boost::numeric::odeint;

    const double dt = 0.01;
    const double tSimulation = 100.0;

    // state_initialization
    dim1 x {0.0, 1.0};

    // containers for time and state of the system
    dim2 x_vec;
    dim1 times;

    //[ define_const_stepper
    euler<dim1> stepper;
    // runge_kutta4<dim1> stepper;

    // integration
    double t = 0.0;
    long unsigned numSteps = (int)(tSimulation / dt);
    x_vec.resize(numSteps);
    for (size_t i = 0; i < numSteps; ++i)
    {
        stepper.do_step(harmonic_oscillator, x, t, dt);
        t += dt;

        times.push_back(t);
        for (int j = 0; j < 2; j++)
            x_vec[i].push_back(x[j]);
    }

    for (size_t i = 0; i < numSteps; i++)
        printf("%15.6f%15.9f%15.9f\n",
               times[i],
               x_vec[i][0],
               x_vec[i][1]);

    return 0;
}
/*-----------------------------------------------------------------*/
