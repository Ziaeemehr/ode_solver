#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>

/*-----------------------------------------------------------------*/
/* rhs_function, The type of container used to hold the state vector */
typedef std::vector<double> dim1;
typedef std::vector<dim1> dim2;

class harm_oscillator
{
    double gamma;

public:
    harm_oscillator(double gam) : gamma(gam) {}

    void operator()(const dim1 &x, dim1 &dxdt, const double /* t */)
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - gamma * x[1];
    }
};
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

void write_stdout(const dim1 &x, const double t)
{
    printf("%15.6f%15.9f%15.9f\n", t, x[0], x[1]);
}

/*-----------------------------------------------------------------*/

int main(int /* argc */, char ** /* argv */)
{
    using namespace std;
    using namespace boost::numeric::odeint;

    const double dt = 0.01;
    const double tSimulation = 100.0;

    // state_initialization
    dim1 x{0.0, 1.0};

    // containers for time and state of the system
    dim2 x_vec;
    dim1 times;

    //[ define_const_stepper
    euler<dim1> stepper;
    // runge_kutta4<dim1> stepper;

    // integration
    size_t steps = integrate_const(stepper,
                                   harm_oscillator(0.05), // ode system
                                   x,                     //initial condition
                                   0.0,                   //initial time
                                   tSimulation,           // final time
                                   dt,                    //time step
                                   write_stdout           //observer
    );
    return 0;
}
/*-----------------------------------------------------------------*/
