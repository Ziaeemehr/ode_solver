#ifndef ODESOLVER_HPP
#define ODESOLVER_HPP

#include <vector>
#include <iostream>
#include <functional>

// using std::cout;
// using std::endl;
using std::vector;

typedef std::vector<double> dim1;
typedef std::vector<vector<double>> dim2;

class HO
{
    double gamma;

public:
    HO(double gamma_) : gamma(gamma_) {}

    void operator()(const vector<double> &y, vector<double> &dydt, const double t)
    {
        dydt[0] = y[1];
        dydt[1] = -y[0] - gamma * y[1];
    }
};

dim2 integrate_euler(
    std::function<void(const dim1 &, dim1 &, const double)> func,
    dim1 &y0,
    double ti,
    double tf,
    double dt);

#endif