#include "odesolver.hpp"


dim2 integrate_euler(
    std::function<void(const dim1 &, dim1 &, const double)> func,
    dim1 &y0,
    double ti,
    double tf,
    double dt)
{
    size_t n = y0.size();
    dim2 state;
    dim1 dydt(n);
    size_t nsteps = int((tf - ti) / dt);

    state.reserve(nsteps);
    for (size_t i = 0; i < nsteps; i++)
        state[i].reserve(n);

    state.push_back(y0);

    for (size_t it = 1; it < nsteps; ++it)
    {
        double t = it * dt;
        func(y0, dydt, t);
        for (size_t i = 0; i < n; ++i)
        {
            y0[i] += dydt[i] * dt;
        }

        state.push_back(y0);
    }

    return state;
}

// int main()
// {
//     double dt = 0.01;
//     HO ho(0.05);
//     dim1 x0{0.5, 1.0};
//     auto s = integrate_euler(ho, x0, 0., 100., dt);
//     // cout << s.size() << " " << s[0].size() << endl;

//     for (size_t i = 0; i < s.size(); ++i)
//     {
//         cout << i * dt << " ";
//         for (size_t j = 0; j < s[0].size(); ++j)
//             cout << s[i][j] << " ";
//         cout << endl;
//     }

//     return 0;
// }
