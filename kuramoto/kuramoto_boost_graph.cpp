#include "omp.h"
#include <cmath>
#include <time.h>
#include <chrono>
#include <vector>
#include <random>
#include <complex>
#include <numeric>
#include <utility>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <sys/time.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>

typedef boost::minstd_rand base_generator_type;
typedef std::vector<double> dim1;
typedef std::vector<std::vector<double>> dim2;
typedef std::vector<int> dim1I;
typedef std::vector<dim1I> dim2I;
typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS>
    Graph;

typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typename boost::graph_traits<Graph>::adjacency_iterator ai, ai_end;

using std::cout;
using std::endl;

//-------------------------------------------------------------------
void graph_from_matrix(dim2I &adj, Graph &graph_t)
{
    const size_t N = adj.size();
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i + 1; j < N; j++)
        {
            if (adj[i][j] != 0)
                boost::add_edge(i, j, graph_t);
        }
    }
}

double get_wall_time()
{
    /*measure real passed time */
    struct timeval time;
    if (gettimeofday(&time, NULL))
    {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

class ODE
{
private:
    int N;
    dim2I adj;
    double dt;
    dim1 omega;
    Graph graph_t;
    int n_steps;
    double t_sim;
    double t_cut;
    double coupling;
    dim1 initial_phases;
    IndexMap vertex_id;
    std::pair<vertex_iter, vertex_iter> vp;

public:
    virtual ~ODE() {}
    //---------------------------------------------------------------
    void set_params(int N,
                    double dt,
                    double t_sim,
                    double t_cut,
                    double coupling,
                    dim1 initial_phase,
                    dim1 omega,
                    Graph graph_t,
                    dim2I adj)
    {
        this->N = N;
        this->dt = dt;
        this->t_sim = t_sim;
        this->t_cut = t_cut;
        this->coupling = coupling;
        this->initial_phases = initial_phase;
        this->omega = omega;
        this->graph_t = graph_t;
        this->adj = adj;
        n_steps = (int)(t_sim / dt);

        // cout << boost::num_vertices(graph_t) << endl;
        vertex_id = get(boost::vertex_index, graph_t);
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
                if ((i != j) && adj[i][j] != 0)
                    sumx += sin(x[j] - x[i]);

            dydt[i] = omega[i] + coupling * sumx;
        }

        return dydt;
    }
    //---------------------------------------------------------------
    dim1 kuramoto_model_bgl(const dim1 &x)
    {
        dim1 dydt(N);

        typename boost::graph_traits<Graph>::adjacency_iterator ai, ai_end;
    
        for (vp = vertices(graph_t); vp.first != vp.second; ++vp.first)
        {
            Vertex v = *vp.first;
            int i = vertex_id[v];
            double sumj = 0.0;
            for (boost::tie(ai, ai_end) = adjacent_vertices(v, graph_t); ai != ai_end; ++ai)
            {
                int j = get(vertex_id, *ai);
                sumj += sin(x[j] - x[i]);
            }
            dydt[i] = omega[i] + coupling * sumj;
        }

        return dydt;
    }
    //---------------------------------------------------------------
    void euler_integrator(dim1 &y)
    {
        dim1 f(N);
        f = kuramoto_model_bgl(y);
        // f = kuramoto_model(y);
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
    dim2I adj = {
        {0, 1, 0, 1, 0, 1, 0, 1, 1, 1},
        {1, 0, 0, 1, 1, 1, 0, 1, 0, 1},
        {0, 0, 0, 1, 1, 1, 1, 0, 1, 1},
        {1, 1, 1, 0, 1, 1, 1, 0, 0, 1},
        {0, 1, 1, 1, 0, 0, 1, 0, 1, 1},
        {1, 1, 1, 1, 0, 0, 0, 1, 0, 1},
        {0, 0, 1, 1, 1, 0, 0, 0, 0, 1},
        {1, 1, 0, 0, 0, 1, 0, 0, 1, 1},
        {1, 0, 1, 0, 1, 0, 0, 1, 0, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 0}};

    ODE ode;
    Graph graph_t(N);
    graph_from_matrix(adj, graph_t);
    ode.set_params(N, 0.001, 500.0, 0.0, 0.0005, initial_phases, omega, graph_t, adj);

    double start = get_wall_time();

    dim1 r = ode.integrate();

    cout << "# Done in : " << get_wall_time() - start << endl;

    // for (auto i : r)
    //     std::cout << i << std::endl;
    std::cout << r[r.size() - 1] << std::endl;

    return 0;
}

// typename boost::graph_traits<Graph>::adjacency_iterator ai, ai_end;
//         for (vp = vertices(graph_t); vp.first != vp.second; ++vp.first)
//         {
//             Vertex v = *vp.first;
//             std::cout << vertex_id[v] << " : ";
//             for (boost::tie(ai, ai_end) = adjacent_vertices(v, graph_t); ai != ai_end; ++ai)
//                 std::cout << get(vertex_id, *ai) << " ";
//             std::cout << std::endl;
//         }