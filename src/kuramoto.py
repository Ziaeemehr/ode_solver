import numpy as np
import pylab as pl
import networkx as nx
from scipy.integrate import odeint


def kuramoto(x, t):
    
    dxdt = np.zeros(N)

    for i in range(N):
        sumj = np.sum(adj[i, :] * np.sin(x-x[i]))
        dxdt[i] = omega[i] + coupling * sumj
    return dxdt
    


N = 10
coupling = 0.01
omega = [1.0] * N
t = np.arange(0, 100, 0.05)
x0 = np.random.uniform(-np.pi, np.pi, N)
adj = nx.to_numpy_array(nx.complete_graph(N))

theta = odeint(kuramoto, x0, t)
print(theta.shape)
# evolution of order parameter
r = np.zeros(len(t))
for i in range(len(t)):
    r[i] = abs(sum(np.exp(1j * theta[i]))) / N
    
pl.plot(t, r, lw=2)
pl.show()
