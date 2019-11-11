import numpy as np
from scipy.integrate import odeint


def f(x, t):
    dxdt = [x[1], - x[0] - 0.05 * x[1]]
    return dxdt

def g(x, t):
    dxdt = 1+x*x
    return dxdt

# x = 0.0
x = [0.0, 1.0]
t = np.arange(0, 100, 0.01)
sol = odeint(f, x, t)
print (t[-1], sol[-1, :])
