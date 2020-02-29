'''
Weak and strong simulation example plot 
Euler -Maruyama approximations 
Example scalar linear SDE
Parameter values
'''

import numpy as np
import pylab as pl


# parameters
a = 3.0         # drift coefficient
b = 1.4         # diffusion coefficient is b
y0 = 1.0        # initial data

P = 10          # total # of sample paths
h = 0.01        # stepsize
T = 1.0         # global time interval
N = int(T / h)  # number of subintervals
PLOT = True
colors = pl.cm.Reds(np.linspace(0, 1, P + 1))

fig, ax = pl.subplots(2)

# Binomial branching process increments
dw = np.zeros((N, P))
w = np.zeros((N + 1, P))
binom = np.random.binomial(1, 1 / 2, size=(N, P))  # Gives 0 or 1 with prob 1/2 %
dw = np.sqrt(h) * (1.0 - 2.0 * binom)              # Binomial increments dw
w[1: N + 1, :] = np.cumsum(dw, axis=1)             # Binomial paths themselves

# Weak solution by Euler-Maruyama approximation
y = np.zeros((N + 1, P))
for p in range(P):
    y[0, p] = y0
    for n in range(N):
        y[n + 1, p] = y[n, p] + a * y[n, p] * h + b * y[n, p] * dw[n, p]

if PLOT:
    for i in range(P):
        ax[0].plot(y[:, i], label=i, color=colors[i])
    ax[0].legend()

# Approximate Wiener path increments
dW = np.zeros((N, P))
W = np.zeros((N + 1, P))
dW = np.sqrt(h) * np.random.normal(size=(N, P))
W[1: N + 1, :] = np.cumsum(dW, axis=1)

# Strong solution by Euler-Maruyama approximation
Y = np.zeros((N + 1, P))
for p in range(P):
    Y[0, p] = y0
    for n in range(N):
        Y[n + 1, p] = Y[n, p] + a * Y[n, p] * h + b * Y[n, p] * dW[n, p]

if PLOT:
    for i in range(P):
        ax[1].plot(Y[:, i], label=i, color=colors[i])
    ax[1].legend()
    pl.show()

# Compute the expectations of y and Y at each timestep
expect_y = np.zeros(N + 1)
expect_Y = np.zeros(N + 1)

for n in range(N + 1):
    expect_Y[n] = np.mean(y[n, :])
    expect_Y[n] = np.mean(Y[n, :])
