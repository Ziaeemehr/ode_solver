import numpy as np
import pylab as pl
from scipy.integrate import odeint


def ode(x0, t):
    x, y = x0
    xp = x - x * y - 0.1 * x * x
    yp = x*y - y - 0.05 * y * y
    return np.asarray([xp, yp])
# ------------------------------------------------------------------#


def euler(x0, t, h, f):

    x0 += h * f(x0, t)
    return x0
# ------------------------------------------------------------------#


def BogackiShampine(x0, t, h, f):

    k1 = f(x0, t)
    k2 = f(x0 + 0.5 * h * k1, t + 0.5 * h)
    k3 = f(x0 + 0.75 * h * k2, t + 0.75 * h)

    x0 += h / 9.0 * (2.0 * k1 + 3.0 * k2 + 4.0 * k3)

    return x0

# ------------------------------------------------------------------#


if __name__ == "__main__":

    dt = 0.01
    # dts = 0.1 * dt
    x0 = [2.0, 1.0]
    tInitial = 0.0
    tFinal = 30.0

    t = np.arange(tInitial, tFinal, dt)
    sol = odeint(ode, x0, t)
    np.savetxt("../../data/odeint.txt",
               np.vstack((t, sol[:, 0], sol[:, 1])).T, fmt="%18.9f")

    numSteps = int((tFinal - tInitial) / dt)
    t = tInitial
    x0 = [2.0, 1.0]
    # --------------------------------------------------------------#
    EULER_FILE = open("../../data/euler.txt", "w")
    for i in range(numSteps):
        t += dt
        x0 = euler(x0, t, dt, ode)
        EULER_FILE.write("%18.9f " % t)
        for j in x0:
            EULER_FILE.write("%18.9f " % j)
        EULER_FILE.write("\n")
    EULER_FILE.close()

    # --------------------------------------------------------------#

    BS_FILE = open("../../data/bs.txt", "w")
    t = tInitial
    x0 = [2.0, 1.0]
    for i in range(numSteps):
        t += dt
        x0 = BogackiShampine(x0, t, dt, ode)
        BS_FILE.write("%18.9f " % t)
        for j in x0:
            BS_FILE.write("%18.9f " % j)
        BS_FILE.write("\n")
    BS_FILE.close()

    # euler(ode, 0, 30, x0, 0.1*dt)
    bs = np.loadtxt("../../data/bs.txt")
    eul = np.loadtxt("../../data/euler.txt")

    pl.plot(sol[:, 0], sol[:, 1], label="odeint", ls=":")
    pl.plot(eul[:, 1], eul[:, 2], label="euler")
    pl.plot(bs[:, 1], bs[:, 2], label="Bogacki-Shampine", ls="--")

    pl.legend()
    pl.savefig("../../data/fig.png", dpi=150)
    pl.show()
