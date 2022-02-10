from solvers import SOLVERS
import numpy as np
import pylab as pl
from copy import copy
from scipy.integrate import odeint

def ode_sys(x0, t):
    x, y = x0
    xp = x - x * y - 0.1 * x * x
    yp = x*y - y - 0.05 * y * y
    return np.asarray([xp, yp])

def integrate(filename, t0, y0, method):
    t = t0
    x0 = copy(y0)
    ode = SOLVERS(dt, ode_sys, method)
    with open(filename, "w") as ofile:
        for i in range(num_steps):
            t = t + dt
            x0 = ode.do_step(x0, t)
            
            ofile.write("{:18.6f}".format(t))
            for j in x0:
                ofile.write("{:25.9f} ".format(j))
            ofile.write("\n")


# ------------------------------------------------------------------#

if __name__ == "__main__":

    dt = 0.01
    init_x0 = [2.0, 1.0]
    t0 = 0.0
    dimension = 2
    simulation_time = 30.0
    num_steps = int((simulation_time - t0) / dt)

    ### ODEINT ------------------------------------------------------
    x0 = copy(init_x0)
    t = np.arange(t0, simulation_time, dt)
    sol = odeint(ode_sys, x0, t)
    np.savetxt("data/odeint.txt",
               np.vstack((t, sol[:, 0], sol[:, 1])).T, 
               fmt="%18.9f")

    # --------------------------------------------------------------#
    integrate("data/bg.txt", t0, init_x0, "bg")
    integrate("data/heun.txt", t0, init_x0, "heun")
    integrate("data/euler.txt", t0, init_x0, "euler")

    bg = np.loadtxt("data/bg.txt")
    heun = np.loadtxt("data/heun.txt")
    euler = np.loadtxt("data/euler.txt")

    pl.plot(sol[::2, 0], sol[::2, 1], "r.", ms=2, alpha=0.5,label="odeint")
    pl.plot(euler[:, 1], euler[:, 2], label="euler", lw=1, alpha=0.5)
    pl.plot(bg[::2, 1], bg[::2, 2], "b*", ms=2, alpha=0.5, label="Bogacki-Shampine")
    pl.plot(heun[::2, 1], heun[::2, 2], "k.:", ms=2, alpha=0.5, label="heun")

    pl.legend()
    pl.savefig("data/fig.png", dpi=150)
    pl.show()
