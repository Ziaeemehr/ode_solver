# file: runme.py

# This file illustrates the proxy class C++ interface generated
# by SWIG.

import odesolver
import pylab as pl
import numpy as np
from time import time
from numpy import pi


dt = 0.01
ode = odesolver.HO(0.05)
x0 = [0.5, 1.0]
s =  odesolver.integrate_euler(ode, x0, 0., 100., dt)
