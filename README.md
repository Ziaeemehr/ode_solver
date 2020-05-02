# ode_solver

-  This repository provide simple examples using ordinary differential equation solvers in python, C/C++.

#### Camparing the results:


|method     |      time        |  x0           | x1        |
| :---      |  :---:           |   :---:       | :---:     |
|intro/euler      |      99.990000   |-0.069638927   |0.117734245|
|intro/rk4        |      99.990000   |-0.044473790   |0.070138180|
|scipy.integrate.odeint    |  99.990000   | -0.04447404   |  0.07013844|
| boost euler  |   99.990000  | -0.070809774  |  0.117084690|
|boost rk4    |   99.990000  | -0.044473790  | 0.070138180|

