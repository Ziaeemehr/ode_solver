import lib

c = lib.my_class()

print c.half([10, 20.0, 30.0])


ode = lib.ODE()
ode.set_params(2, 0.1, 10.0, 5.0, 0.1, [1.0, 1.0], [0.1, -0.1])