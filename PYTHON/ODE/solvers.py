
class SOLVERS:

    def __init__(self, dt, derivative, method="euler"):
        self.dt = dt
        self.derivative = derivative

        if method == "euler":
            self.integrator = self.euler
        elif method == "bg":
            self.integrator = self.BogackiShampine
        elif method == "heun":
            self.integrator = self.heun

    def do_step(self, x0, t):
        return self.integrator(x0, t, self.derivative)

    def euler(self, x0, t, derivative):
        x0 += self.dt * derivative(x0, t)
        return x0

    def heun(self, x0, t, derivative):

        k1 = derivative(x0, t)
        k2 = derivative(x0 + self.dt * k1, t + self.dt)
        x0 += 0.5*self.dt * (k1 + k2)

        return x0

    def BogackiShampine(self, x0, t, derivative):

        k1 = derivative(x0, t)
        k2 = derivative(x0 + 0.5 * self.dt * k1, t + 0.5 * self.dt)
        k3 = derivative(x0 + 0.75 * self.dt * k2, t + 0.75 * self.dt)
        x0 += self.dt / 9.0 * (2.0 * k1 + 3.0 * k2 + 4.0 * k3)

        return x0
