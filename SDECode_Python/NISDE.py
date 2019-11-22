from SDE import *


class EulerMaruyama:
    def __init__(self, sde):
        self.sde = sde

    def solve(self, ic, n, t, p, seed=None):
        local_state = np.random.RandomState(seed)
        y = np.zeros(n+1)
        y[0] = ic
        h = t/n
        dw = np.multiply(np.sqrt(h), local_state.normal(0, 1, n))
        for i in range(0, n):
            y[i+1] = y[i] + self.sde.f(y[i], p)*h + self.sde.g(y[i], p)*dw[i]
        return y, dw


class Milstein:
    def __init__(self, sde):
        self.sde = sde

    def solve(self, ic, n, t, p):
        y = np.zeros(n+1)
        y[0] = ic
        h = t/n
        sqrt_h = np.sqrt(h)
        dw = np.multiply(sqrt_h, np.random.normal(0, 1, n))
        for i in range(0, n):
            f = self.sde.f(y[i], p)
            g = self.sde.g(y[i], p)
            z = y[i] + f*h + g*sqrt_h
            y[i+1] = y[i] + f*h + g*dw[i] + (self.sde.g(z, p) - g) * (dw[i]*dw[i]-h) / (2.0 * sqrt_h)
        return y