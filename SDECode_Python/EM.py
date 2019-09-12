from SDE import *


class EulerMaruyama:
    def __init__(self, sde):
        self.sde = sde

    def solve(self, ic, n, t, p):
        y = np.zeros(n+1)
        y[0] = ic
        h = t/n
        dw = np.multiply(np.sqrt(h), np.random.normal(0, 1, n))
        for i in range(0, n):
            y[i+1] = y[i] + self.sde.f(y[i], p)*h + self.sde.g(y[i], p)*dw[i]
        return y
