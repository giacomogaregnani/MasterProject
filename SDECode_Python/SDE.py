import numpy as np


class SDE:
    def __init__(self):
        pass

    def f(self, x, p):
        return float()

    def g(self, x, p):
        return float()


class Langevin(SDE):
    def __init__(self):
        SDE.__init__(self)

    def v(self, x):
        return float()

    def gradv(self, x):
        return float()

    def f(self, x, p):
        return -1.0 * p[0] * self.gradv(x)

    def g(self, x, p):
        return np.sqrt(2.0 * p[np.size(p)-1])


class OrnUhl(Langevin):
    def __init__(self):
        Langevin.__init__(self)

    def v(self, x):
        return x*x / 2

    def gradv(self, x):
        return x
