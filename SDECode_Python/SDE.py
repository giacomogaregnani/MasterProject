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

    def grad_v(self, x):
        return float()

    def f(self, x, p):
        return -1.0 * p[0] * self.grad_v(x)

    def g(self, x, p):
        return np.sqrt(2.0 * p[np.size(p)-1])


class OrnUhl(Langevin):
    def __init__(self):
        Langevin.__init__(self)

    def v(self, x):
        return x*x / 2

    def grad_v(self, x):
        return x


class Bistable(Langevin):
    def __init__(self):
        Langevin.__init__(self)

    def v(self, x):
        x_sqd = x * x
        return x_sqd * x_sqd / 4.0 - x_sqd / 2

    def grad_v(self, x):
        return x * x * x - x


class MSLangevin(SDE):
    def __init__(self, eps):
        SDE.__init__(self)
        self.eps = eps

    def v(self, x):
        return float()

    def grad_v(self, x):
        return float()

    def p(self, x):
        return float()

    def grad_p(self, x):
        return float()

    def f(self, x, p):
        return -1.0 * p[0] * self.grad_v(x) - 1.0 / self.eps * self.grad_p(x / self.eps)

    def g(self, x, p):
        return np.sqrt(2.0 * p[np.size(p)-1])


class MSOrnUhl(MSLangevin):
    def __init__(self, eps):
        MSLangevin.__init__(self, eps)

    def v(self, x):
        return x*x / 2.0

    def grad_v(self, x):
        return x

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)


class MSBistable(MSLangevin):
    def __init__(self, eps):
        MSLangevin.__init__(self, eps)

    def v(self, x):
        x_sqd = x * x
        return x_sqd * x_sqd / 4.0 - x_sqd / 2

    def grad_v(self, x):
        return x * x * x - x

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)
