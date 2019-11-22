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

    def grad_v_vect(self, x):
        pass

    def lapl_v_vect(self, x):
        pass

    def f(self, x, p):
        return -1.0 * p[0] * self.grad_v(x)

    def g(self, x, p):
        return np.sqrt(2.0 * p[-1])


class OrnUhl(Langevin):
    def __init__(self):
        Langevin.__init__(self)

    def v(self, x):
        return x*x / 2

    def grad_v(self, x):
        return x

    def lapl_v_vect(self, x):
        return np.ones(np.size(x))

    def grad_v_vect(self, x):
        return x


class Quartic(Langevin):
    def __init__(self):
        Langevin.__init__(self)

    def v(self, x):
        x_sqd = x*x
        return x_sqd*x_sqd / 4

    def grad_v(self, x):
        return x*x*x

    def grad_v_vect(self, x):
        return np.multiply(x, np.multiply(x, x))

    def lapl_v_vect(self, x):
        return 3.0*np.multiply(x, x)


class Sixth(Langevin):
    def __init__(self):
        Langevin.__init__(self)

    def v(self, x):
        x_cube = x*x*x
        return x_cube*x_cube / 6

    def grad_v(self, x):
        x_sqd = x*x
        return x_sqd*x_sqd * x

    def grad_v_vect(self, x):
        x_sqd = np.multiply(x, x)
        return np.multiply(np.multiply(x, x_sqd), x_sqd)

    def lapl_v_vect(self, x):
        x_sqd = np.multiply(x, x)
        return 5.0 * np.multiply(x_sqd, x_sqd)


class Bistable(Langevin):
    def __init__(self):
        Langevin.__init__(self)

    def v(self, x):
        x_sqd = x * x
        return x_sqd * x_sqd / 4.0 - x_sqd / 2

    def grad_v(self, x):
        return x * x * x - x

    def grad_v_vect(self, x):
        return np.multiply(x, np.multiply(x, x)) - x

    def lapl_v_vect(self, x):
        return 3 * np.multiply(x, x) - 1.0

class MSLangevin(SDE):
    def __init__(self, eps):
        SDE.__init__(self)
        self.eps = eps

    def v(self, x):
        return float()

    def grad_v(self, x):
        return float()

    def lapl_v(self, x):
        return float()

    def p(self, x):
        return float()

    def grad_p(self, x):
        return float()

    def f(self, x, p):
        return -1.0 * p[0] * self.grad_v(x) - 1.0 / self.eps * self.grad_p(x / self.eps)

    def g(self, x, p):
        return np.sqrt(2.0 * p[-1])


class MSOrnUhl(MSLangevin):
    def __init__(self, eps):
        MSLangevin.__init__(self, eps)

    def v(self, x):
        return x*x / 2.0

    def grad_v(self, x):
        return x

    def lapl_v(self, x):
        return 1

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)


class MSQuartic(MSLangevin):
    def __init__(self, eps):
        MSLangevin.__init__(self, eps)

    def v(self, x):
        x_sqd = x*x
        return x_sqd*x_sqd / 4.0

    def grad_v(self, x):
        return x*x*x

    def lapl_v(self, x):
        return 3*x*x

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)


class MSSixth(MSLangevin):
    def __init__(self, eps):
        MSLangevin.__init__(self, eps)

    def v(self, x):
        x_sqd = x*x
        return x_sqd*x_sqd*x_sqd / 6.0

    def grad_v(self, x):
        x_sqd = x*x
        return x_sqd*x_sqd*x

    def lapl_v(self, x):
        x_sqd = x*x
        return 5.0*x_sqd*x_sqd

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

    def lapl_v(self, x):
        return 3.0*x*x - 1.0

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)
