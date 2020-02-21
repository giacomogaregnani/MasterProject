import numpy as np


class SDE:
    def __init__(self):
        pass

    def f(self, x, alpha):
        return float()

    def g(self, x, sigma):
        return float()


class Langevin(SDE):
    def __init__(self, n):
        SDE.__init__(self)
        self.n = n

    def v(self, x):
        return float()

    def grad_v(self, x):
        return np.empty(self.n)

    def grad_v_vect(self, x):
        pass

    def lapl_v_vect(self, x):
        pass

    def f(self, x, alpha):
        return -1.0 * np.dot(alpha, self.grad_v(x))

    def g(self, x, sigma):
        return np.sqrt(2.0 * sigma)


class OrnUhl(Langevin):
    def __init__(self):
        Langevin.__init__(self, 1)

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
        Langevin.__init__(self, 1)

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
        Langevin.__init__(self, 1)

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
        Langevin.__init__(self, 1)

    def v(self, x):
        x_sqd = x * x
        return x_sqd * x_sqd / 4.0 - x_sqd / 2

    def grad_v(self, x):
        return x * x * x - x

    def grad_v_vect(self, x):
        return np.multiply(x, np.multiply(x, x)) - x

    def lapl_v_vect(self, x):
        return 3 * np.multiply(x, x) - 1.0


class Bistable2(Langevin):
    def __init__(self):
        Langevin.__init__(self, 2)

    def v(self, x):
        ret_val = np.zeros(2)
        x_sqd = x * x
        ret_val[0] = x_sqd * x_sqd / 4.0
        ret_val[1] = -1.0 * x_sqd / 2
        return ret_val

    def grad_v(self, x):
        ret_val = np.zeros(2)
        ret_val[0] = x * x * x
        ret_val[1] = -1.0 * x
        return ret_val


class KLStyle(Langevin):
    def __init__(self, n):
        Langevin.__init__(self, n)
        self.n = n

    def v(self, x):
        ret_val = np.zeros(self.n)
        for i in range(0, self.n):
            ret_val[i] = pow(x, 2*i+2)
        return ret_val

    def grad_v(self, x):
        ret_val = np.zeros(self.n)
        for i in range(0, self.n):
            ret_val[i] = pow(x, 2*i+1)
        return ret_val


class Cheb(Langevin):
    def __init__(self, n):
        Langevin.__init__(self, n)
        self.n = n

    # TODO: Code this but it is useless
    def v(self, x):
        ret_val = np.zeros(self.n)
        return ret_val

    def grad_v(self, x):
        pol = np.zeros(self.n+1)
        d_pol = np.zeros(self.n+1)
        d_pol[0] = 0.0
        pol[0] = 1.0
        d_pol[1] = 1.0
        pol[1] = x
        for i in range(2, self.n+1):
            pol[i] = 2 * x * pol[i-1] - pol[i-2]
            d_pol[i] = 2 * pol[i-1] + 2 * x * d_pol[i-1] - d_pol[i-2]
        ret_val = d_pol[1:]
        return ret_val


class MSLangevin(SDE):
    def __init__(self, eps, n):
        SDE.__init__(self)
        self.eps = eps
        self.n = n

    def v(self, x):
        return float()

    def grad_v(self, x):
        return np.empty(self.n)

    def lapl_v(self, x):
        return float()

    def p(self, x):
        return float()

    def grad_p(self, x):
        return float()

    def f(self, x, alpha):
        return -1.0 * np.dot(alpha, self.grad_v(x)) - 1.0 / self.eps * self.grad_p(x / self.eps)

    def g(self, x, sigma):
        return np.sqrt(2.0 * sigma)


class MSOrnUhl(MSLangevin):
    def __init__(self, eps):
        MSLangevin.__init__(self, eps, 1)

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
        MSLangevin.__init__(self, eps, 1)

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
        MSLangevin.__init__(self, eps, 1)

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
        MSLangevin.__init__(self, eps, 1)

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


class MSBistable2(MSLangevin):
    def __init__(self, eps):
        MSLangevin.__init__(self, eps, 2)

    def v(self, x):
        ret_val = np.empty(2)
        x_sqd = x * x
        ret_val[0] = x_sqd * x_sqd / 4.0
        ret_val[1] = -1.0 * x_sqd / 2
        return ret_val

    def grad_v(self, x):
        ret_val = np.empty(2)
        ret_val[0] = x * x * x
        ret_val[1] = -1.0 * x
        return ret_val

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)


class MSKLStyle(MSLangevin):
    def __init__(self, eps, n):
        MSLangevin.__init__(self, eps, n)
        self.n = n

    def v(self, x):
        ret_val = np.zeros(self.n)
        for i in range(0, self.n):
            ret_val[i] = pow(x, 2*i + 2) / (2*i + 2)
        return ret_val

    def grad_v(self, x):
        ret_val = np.zeros(self.n)
        for i in range(0, self.n):
            ret_val[i] = pow(x, 2*i + 1)
        return ret_val

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)


class MSCheb(MSLangevin):
    def __init__(self, eps, n):
        MSLangevin.__init__(self, eps, n)
        self.n = n

    # TODO: Code this but it is useless
    def v(self, x):
        ret_val = np.zeros(self.n)
        return ret_val

    def grad_v(self, x):
        pol = np.zeros(self.n+1)
        d_pol = np.zeros(self.n+1)
        d_pol[0] = 0.0
        pol[0] = 1.0
        d_pol[1] = 1.0
        pol[1] = x
        for i in range(2, self.n+1):
            pol[i] = 2 * x * pol[i-1] - pol[i-2]
            d_pol[i] = 2 * pol[i-1] + 2 * x * d_pol[i-1] - d_pol[i-2]
        ret_val = d_pol[1:]
        return ret_val

    def p(self, x):
        return np.cos(x)

    def grad_p(self, x):
        return -1.0 * np.sin(x)
