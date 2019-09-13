import numpy as np


def moving_average(a, n=5):
    ret = np.convolve(a, np.ones(n) / n, mode='valid')
    return ret


class ParEst:
    def __init__(self, x, grad_v, h):
        self.x = x
        self.grad_v = grad_v
        self.h = h
        self.x_diff = np.diff(x)
        self.n = np.size(x)

    def drift(self):
        grad_v_eval = np.zeros(self.n - 1)
        for i in range(0, self.n - 1):
            grad_v_eval[i] = self.grad_v(self.x[i])
        num = np.sum(np.multiply(self.x_diff, grad_v_eval))
        den = np.sum(np.square(grad_v_eval))
        return -1.0 * num / (self.h * den)

    def diffusion(self):
        return np.sum(np.square(self.x_diff)) / (2.0 * self.h * self.n)
