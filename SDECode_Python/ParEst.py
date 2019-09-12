import numpy as np


class ParEst:
    def __init__(self, x, v, h):
        self.x = x
        self.v = v
        self.h = h
        self.xdiff = np.diff(x)
        self.n = np.size(x)

    def drift(self):
        eval = np.zeros(self.n-1)
        for i in range(0, self.n-1):
            eval[i] = self.v(self.x[i])
        num = np.sum(np.multiply(self.xdiff, eval))
        denom = np.sum(np.square(eval))
        return -1.0*num / (self.h*denom)

    def diffusion(self):
        return np.sum(np.square(self.xdiff)) / (2.0 * self.h * self.n)
