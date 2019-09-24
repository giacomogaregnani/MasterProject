import numpy as np
import scipy.stats as stats


def moving_average(a, n):
    ret = np.zeros(np.size(a))
    for i in range(0, n):
        ret[i] = ret[i-1] + (a[i] - a[0]) / n
    for i in range(n, np.size(ret)):
        ret[i] = ret[i-1] + (a[i] - a[i-n]) / float(n)
    return ret


def get_density(x):
    kde = stats.gaussian_kde(x)
    x_sorted = np.sort(x)
    density = kde.evaluate(x_sorted)
    return x_sorted, density


class ParEst:
    def __init__(self, x, grad_v, h, is_vect=False):
        self.x = x
        self.grad_v = grad_v
        self.h = h
        self.x_diff = np.diff(x)
        self.n = np.size(x)
        self.is_vect = is_vect

    def drift(self, x_alt=np.zeros(0)):

        if np.size(x_alt) > 1:
            grad_v_eval = np.zeros(self.n - 1)
            if self.is_vect:
                grad_v_eval = self.grad_v(x_alt[0:-1])
            else:
                for i in range(0, self.n - 1):
                    grad_v_eval[i] = self.grad_v(x_alt[i])
        else:
            grad_v_eval = np.zeros(self.n - 1)
            if self.is_vect:
                grad_v_eval = self.grad_v(self.x[0:-1])
            else:
                for i in range(0, self.n - 1):
                    grad_v_eval[i] = self.grad_v(self.x[i])

        num = np.sum(np.multiply(self.x_diff, grad_v_eval))
        # print([num, num*40])
        den = np.sum(np.square(grad_v_eval))
        # print(den*self.h)
        return -1.0 * num / (self.h * den)

    def drift_subs(self, delta):
        drift_est = 0.0
        for i in range(0, delta):
            x_subs = self.x[i::delta]
            x_subs_diff = np.diff(x_subs)
            grad_v_eval = np.zeros(np.size(x_subs)-1)
            if self.is_vect:
                grad_v_eval = self.grad_v(x_subs[0:-1])
            else:
                for j in range(0, np.size(grad_v_eval)):
                    grad_v_eval[j] = self.grad_v(x_subs[j])
            num = np.sum(np.multiply(x_subs_diff, grad_v_eval))
            den = np.sum(np.square(grad_v_eval))
            drift_est_local = -1.0 * num / (delta * self.h * den)
            drift_est += drift_est_local / delta
        return drift_est

    def diffusion(self):
        return np.sum(np.square(self.x_diff)) / (2.0 * self.h * self.n)
