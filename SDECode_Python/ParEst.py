import numpy as np
import scipy.stats as stats
from scipy.signal import convolve
from scipy.special import gamma


def filter_trajectory(a, delta, T, beta=1, type=0):
    n = np.size(a)
    h = T / (n-1)
    t_vec = np.arange(n) * h
    c_beta = beta / gamma(1.0/beta)

    if type == 0:
        k = c_beta * np.exp(-1.0 * (t_vec / delta) ** beta) / delta
    elif type == 1:
        k = c_beta * np.exp(-1.0 * (t_vec ** beta) / delta) / (delta ** (1.0 / beta))
    else:
        k = c_beta/delta - c_beta/(delta**(beta+1))*(t_vec**beta)
        k[k < 0] = 0

    # The two lines below are equivalent to
    # ret = np.zeros(n)
    # for i in range(0, n):
    #     for j in range(0, i):
    #         ret[i] += k[i-j] * a[j]
    # ret *= h

    ret = convolve(a, k, mode='full') * h
    ret = ret[0:n]

    return ret, k[0]


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

    def formula_verification(self, data=None, grad_p=None, eps=0, dw=None, sigma=0):
        grad_v_eval = self.grad_v(self.x[0:-1])
        num = self.h/eps * np.sum(np.multiply(grad_p(data[0:-1]/eps), grad_v_eval))
        den = self.h * np.sum(np.multiply(grad_v_eval, self.grad_v(data[0:-1])))
        num2 = -np.sqrt(2.0*sigma) * np.sum(np.multiply(grad_v_eval, dw))
        return num/den, num2/den

    def drift(self):
        grad_v_eval = self.grad_v(self.x[0:-1])
        num_summand = np.multiply(self.x_diff, grad_v_eval)
        num = np.sum(num_summand)
        den = np.sum(np.square(grad_v_eval))
        return -1.0 * num / (self.h * den)

    def drift_mixed(self, data=None):
        grad_v_eval = self.grad_v(self.x[0:-1])
        num_summand = np.multiply(np.diff(data), grad_v_eval)
        num = np.sum(num_summand)
        # den = np.sum(np.square(grad_v_eval))
        den = np.sum(np.multiply(grad_v_eval, self.grad_v(data[0:-1])))
        return -1.0 * num / (self.h * den)

    def drift_tilde(self, lapl, delta=1.0, sigma=0):
        if sigma == 0:
            sigma = self.diffusion() / delta
        grad_v_eval = self.grad_v(self.x[0:-1])
        den = np.sum(np.square(grad_v_eval)) * self.h
        num = sigma * np.sum(lapl(self.x[0:-1])) * self.h
        return num / den, sigma

    def diffusion(self):
        return np.sum(np.square(self.x_diff)) / (2.0 * self.h * self.n)

    def drift_bayesian(self, a_pr, var_pr):
        grad_v_eval = self.grad_v(self.x[0:-1])
        int_num = np.sum(np.multiply(self.x_diff, grad_v_eval))
        int_denom = self.h * np.sum(np.square(grad_v_eval))
        den = 1.0 + var_pr**2.0 * int_denom
        num = a_pr - var_pr**2.0 * int_num
        a_post = num / den
        var_post = 1.0 / (var_pr**(-2.0) + int_denom)
        return a_post, var_post
